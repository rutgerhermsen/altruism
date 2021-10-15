module collectives
  use parameters
  use shorthands
  use data_structures
  use fft
  use kernels
  use state
  use utils

  implicit none

  integer :: rep_c = 0            ! stores how many times the collectives have been analyzed

contains

  subroutine consider_collective_analysis(field, o_list, o_stats, t)
    type(fld), intent(inout) :: field(2)
    type(org), intent(inout) :: o_list(MAX_POP, 2)
    type(org_sts), intent(in) :: o_stats(2)
    integer, intent(in) :: t
    integer :: j = 0
    integer :: x = 0
    complex(DP) :: fft_K_smooth(N)
    complex(DP) :: K_smooth_c(N)
    logical :: exist
    real(DP) :: offs = real(RES_SMOOTH - N, kind = DP)/real(RES_SMOOTH - 1, kind = DP)
    real(DP) :: slope = real(N - 1, kind = DP)/real(RES_SMOOTH - 1, kind = DP)
    logical :: minima(N)
    logical :: tentative

    ! first update the smoothed field
    fft_K_smooth = cmplx(field(py)%occ, kind = DP)
    call cfft1f(N, 1, fft_K_smooth, N, wsave, lensav, work, lenwrk, ier)
    K_smooth_c = fft_K_smooth*fft_smooth_kernel*cmplx(N, kind = DP)
    call cfft1b(N, 1, K_smooth_c, N, wsave, lensav, work, lenwrk, ier)


    if (t .ne. 0) then
      field(py)%K_smooth = (1d0 - TIME_SMOOTHING)*real(K_smooth_c, kind = DP) + TIME_SMOOTHING*field(py_N)%K_smooth
    else
      field(py)%K_smooth = real(K_smooth_c, kind = DP)
    end if

    if (ANALYZE_COLLECTIVES .and. (t >=  rep_c*ANALYSIS_INTERVAL)) then
      rep_c = rep_c + 1
      if (OUTPUT_SMOOTHED_FIELD) then
        inquire(file="smoothed_field.txt", exist=exist)
        if (exist) then
          open(unit=13, file="smoothed_field.txt", form="formatted", action="write", position="append")
        else
          open(unit=13, file="smoothed_field.txt", form="formatted", action="write")
        end if

        write(13,*) (real(K_smooth_c(nint(offs + slope*real(x, kind = DP))), kind = DP), x = 1,RES_SMOOTH)

        close(unit=13)
      end if
      ! first, find borders and characterize collectives tentatively
      tentative = .true.
      call find_collective_borders(field(py), c_list, minima, tentative)
      call characterize_collectives(tentative)
      ! then, quality-check the borders
      call check_borders(minima)
      ! now, characterize collectives definitively
      tentative = .false.
      call find_collective_borders(field(py), c_list, minima, tentative)
      call characterize_collectives(tentative)

      if (OUTPUT_COLL) then
        call output_collectives(t)
      end if
      call find_w_list_MLS1()
      call find_w_list_coll_MLS1()
      call find_w_list_coll_MLS2(t)
      call calculate_coll_selection(t)

      if (pyC_N .eq. 1) then
        pyC_N = 2
        pyC = 1
      else
        pyC_N = 1
        pyC = 2
      end if
      o_stats_prev = o_stats(py)
      !mean_p_prev = o_stats(py)%mean_p ! store the mean of this state
      !n_prev = o_stats(py)%n
      o_list(1:o_stats(py)%n, py)%a = (/(j, j=1,o_stats(py)%n, 1)/)
      o_list_prev = o_list(:, py)
    end if
  end subroutine consider_collective_analysis


  subroutine find_collective_borders(field, c_list, minima, tentative)
    type(fld), intent(in) :: field
    type(coll), intent(inout) :: c_list(MAX_COLL, 2)
    logical, intent(inout) :: minima(N)
    logical, intent(in) :: tentative
    real(DP) :: deriv(N)
    logical :: before(N), after(N)
    integer :: ones(N) = 1
    logical :: endpoints(N), spoints(N), shifted_spoints(N)
    logical :: lowDensity(N)
    integer :: i, count_c = 1


    deriv = (cshift(field%K_smooth, -1) - cshift(field%K_smooth, 1))

    if (tentative) then
      ! find list of minima
      lowDensity = (field%K_smooth < K_THRES_WEAK)

      before = (cshift(field%K_smooth, 1) > field%K_smooth)
      after =  (cshift(field%K_smooth, -1) > field%K_smooth)

      minima = ((before .and. after) .and. lowDensity)
    end if

    c_stats(pyC_N)%n = sum(ones, minima)

    if (c_stats(pyC_N)%n .eq. 0) then
      c_stats(pyC_N)%n = 1
    end if

    c_list(1:c_stats(pyC_N)%n, pyC_N)%l = 1
    c_list(1:c_stats(pyC_N)%n, pyC_N)%r = 1
    c_list(1:c_stats(pyC_N)%n, pyC_N)%minr = 1
    if (c_stats(pyC_N)%n > 1) then
      endpoints = minima
      spoints = minima .and. (deriv > 0d0)
      shifted_spoints = cshift(spoints, -1)
      where (spoints)
        endpoints = .false.
      end where
      where (shifted_spoints)
        endpoints = .true.
      end where
      count_c = 1
      do i = 1,N
        if (endpoints(i)) then
          c_list(count_c, pyC_N)%r = i
          if (shifted_spoints(i)) then
            c_list(count_c, pyC_N)%minr = mm(i + 1)
          else
            c_list(count_c, pyC_N)%minr = i
          end if
          count_c = count_c + 1
          c_list(modulo(count_c - 1, c_stats(pyC_N)%n) + 1, pyC_N)%l = mm(i + 1)
        end if
      end do
    else
      c_list(1, pyC_N)%l = 1
      c_list(1, pyC_N)%r = N
      c_list(1, pyC_N)%minr = N ! this case is rather ill-defined. Should not happen where we care about it.
    end if

  end subroutine find_collective_borders


  subroutine characterize_collectives(tentative)
    logical, intent(in) :: tentative
    integer :: myc
    integer :: jmax = 0
    integer :: maxi = 0
    integer :: o = 0
    integer :: i = 0
    integer :: j = 0

    c_list(1:c_stats(pyC_N)%n, pyC_N)%p = 0d0
    c_list(1:c_stats(pyC_N)%n, pyC_N)%x = 0d0
    c_list(1:c_stats(pyC_N)%n, pyC_N)%n = 0
    origins(:,:) = 0
    do o = 1,o_stats(py)%n
      myc = assign_to_collective(o)
      origins(myc, o_list(o, py)%c) = origins(myc, o_list(o, py)%c) + 1
      if (.not. tentative) then
        o_list(o, py)%c = myc
      end if
      c_list(myc, pyC_N)%p = c_list(myc, pyC_N)%p + o_list(o, py)%p
      if (myc .ne. 1) then ! periodic boundary causes problems when averaging position of collective 1
        c_list(myc, pyC_N)%x = c_list(myc, pyC_N)%x + real(o_list(o, py)%x, kind = DP)
      else
        if (o_list(o, py)%x <= c_list(1, pyC_N)%r) then
          c_list(myc, pyC_N)%x = c_list(myc, pyC_N)%x + real(o_list(o, py)%x + N, kind = DP)
        else
          c_list(myc, pyC_N)%x = c_list(myc, pyC_N)%x + real(o_list(o, py)%x, kind = DP)
        end if
      end if
      c_list(myc, pyC_N)%n = c_list(myc, pyC_N)%n + 1
    end do

    c_list(1:c_stats(pyC_N)%n, pyC_N)%p = c_list(1:c_stats(pyC_N)%n, pyC_N)%p / &
    & real(c_list(1:c_stats(pyC_N)%n, pyC_N)%n, kind = DP)
    c_list(1:c_stats(pyC_N)%n, pyC_N)%x = c_list(1:c_stats(pyC_N)%n, pyC_N)%x / &
    & real(c_list(1:c_stats(pyC_N)%n, pyC_N)%n, kind = DP)
    c_list(1, pyC_N)%x = modulo(c_list(1, pyC_N)%x - 5d-1, N_REAL) + 5d-1

    ! find, for each collective, the parent collective
    ! we assume that this is the collective that contributed most
    ! in terms of the organisms
    c_list(1:c_stats(pyC)%n, pyC)%offspring_MLS2 = 0
    c_list(1:c_stats(pyC_N)%n, pyC_N)%parent= 0
    do i = 1, c_stats(pyC_N)%n
      maxi = 0
      jmax = 1
      do j = 1, c_stats(pyC)%n
        if (origins(i,j) > maxi ) then ! I accept a bias in case of a draw...
          maxi = origins(i,j)
          jmax = j
        end if
      end do
      c_list(jmax, pyC)%offspring_MLS2 = c_list(jmax, pyC)%offspring_MLS2 + 1
      c_list(i, pyC_N)%parent = jmax
    end do

  end subroutine characterize_collectives


  subroutine check_borders(minima)
    logical, intent(inout) :: minima(N)
    integer :: i
    logical :: strictdens(N)
    integer :: countrejected

    countrejected = 0

    strictdens = (field(py)%K_smooth < K_THRES_STRICT)

    if ( (.not. strictdens(c_list(c_stats(pyC_N)%n, pyC_N)%minr)) .and. &
    & ( c_list(c_stats(pyC_N)%n, pyC_N)%parent == c_list(1, pyC_N)%parent)) &
    & then
    minima(c_list(c_stats(pyC_N)%n, pyC_N)%minr) = .false.
    countrejected = countrejected + 1
  end if

  do i = 1,(c_stats(pyC_N)%n -1)

    if ( (.not. strictdens(c_list(i, pyC_N)%minr)) .and. &
    & ( c_list(i, pyC_N)%parent == c_list(i + 1, pyC_N)%parent)) &
    & then
    minima(c_list(i, pyC_N)%minr) = .false.
    countrejected = countrejected + 1
  end if
end do

end subroutine


subroutine output_collectives(t)
  integer, intent(in) :: t
  logical :: exist
  integer :: j

  inquire(file="coll_list_time.txt", exist=exist)
  if (exist) then
    open(unit=13, file="coll_list_time.txt", form="formatted", action="write", position="append")
  else
    open(unit=13, file="coll_list_time.txt", form="formatted", action="write")
    write(13, "(g0)") "timestep, time, left_border, right_border, center_of_mass, mean_altruism"
  end if
  do j = 1,c_stats(pyC_N)%n
    write(13, "(*(g0))" &!, advance = "no"
    & ) t, ', ', real(t, kind = DP) * DT, ', ', &
    &c_list(j, pyC_N)%l, ', ', c_list(j, pyC_N)%r, ', ', &
    &c_list(j, pyC_N)%x, ', ', c_list(j, pyC_N)%p
  end do
  close(unit=13)

end subroutine


subroutine find_w_list_coll_MLS1()
  integer :: parent_c
  real(DP) :: mean_offspring_MLS1

  if (rep_c > 1) then
    ! find, for each parent collective
    ! the number of organisms it produced in the new collectives
    do parent_c = 1, c_stats(pyC)%n
      c_list(parent_c, pyC)%offspring_MLS1 = sum(origins(1:c_stats(pyC)%n, parent_c))
    end do
    mean_offspring_MLS1 = real(o_stats(py)%n, kind = DP)/real(o_stats_prev%n, kind = DP)

    c_list(1:c_stats(pyC)%n , pyC)%w_MLS1 = (real(c_list(1:c_stats(pyC)%n, pyC)%offspring_MLS1, kind = DP) &
    & / real(c_list(1:c_stats(pyC)%n, pyC)%n, kind = DP))/mean_offspring_MLS1
  end if
end subroutine find_w_list_coll_MLS1

subroutine find_w_list_MLS1()
  integer :: o
  real(DP) :: mean_offspring_MLS1

  if (rep_c > 1) then
    ! find, for each organism in the parent population
    ! the number of organisms it produced in the new collectives
    o_list_prev(1:o_stats_prev%n)%offspring_MLS1 = 0
    do o = 1, o_stats(py)%n
      o_list_prev(o_list(o, py)%a)%offspring_MLS1 = o_list_prev(o_list(o, py)%a)%offspring_MLS1 + 1
    end do
    mean_offspring_MLS1 = real(o_stats(py)%n, kind = DP)/real(o_stats_prev%n, kind = DP)
    o_list_prev(1:o_stats_prev%n)%w_MLS1 = real(o_list_prev(1:o_stats_prev%n)%offspring_MLS1, kind = DP) / mean_offspring_MLS1

  end if
end subroutine find_w_list_MLS1

subroutine find_w_list_coll_MLS2(t)
  integer, intent(in) :: t
  integer :: line_id = 0
  logical :: exist = .false.
  real(DP) :: dis = 0d0
  real(DP) :: xo, xp = 0d0
  real(DP) :: deltat = 0d0
  real(DP) :: mean_offspring_MLS2 = 0d0
  integer :: i = 0
  integer :: j = 0
  integer :: c_id_counter = 0                     !

  if (rep_c < 2) then
    do i = 1, c_stats(pyC_N)%n
      c_list(i, pyC_N)%c_id = c_id_counter
      c_id_counter = c_id_counter + 1
    end do
  else
    if (c_stats(pyC_N)%n < 2) then
      print *, "Number of collectives is smaller than 2 at time ", t
    end if
    mean_offspring_MLS2 = real(c_stats(pyC_N)%n, kind = DP) / &
    & real(c_stats(pyC)%n, kind = DP)
    c_list(1:c_stats(pyC)%n , pyC)%w_MLS2 = real(c_list(1:c_stats(pyC)%n, pyC)%offspring_MLS2,kind = DP) / &
    & mean_offspring_MLS2

    ! identify whether the collections were just born and give them an id
    do i = 1, c_stats(pyC_N)%n
      if (c_list(c_list(i, pyC_N)%parent, pyC)%offspring_MLS2 > 1) then
        c_list(i, pyC_N)%c_id = c_id_counter
        c_id_counter = c_id_counter + 1
      else
        c_list(i, pyC_N)%c_id = c_list(c_list(i, pyC_N)%parent, pyC)%c_id
      end if
    end do

    if (OUTPUT_COLL_DEMOGRAPHICS) then
      inquire(file="coll_demographics.txt", exist=exist)
      if (exist) then
        open(unit=13, file="coll_demographics.txt", form="formatted", action="write", position="append")
      else
        open(unit=13, file="coll_demographics.txt", form="formatted", action="write")
        write(13, "(g0)") "timestep, time, center_of_mass, type"
      end if
      do i = 1, c_stats(pyC)%n
        if (c_list(i, pyC)%offspring_MLS2 > 1) then
          write(13, "(*(g0))") real(t, kind = DP) - ANALYSIS_INTERVAL_REAL/2d0, ', ', &
          & (real(t, kind = DP) - ANALYSIS_INTERVAL_REAL/2d0) * DT, ', ', &
          & c_list(i, pyC)%x, ', ', "b"
        else
          if (c_list(i, pyC)%offspring_MLS2 < 1) then
            write(13, "(*(g0))") real(t, kind = DP) - ANALYSIS_INTERVAL_REAL/2d0, ', ', &
            & (real(t, kind = DP) - ANALYSIS_INTERVAL_REAL/2d0) * DT, ', ', &
            & c_list(i, pyC)%x, ', ', "d"
          end if
        end if
      end do
      close(unit=13)
    end if

    if (OUTPUT_COLL_LIST_TIME_LINES) then
      inquire(file="coll_list_time_lines.txt", exist=exist)
      if (exist) then
        open(unit=13, file="coll_list_time_lines.txt", form="formatted", action="write", position="append")
      else
        open(unit=13, file="coll_list_time_lines.txt", form="formatted", action="write")
        write(13, "(g0)") "# specifies line segments from parent to offspring collectives in space-time."
        write(13, "(g0)") "# Start and end points of line segments are given by rows with the same 'line id'"
        write(13, "(g0)") "# Given the periodic boundaries, the shortes route between two positions"
        write(13, "(g0)") "# can cross the boundary; in that case two lines are drawn"
        write(13, "(g0)") "timestep, time, position, phenotype offspring, line id, id of offspring collective"
      end if
      do j = 1,c_stats(pyC_N)%n
        xp = c_list(c_list(j, pyC_N)%parent, pyC)%x
        xo = c_list(j, pyC_N)%x

        dis = real(abs(xp - xo), kind = DP)
        ! shortest route does not cross the boundary:
        if (dis < (N_REAL/2d0)) then
          write(13, "(*(g0))") real(t, kind = DP) - ANALYSIS_INTERVAL_REAL, ', ', &
          & (real(t, kind = DP) - ANALYSIS_INTERVAL_REAL)*DT, ', ', &
          & xp, ', ', c_list(j, pyC_N)%p, ', ', line_id, ', ', c_list(j, pyC_N)%c_id
          write(13, "(*(g0))") real(t, kind = DP), ', ', &
          & real(t, kind = DP)*DT, ', ', &
          & xo, ', ', c_list(j, pyC_N)%p, ', ', &
          & line_id, ', ', c_list(j, pyC_N)%c_id
          line_id = line_id + 1
        ! shortest route does cross the boundary:
        else
          if ( xp > xo) then ! parent position is higher than offspring position; line crosses top border:
            deltat = (N_REAL + 5d-1 - xp) * &
            & ANALYSIS_INTERVAL_REAL / &
            & (xo - xp + N_REAL)
            write(13, "(*(g0))") real(t, kind = DP) - ANALYSIS_INTERVAL_REAL, ', ', &
            & (real(t, kind = DP) - ANALYSIS_INTERVAL_REAL)*DT, ', ', &
            & xp, ', ', c_list(j, pyC_N)%p, ', ', line_id, ', ', c_list(j, pyC_N)%c_id
            write(13, "(*(g0))") real(t, kind = DP) - ANALYSIS_INTERVAL_REAL + deltat, ', ', &
            & (real(t, kind = DP) - ANALYSIS_INTERVAL_REAL + deltat)*DT, ', ', &
            & N_REAL + 5d-1, ', ', c_list(j, pyC_N)%p, ', ', &
            & line_id, ', ', c_list(j, pyC_N)%c_id

            line_id = line_id + 1

            write(13, "(*(g0))") real(t, kind = DP) - ANALYSIS_INTERVAL_REAL + deltat, ', ', &
            & (real(t, kind = DP) - ANALYSIS_INTERVAL_REAL + deltat)*DT, ', ', &
            & 5d-1, ', ', c_list(j, pyC_N)%p, ', ', line_id, ', ', c_list(j, pyC_N)%c_id
            write(13, "(*(g0))") real(t, kind = DP), ', ', &
            & real(t, kind = DP)*DT, ', ', &
            & xo, ', ', c_list(j, pyC_N)%p, ', ', &
            & line_id, ', ', c_list(j, pyC_N)%c_id
            line_id = line_id + 1
          else ! parent position is lower than offspring position; line crosses top border:
            deltat = (xp - 5d-1) * &
            & ANALYSIS_INTERVAL_REAL / &
            & (N_REAL + xp - xo)
            write(13, "(*(g0))") real(t, kind = DP) - ANALYSIS_INTERVAL_REAL, ', ', &
            & (real(t, kind = DP) - ANALYSIS_INTERVAL_REAL)*DT, ', ', &
            & xp, ', ', c_list(j, pyC_N)%p, ', ', line_id, ', ', c_list(j, pyC_N)%c_id
            write(13, "(*(g0))") real(t, kind = DP) - ANALYSIS_INTERVAL_REAL + deltat, ', ', &
            & (real(t, kind = DP) - ANALYSIS_INTERVAL_REAL + deltat)*DT, ', ', &
            & 5d-1, ', ', c_list(j, pyC_N)%p, ', ', &
            & line_id, ', ', c_list(j, pyC_N)%c_id
            line_id = line_id + 1

            write(13, "(*(g0))") real(t, kind = DP) - ANALYSIS_INTERVAL_REAL + deltat, ', ', &
            & (real(t, kind = DP) - ANALYSIS_INTERVAL_REAL + deltat)*DT, ', ', &
            & N_REAl + 5d-1, ', ', c_list(j, pyC_N)%p, ', ', line_id, ', ', c_list(j, pyC_N)%c_id
            write(13, "(*(g0))") real(t, kind = DP), ', ', &
            & real(t, kind = DP)*DT, ', ', &
            & xo, ', ', c_list(j, pyC_N)%p, ', ', &
            & line_id, ', ', c_list(j, pyC_N)%c_id
            line_id = line_id + 1
          end if
        end if
      end do
      close(unit=13)
    end if

  end if
end subroutine find_w_list_coll_MLS2


integer function assign_to_collective(oA) result(guess)
  integer, intent(in) :: oA
  logical :: found

  found = .false.
  guess = min(o_list(oA, py)%c, c_stats(pyC_N)%n)

  do while (.not. found)
    if (guess .eq. 1) then
      if ( &
      & (c_list(guess, pyC_N)%r >= o_list(oA, py)%x) .or. &
      & ((c_list(guess, pyC_N)%l <= o_list(oA, py)%x) &
      &     .and. (c_list(guess, pyC_N)%r < c_list(guess, pyc_N)%l) &
      &    )) then
      found = .true.
    else
      guess = guess + 1
    end if
  else
    if ((c_list(guess, pyC_N)%l <= o_list(oA, py)%x) &
    & .and. (c_list(guess, pyC_N)%r >= o_list(oA, py)%x) ) then
    found = .true.
  else
    if (c_list(guess, pyC_N)%r < o_list(oA, py)%x) then
      guess = modulo(guess, c_stats(pyC_N)%n) + 1
    else
      guess = guess - 1
    end if
  end if
end if
end do

end function assign_to_collective


subroutine calculate_coll_selection(t)
  integer, intent(in) :: t
  real(DP) :: S_MLS1 = 0d0
  real(DP) :: transmission_MLS1 = 0d0
  real(DP) :: S_among_MLS1 = 0d0
  real(DP) :: S_within_MLS1 = 0d0
  real(DP) :: S_c_MLS2 = 0d0
  real(DP) :: Delta_mean_c = 0d0
  real(DP) :: meanAfter_c = 0d0
  real(DP) :: meanBefore_c = 0d0
  logical :: exist = .false.

  if (rep_c > 1) then
    !  first: selection MLS2
    S_c_MLS2 = simple_covariance( &
    & c_list(1:c_stats(pyC)%n, pyC)%w_MLS2 - c_stats(pyC)%mean_w_MLS2, &
    & c_list(1:c_stats(pyC)%n, pyC)%p - o_stats(py)%mean_p, &
    & c_stats(pyC)%n)

    meanAfter_c = pSum( &
    & c_list(1:c_stats(pyC_N)%n, pyC_N)%p) / real(c_stats(pyC_N)%n, kind = DP)
    meanBefore_c = pSum( &
    & c_list(1:c_stats(pyC)%n, pyC)%p) / real(c_stats(pyC)%n, kind = DP)


    Delta_mean_c = meanAfter_c - meanBefore_c

    ! then: selection MLS1, Price Method

    ! selection:
    S_MLS1 = simple_covariance( &
    & o_list_prev(1:o_stats_prev%n)%w_MLS1 - 1d0, &
    & o_list_prev(1:o_stats_prev%n)%p - o_stats_prev%mean_p, &
    & o_stats_prev%n)
    transmission_MLS1 = o_stats(py)%mean_p - o_stats_prev%mean_p - S_MLS1


    ! selection among collectives:
    S_among_MLS1 = pSum( &
    & (c_list(1:c_stats(pyC)%n, pyC)%w_MLS1 - c_stats(pyC)%mean_w_MLS1 ) * &
    & (c_list(1:c_stats(pyC)%n, pyC)%p - o_stats_prev%mean_p) * &
    & real(c_list(1:c_stats(pyC)%n, pyC)%n, kind = DP) &
    & ) &
    & / real(o_stats_prev%n, kind = DP)

    S_within_MLS1 = S_MLS1 - S_among_MLS1

    !        S_among_MLS1 = pSum( &
    !            & c_list(1:n_c(pyC) , pyC)%w_MLS1 * &
    !            & c_list(1:n_c(pyC), pyC)%p * &
    !            & real(c_list(1:n_c(pyC), pyC)%n, kind = DP) &
    !            & ) &
    !           & / real(n_prev, kind = DP) - &
    !           & (pSum( &
    !            & c_list(1:n_c(pyC), pyC)%p * &
    !            & real(c_list(1:n_c(pyC), pyC)%n, kind = DP) &
    !            & ) &
    !           & / real(n_prev, kind = DP)) * &
    !            & (pSum( &
    !            & c_list(1:n_c(pyC), pyC)%w_MLS1 * &
    !            & real(c_list(1:n_c(pyC), pyC)%n, kind = DP) &
    !            & ) &
    !           & / real(n_prev, kind = DP))

    inquire(file="selection_coll.txt", exist=exist)
    if (exist) then
      open(unit=13, file="selection_coll.txt", form="formatted", action="write", position="append")
    else
      open(unit=13, file="selection_coll.txt", form="formatted", action="write")
      write(13, "(g0)") "timestep, time, &
      & meanOfOrganisms, changeInMeanOrganisms, selectionOrganisms, transmissionOrganisms, SamongMLS1, SwithinMLS1, &
      & meanOfCollectives, changeInMeanCollectives, selectionCollectives, TransmissionCollectives"
    end if
    write(13, "(*(g0))") t, ', ', real(t, kind = DP)*DT, ', ',&
    & o_stats(py)%mean_p, ', ', o_stats(py)%mean_p - o_stats_prev%mean_p, ', ', S_MLS1, ', ', transmission_MLS1, ', ', &
    & S_among_MLS1, ', ', S_within_MLS1, ', ', &
    & meanAfter_c, ', ', Delta_mean_c, ', ', S_c_MLS2, ', ', Delta_mean_c - S_c_MLS2
    close(unit=13)
  end if

end subroutine calculate_coll_selection



end module collectives
