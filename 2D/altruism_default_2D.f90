program altruism_default_2D

  ! modules providing...
  use ziggurat        ! ... cheap Gaussian random numbers
  use pairwise        ! ... recursive pairwise summation
  use parameters      ! ... parameters
  use shorthands      ! ... convenient shorthands
  use data_structures ! ... data structures
  use utils           ! ... useful utility functions
  use fft             ! ... routines preparing memory for FFT
  use kernels         ! ... routines to construct (FFT of) kernel functions
  use g_of_r          ! ... routines to calculate radial distribution function
  use state           ! ... structures storing representations of the state of the simulation
  use price           ! ... routines for calculating terms of the Price equation
  use multiscale      ! ... routines for calculating multiscale selection

  implicit none

  ! bookkeeping
  integer :: t = 0     ! time
  integer :: rep = 0   ! stores how many times the state of the system has been reported
  integer :: o = 0     ! loop variable organism

  !preparation
  call prepare_fft()              ! prepare (memory for) FFT
  call set_shorthands()           ! define shorthands used throughout
  call prepare_state()            ! allocate memory for list of organisms
  call assemble_kernels()         ! pre-assemble FFTs of kernels

  ! prepare calculation of radial distribution function
  if (OUTPUT_G_OF_R) then
    call prepare_distances()
  end if

  ! set random seed (for reproducibility)
  call random_seed(put = seed)

  ! set initial conditions
  call initial_conditions()

  ! >>> GO!
  write(*, '(a)') "# Starting the simulation."
  timesteps: &
  do while (t < T_MAX .and. o_stats(py)%n > 0)
    ! wipe variables that wil be used to store
    ! the new state to be constructed
    call cleanup()

    ! update density etc. of current ("old") state
    call update_quantities()

    demographics: &
    ! for all organisms in the old state
    do o = 1,o_stats(py)%n

      ! ...  determine their fate ...
      mother: &
      if (.not. consider_death()) then
        ! keep in the next time step ...
        call keep(o)
        ! ... but disperse
        call disperse(o_stats(py_N)%n)
        ! keep track of occupancy of gridpoints etc.
        call update_occ_p_field(o_stats(py_N)%n)
      end if mother

      ! ... and determine the fate of their children.
      reproduction: &
      if (consider_reproduction(o)) then
        child: &
        if (.not. consider_death()) then
          call keep(o)
          ! children can obtain mutations
          call consider_mutation(o, o_stats(py_N)%n)
          call disperse(o_stats(py_N)%n) !
          call update_occ_p_field(o_stats(py_N)%n)
        end if child
      end if reproduction
    end do demographics

    ! is the total population out of bound?
    call check_n()

    !!!! IMPORTANT!
    !!!! Even though we just constructed the new state,
    !!!! we will now calculate S-curves and output for the *old* state.
    !!!! We do this this late in the process because
    !!!! we can only know the "realize" fitness
    !!!! after constructing the new state.

    ! calculate directional selection, purifying selection, drift, transmission
    ! using the Price Equation.
    call calculate_price(t)

    ! calculate local and interlocal selection for chosen scales,
    ! at certain times only
    call consider_S_curve()

    ! report progress & calculate radial distribution function, at certain times.
    call consider_progress_report(t, rep, o_stats(py)%n)

    ! update time
    t = t + 1

    ! swap data arrays for odd and even times.
    ! the "old" state is pointed to by
    py = modulo(t + 1, 2) + 1
    ! the "new" state is pointed to by
    py_N = modulo(t, 2) + 1

  end do timesteps

  ! If the simulation made it all the way to T_MAX
  if (t == T_MAX) then
    ! calculate statistics
    call compute_means()
    ! output the local and interlocal selection at the chosen scales
    if (OUTPUT_ALL_S_CURVES) then
      call output_S_curves()
    end if
    if (OUTPUT_MEAN_S_CURVE) then
      call output_S_curve_mean()
    end if
  end if

  call cleanup_fft()
  call cleanup_state()

contains

  subroutine update_quantities()
    ! updates the local densities that determine the competition
    call update_K_dens()
    ! updates the local altruism
    call update_K_altr()
    ! update the growth-rate list
    call update_growth_rates()
    ! update the mean phenotype
    o_stats(py)%mean_p = &
    & pSum(o_list(1:o_stats(py)%n, py)%p) / real(o_stats(py)%n, kind = DP)
  end subroutine update_quantities

  ! S_curves are calculated only at certain times
  subroutine consider_S_curve()
    integer :: i
    logical :: exis
    character(len=10) :: fileRep
    
    if (&
    & NR_AV > 0 .and. t >=  (T_MAX - NR_AV) &
    & .and. (modulo(T_MAX - t - 1, SC_IN) == 0)&
    & ) then
    print *, "# Start calculating an S-curve at time", t
    call update_w_field()
    call compute_S_curve(field(py), o_stats(py))
    if (OUTPUT_SELECTION_LIST) then
        inquire(file="selection_Hilje.txt", exist=exis)
        if (exis) then
            open(12, file="selection_Hilje.txt", status="old", position="append", action="write")
        else
            open(12, file="selection_Hilje.txt", status="new", action="write")
            write(12, "(*(g0))", advance = "no") &
                & "time", TAB, "Stot"
            do i = 1, (NR_SCALES + 1)
                write (fileRep, "(F10.4)") S_curve_vals(nr_SC_C, i, 1)
                write(12, "(*(g0))", advance = "no") &
                & TAB, "Slocal_"//adjustl(trim(fileRep)), TAB,  &
                & "Sinterlocal_"//adjustl(trim(fileRep))
            end do
            write(12, *)
        end if
        write(12, "(E23.16, a, E23.16)", advance = "no") &
                & t*DT, TAB, S
        do i = 1, (NR_SCALES + 1)
            write(12, "(a, E23.16, a, E23.16)", advance = "no") &
                & TAB, S_curve_vals(nr_SC_C, i, 3), TAB,  S_curve_vals(nr_SC_C, i, 2)
        end do
        write(12, *)
        close(unit=12)
    end if
  end if
end subroutine consider_S_curve

subroutine cleanup()
  o_stats(py_N)%n = 0
  field(py_N)%occ = 0
  field(py_N)%p = 0
end subroutine cleanup

! wipe memory of the variables that will be used to construct the new state
subroutine check_n()
  if (o_stats(py_N)%n > MAX_POP) then
    print *, "Fatal problem!!!"
    print *, "n too large!!!"
    print *, "MAX_POP:", MAX_POP
    print *, "o_stats(py_N)%n:", o_stats(py_N)%n
    call exit(0)
  end if
end subroutine check_n

! does the individual die?
logical function consider_death() result (died)
  real(DP) :: rnd

  call random_number(rnd)
  died = (rnd < P_D)
end function consider_death

! calculate the growth rate for all individuals
subroutine update_growth_rates()
  integer i

  do i = 1, o_stats(py)%n
    o_list(i, py)%gr = growthrate(i)
  end do
end subroutine

! consider whether the individual oA reproduces
logical function consider_reproduction(oA) result (reproduces)
  integer, intent(in) :: oA
  real(DP) :: rnd

  call random_number(rnd)
  reproduces = (rnd < o_list(oA, py)%gr)
end function consider_reproduction

! update the local density that determines the competition
subroutine update_K_dens()
  complex(DP) :: fft_occ_field(N,N)

  fft_occ_field = cmplx(field(py)%occ, kind = DP)
  call cfft2f(N, N, N, fft_occ_field, wsave, lensav, work, lenwrk, ier)
  field(py)%K_dens = fft_occ_field*fft_kernel*cmplx(N**2, kind = DP)
  call cfft2b(N, N, N, field(py)%K_dens, wsave, lensav, work, lenwrk, ier)
end subroutine update_K_dens

! update the level of altruism
subroutine update_K_altr()
  complex(DP) :: fft_p_field(N,N)

  fft_p_field = cmplx(field(py)%p, kind = DP)
  call cfft2f(N, N, N, fft_p_field, wsave, lensav, work, lenwrk, ier)
  field(py)%K_altr = fft_p_field*fft_kernelAltr*cmplx(N**2, kind = DP)
  call cfft2b(N, N, N, field(py)%K_altr, wsave, lensav, work, lenwrk, ier)
end subroutine update_K_altr

! after each REPORT_INTERVAL the state of the system is reported
! and certain properties written to file
subroutine consider_progress_report(t, rep, nA)
  integer, intent(inout) :: nA
  integer, intent(in) :: t
  integer, intent(inout) :: rep
  complex(DP) :: autocorr(N,N)
  real(DP) :: g_of_r(0:HALF_N**2)
  complex(DP) :: order(0:HALF_N**2)

  if ( t >= rep*REPORT_INTERVAL) then
    rep = rep + 1;
    print *, &
      & "Progress: Timestep = ", t, "; Population size: ", nA, "; Mean level altruism = ", o_stats(py)%mean_p

    ! If desired, the local selection estimate is calculated
    if (OUTPUT_CONTR_SB_FIELD) then
      call update_w_field()
      call update_w_p_field()
      call compute_contr_Slocal(field(py))
    end if

    ! if desired, the pair-distribution function g(r) is calculated
    if (OUTPUT_G_OF_R) then
      call compute_autocorrelation(autocorr, field(py)%occ, o_stats(py)%n)
      call compute_g_of_r(autocorr, g_of_r)
    end if

    ! if desired, the hexagonal order parameter is calculated
    if (OUTPUT_ORDER) then
      call compute_autocorrelation(autocorr, field(py)%occ, o_stats(py)%n)
      call compute_order(autocorr, order)
    end if

    ! the actual output happens here
    call save_snapshot(rep, g_of_r, order)
  end if

end subroutine consider_progress_report

subroutine compute_order(autocorr, order)
  complex(DP), intent(in) :: autocorr(N,N)
  complex(DP), intent(out):: order(0:HALF_N**2)
  integer :: col, row
  integer :: dcol, drow
  integer :: dis_sq
  real(DP) :: angle

  order = (0d0, 0d0)
  do col = 1,N
    do row = 1, N
      if (row - 1 < N - row + 1) then
        drow = row -1
        dcol = col -1
      else
        drow = row - 1 - N
        dcol = col - 1 - N
      end if
      dis_sq = drow**2 + dcol**2

      if (dis_sq < HALF_N**2 + 1) then
        angle = atan2(real(drow, kind = DP), real(dcol, kind = DP))
        order(dis_sq) = order(dis_sq) + autocorr(row, col)*exp((0d0, 6d0)*angle)
      end if
    end do
  end do

end subroutine

subroutine disperse(oN)
  integer, intent(in) :: oN
  integer :: dx
  integer :: dy
  real(DP) :: rndx
  real(DP) :: rndy

  dy = nint(rnor()*RANGE_DISPERSION*RESOLUTION)
  dx = nint(rnor()*RANGE_DISPERSION*RESOLUTION)

  if (WELL_MIXED) then
    call random_number(rndy)
    call random_number(rndx)
    o_list(oN, py_N)%y = floor(rndy*N_REAL) + 1
    o_list(oN, py_N)%x = floor(rndx*N_REAL) + 1
  else
    o_list(oN, py_N)%y = mm(o_list(oN ,py_N)%y + dy)
    o_list(oN, py_N)%x = mm(o_list(oN, py_N)%x + dx)
  end if

end subroutine

subroutine update_occ_p_field(oN)
  integer, intent(in) :: oN

  field(py_N)%occ(o_list(oN, py_N)%y, o_list(oN, py_N)%x) = &
  & field(py_N)%occ(o_list(oN, py_N)%y, o_list(oN, py_N)%x) + 1
  field(py_N)%p(o_list(oN, py_N)%y, o_list(oN, py_N)%x) = &
  & field(py_N)%p(o_list(oN, py_N)%y, o_list(oN, py_N)%x) + (o_list(oN, py_N)%p)
end subroutine update_occ_p_field

subroutine consider_mutation(o, oN)
  integer, intent(in) :: o
  integer, intent(in) :: oN
  real(DP) :: delta_p
  real(DP) :: rnd1, rnd2

  call random_number(rnd1)

  if (rnd1 < P_MUT) then
    call random_number(rnd2)
    if (rnd2 < 5d-1) then
      delta_p = -MEAN_SIZE_MUT*random_exponential()
    else
      delta_p = MEAN_SIZE_MUT*random_exponential()
    end if

    ! if the phenotype becomes negative, set it to zero
    if (o_list(oN, py_N)%p + delta_p > 0d0) then
      o_list(oN, py_N)%p = o_list(oN, py_N)%p + delta_p
    else
      delta_p = -o_list(oN, py_N)%p
      o_list(oN, py_N)%p = 0d0
    end if

    o_list(o, py)%tot_Delta_p = o_list(o, py)%tot_Delta_p + delta_p
  end if

end subroutine consider_mutation

subroutine save_snapshot(rep, g_of_r, order)
  integer, intent(in) :: rep
  real(DP), intent(in) :: g_of_r(0:HALF_N**2)
  complex(DP), intent(in) :: order(0:HALF_N**2)
  integer :: i,j
  logical :: exis
  integer :: histo(HIST_NR_BINS)
  character(len=4) :: fileRep

  write (fileRep, "(I4.4)") rep

  ! always export certain statistics to the file "selection.txt"
  inquire(file="selection.txt", exist=exis)
  if (exis) then
    open(12, file="selection.txt", status="old", position="append", action="write")
  else
    open(12, file="selection.txt", status="new", action="write")
    write(12, "(g0)") "timestep, time, population size, mean phenotype, &
    & mean phenotype next step, &
    & mean absolute fitness, &
    & selection differential S, rolling mean of S, transmission T, &
    & rolling mean of T, drift D, purifying selection, &
    & rolling mean of purifying selection, &
    & purifying selection (alternative definition), &
    & selection differential S (cumulative), &
    & transmission (cumulative), drift (cumulative), &
    & purifying selection (cumulative)"
  end if
  write(12, "(*(g0))") &
  & t, ', ', real(t, kind = DP)*DT, ', ', &
  & o_stats(py)%n, ', ', &
  & o_stats(py)%mean_p, ', ', &
  & pSum(o_list(1:o_stats(py_N)%n, py_N)%p) &
  & / real(o_stats(py_N)%n, kind = DP), ', ', &
  & pSum(o_list(1:o_stats(py)%n, py)%W_abs) &
  & / real(o_stats(py)%n, kind = DP), ', ', &
  & S, ', ', &
  & pSum(S_r_mean)/real(min(MEAN_INTERVAL, t + 1), kind = DP), ', ', &
  & transmission, ', ', &
  & pSum(transmission_r_mean)/real(min(MEAN_INTERVAL, t + 1), kind = DP), ', ', &
  & drift, ', ', &
  & S_pur, ', ', &
  & pSum(S_pur_r_mean)/real(min(MEAN_INTERVAL, t + 1), kind = DP), ', ', &
  & S_pur_alt, ', ', S_cum, ', ', transmission_cum, ', ', &
  & drift_cum,  ', ', S_pur_cum
  close(unit=12)

  ! if desired, outputs a list with properties of each organism
  if (OUTPUT_O_LIST) then
    open(unit=13, file="o_list"//fileRep//".txt", form="formatted", action="write", status="replace")
    write(13, "(*(g0))") "# organisms at timestep ", t, ", time ", real(t, kind = DP)*DT
    write(13, "(g0)") "position x, position y, grid-point x, grid-point y, phenotype, relative fitness"
    do i=1,o_stats(py)%n
      write(13, "(*(g0))") real(o_list(i, py)%x, kind = DP)/RESOLUTION, ', ', &
      & real(o_list(i, py)%y, kind = DP)/RESOLUTION, ', ',&
      & o_list(i, py)%x, ', ', &
      & o_list(i, py)%y, ', ', &
      & o_list(i, py)%p, ', ', &
      & o_list(i, py)%w
    end do
    close(unit=13)
  end if

  ! if desired, outputs basic population-level statistics intended for plotting
  if (OUTPUT_COUNT_GLOBAL_LIST) then
    inquire(file="countglobal.txt", exist=exis)
    if (exis) then
      open(12, file="countglobal.txt", status="old", position="append", action="write")
    else
      open(12, file="countglobal.txt", status="new", action="write")
      write(12, "(*(g0))") "time", TAB, "population_size", TAB, "mean_phenotype"
    end if
    write(12, "(*(g0))") real(t, kind = DP)*DT, TAB, o_stats(py)%n, TAB, o_stats(py)%mean_p
    close(unit=12)
  end if

  ! if desired, outputs the occupancy (number of individuals) on each position
  if (OUTPUT_OCC_FIELD) then
    open(unit=13, file="occ_field"//fileRep//".txt", form="formatted", action="write", status="replace")
    write(13, "(*(g0))") "# occupancy at timestep ", t, ", time ", real(t, kind = DP)*DT
    write(13, "(g0)") "# matrix format."
    do j=1,N
      write(13,*) (field(py)%occ(i, j), i= 1,N)
    end do
    close(unit=13)
  end if

  ! if desired, outputs the total value of the phenotypes at each position
  if (OUTPUT_P_FIELD) then
    open(unit=13, file="p_field"//fileRep//".txt", form="formatted", action="write", status="replace")
    write(13, "(*(g0))") "# total phenotype at timestep ", t, ", time ", real(t, kind = DP)*DT
    write(13, "(g0)") "# matrix format."
    do j=1,N
      write(13,*) (real(field(py)%p(i, j), kind = 4), i= 1,N)
    end do
    close(unit=13)
  end if

  ! if desired, outputs the level of altruism experienced at each position
  if (OUTPUT_ALTR_FIELD) then
    open(unit=13, file="altr"//fileRep//".txt", form="formatted", action="write", status="replace")
    write(13, "(*(g0))") "# level of altruism experienced at timestep ", t, ", time ", real(t, kind = DP)*DT
    write(13, "(g0)") "# matrix form"
    do j=1,N
      write(13,*) (real(field(py)%K_altr(i, j), kind = 4), i= 1,N)
    end do
    close(unit=13)
  end if

  ! if desired, outputs the contribution of each position to local selection
  ! (the LSD multiplied with local density)
  ! based on scale/band width RANGE_CONTR_SB
  if (OUTPUT_CONTR_SB_FIELD) then
    open(unit=13, file="contrSlocal"//fileRep//".txt", form="formatted", action="write", status="replace")
    write(13, "(*(g0))") "# contribution of each position to local selection at timestep ", t, ", time ", real(t, kind = DP)*DT
    write(13, "(g0)") "# matrix form"
    do j=1,N
      write(13,*) (real(field(py)%contrib_Slocal(i, j), kind = 4), i= 1,N)
    end do
    close(unit=13)
  end if

  ! if desired, outputs a frequency table of the phenotypes
  ! each table is a row;
  ! each new table is added as a new row
  if (OUTPUT_HISTOGRAMS) then
    call bin_p(histo)

    inquire(file="histograms.txt", exist=exis)
    if (exis) then
      open(12, file="histograms.txt", status="old", position="append", action="write")
    else
      open(12, file="histograms.txt", status="new", action="write")
    end if
    write(12, *) histo
    close(unit=12)
  end if

  if (OUTPUT_G_OF_R) then
    open(unit=13, file="g_of_r"//fileRep//".txt", form="formatted", action="write", status="replace")
    write(13, "(*(g0))") "# pair-correlation function at timestep ", t, ", time ", real(t, kind = DP)*DT
    write(13, "(g0)") "position, g(r)"
    do j=0, HALF_N**2
      if (distances(j) > 0) then
        write(13, "(*(g0))") sqrt(real(j, kind = DP))/RESOLUTION, ', ', &
        & g_of_r(j)/real(distances(j), kind = DP)
      end if
    end do
    close(unit=13)
  end if

  if (OUTPUT_ORDER) then
    open(unit=13, file="order"//fileRep//".txt", form="formatted", action="write", status="replace")
    write(13, "(*(g0))") "# hexagonal order parameer at timestep ", t, ", time ", real(t, kind = DP)*DT
    write(13, "(g0)") "position, order_parameter"
    do j=0, HALF_N**2
      if (distances(j) > 0) then
        write(13, "(*(g0))") sqrt(real(j, kind = DP))/RESOLUTION, ', ', &
        & abs(order(j))/real(distances(j), kind = DP)
      end if
    end do
    close(unit=13)
  end if

end subroutine save_snapshot

! creates a frequency table for the phenotypes.
! negative phenotypes are treated as 0
! phenotypes above HIST_MAX are put in the last bin
subroutine bin_p(histo)
  integer, intent(inout) :: histo(HIST_NR_BINS)
  integer :: o
  integer :: bin

  histo = 0
  do o = 1, o_stats(py)%n
    bin = floor(o_list(o, py)%p*real(HIST_NR_BINS, kind = DP)/HIST_MAX) + 1
    if (bin > HIST_NR_BINS) then
      histo(HIST_NR_BINS) = histo(HIST_NR_BINS) + 1
    else
      histo(bin) = histo(bin) + 1
    end if
  end do
end subroutine bin_p

! calculates the growth rate of individual oA in the old field
real(DP) function growthrate(oA) result (gr)
  integer, intent(in) :: oA
  real(DP) :: ben

  ben = real(field(py)%K_altr(o_list(oA, py)%y, o_list(oA, py)%x), kind = DP)

  gr = P_R * min( max( &
  & (1d0 - COST*o_list(oA, py)%p + MAX_BENEFIT*ben / (ben + MAX_BENEFIT/BENEFIT)) * &
  & (1d0 - real(field(py)%K_dens(o_list(oA, py)%y, o_list(oA, py)%x), kind = DP)/K_COMPETITION), 0d0), 1d0)
end function growthrate

! installs individuals at uniformly random positions
! sets phenotypes to 0
subroutine initial_conditions()
  real(DP) :: rnd = 0d0

  write(*, '(a)', advance = "no") "# Initializing..."
  ! set initial population size
  o_stats(py)%n = INIT_TOTPOP

  do o = 1, o_stats(py)%n
    call random_number(rnd)
    o_list(o, py)%x = ceiling(rnd*N_REAL)
    call random_number(rnd)
    o_list(o, py)%y = ceiling(rnd*N_REAL)
    o_list(o, py)%p = 0d0

    field(py)%occ(o_list(o, py)%y, o_list(o, py)%x) = field(py)%occ(o_list(o, py)%y, o_list(o, py)%x) + 1
    field(py)%p(o_list(o, py)%y, o_list(o, py)%x) = field(py)%p(o_list(o, py)%y, o_list(o, py)%x) + o_list(o, py)%p
  end do

  print *, "Done."
end subroutine initial_conditions


subroutine keep(o)
  integer, intent(in) :: o

  o_stats(py_N)%n = o_stats(py_N)%n + 1
  o_list(o_stats(py_N)%n, py_N) = o_list(o, py)
  o_list(o, py)%offspring = o_list(o, py)%offspring + 1

  o_list(o_stats(py_N)%n, py_N)%offspring = 0
  o_list(o_stats(py_N)%n, py_N)%tot_Delta_p = 0d0
  o_list(o_stats(py_N)%n, py_N)%w = 0d0
  o_list(o_stats(py_N)%n, py_N)%W_abs = 0d0

end subroutine keep

end program altruism_default_2D
