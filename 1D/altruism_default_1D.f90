program altruism_default_1D

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
  use collectives     ! ... routines for calculating multilevel selection (MLS1 and MLS2)

  implicit none

  ! bookkeeping
  integer :: t = 0    ! time in computational timesteps
  integer :: rep = 0  ! stores how many times the state of the system has been reported
  integer :: o        ! loop variable organism

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
  call initial_conditions(o_list, o_stats, field)

  ! >>> GO!
  write(*, '(a)') "# Starting the simulation."
  timesteps: &
  do while (t < T_MAX .and. o_stats(py)%n > 0)
    ! wipe variables that wil be used to store
    ! the new state to be constructed
    call cleanup()

    ! update density etc. of current ("old") state
    call update_quantities()
    
    ! perform collective-level analysis, at certain times
    call consider_collective_analysis(field, o_list, o_stats, t)

    demographics: &
    ! for all organisms in old state ...
    do o = 1,o_stats(py)%n

      ! ...  determine their fate ...
      mother: &
      if (.not. consider_death()) then
        ! keep in the next time step ...
        call keep(o, o_list, o_stats(py_N))
        ! ... but disperse
        call disperse(o_list(:, py_N), o_stats(py_N)%n)
        ! keep track of occupancy of gridpoints etc.
        call update_occ_p_field(o_stats(py_N)%n)
      end if mother

      ! ... and determine the fate of their children.
      reproduction: &
      if (consider_reproduction(o)) then
        child: &
        if (.not. consider_death()) then
          call keep(o, o_list, o_stats(py_N))
          ! children can obtain mutations
          call consider_mutation(o, o_stats(py_N)%n)
          call disperse(o_list(:, py_N), o_stats(py_N)%n)
          call update_occ_p_field(o_stats(py_N)%n)
        end if child
      end if reproduction
    end do demographics

    ! is the total population size out of bound?
    call check_n()

    !!!! IMPORTANT!
    !!!! Even though we just constructed the new state,
    !!!! we will now calculate S-curves and output for the *old* state.
    !!!! We do this this late in the process because
    !!!! we can only know the "realized" fitness
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

  ! If the simulation made it all the way to T_MAX,
  if (t == T_MAX) then
    ! calculate statistics
    call compute_means()
    ! output the local and interlocal selection at the chosen scales
    if (OUTPUT_ALL_S_CURVES) then
      call output_S_curves()
    end if
  end if

  call cleanup_fft()
  call cleanup_state()

contains

  ! Subroutine updates quantities of current state
  subroutine update_quantities()
    ! updates the local densities that determine the competition
    call update_K_dens()
    ! updates the local level of altruism
    call update_K_altr()
    ! updates the growth rates
    call update_growth_rates(o_list, o_stats, field, py)
    ! update the mean phenotype
    o_stats(py)%mean_p = &
    & pSum(o_list(1:o_stats(py)%n, py)%p) / real(o_stats(py)%n, kind = DP)
  end subroutine update_quantities

  ! S_curves are calculated only at certain times
  subroutine consider_S_curve()
    if (&
    & NR_AV > 0 .and. t >=  (T_MAX - NR_AV) .and. &
    & (modulo(T_MAX - t - 1, SC_IN) == 0) &
    & ) then
      print *, "# Start calculating an S-curve at time", t
      call update_w_field()
      call compute_S_curve(field(py), o_stats(py))
    end if
  end subroutine consider_S_curve

  ! wipe memory of the variables that will be used to construct the new state
  subroutine cleanup()
    o_stats(py_N)%n = 0
    field(py_N)%occ = 0
    field(py_N)%p = 0d0
    field(py_N)%p_e = 0d0
  end subroutine cleanup

  ! we assumed throughout that the population size does not exceed MAX_POP
  ! we check this here
  subroutine check_n()
    if (o_stats(py_N)%n > MAX_POP) then
      print *, "Fatal problem!!!"
      print *, "n too large!!!"
      print *, "MAX_POP:", MAX_POP
      print *, "n(py_N):", o_stats(py_N)%n
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
  ! in o_list pyA
  subroutine update_growth_rates(o_listA, o_statsA, fieldA, pyA)
    type(org), intent(inout) :: o_listA(MAX_POP, 2)
    type(org_sts), intent(inout) :: o_statsA(2)
    type(fld), intent(inout) :: fieldA(2)
    integer, intent(in) :: pyA

    integer i

    do i = 1, o_statsA(pyA)%n
      o_listA(i, pyA)%gr = growthrate(i, pyA, o_listA, fieldA)
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
    complex(DP) :: fft_occ_field(N)

    fft_occ_field = cmplx(field(py)%occ, kind = DP)
    call cfft1f(N, 1, fft_occ_field, N, wsave, lensav, work, lenwrk, ier)
    field(py)%K_dens = fft_occ_field*fft_kernel*cmplx(N, kind = DP)
    call cfft1b(N, 1, field(py)%K_dens, N, wsave, lensav, work, lenwrk, ier)
  end subroutine update_K_dens

  ! update the level of altruism
  subroutine update_K_altr()
    complex(DP) :: fft_p_e_field(N)

    fft_p_e_field = cmplx(field(py)%p_e, kind = DP)
    call cfft1f(N, 1, fft_p_e_field, N, wsave, lensav, work, lenwrk, ier)
    field(py)%K_altr = fft_p_e_field*fft_kernelAltr*cmplx(N, kind = DP)
    call cfft1b(N, 1, field(py)%K_altr, N, wsave, lensav, work, lenwrk, ier)
  end subroutine update_K_altr

  ! after each REPORT_INTERVAL the state of the system is reported
  ! and certain properties written to file
  subroutine consider_progress_report(t, rep, nA)
    integer, intent(inout) :: nA
    integer, intent(in) :: t
    integer, intent(inout) :: rep
    complex(DP) :: autocorr(N)
    real(DP) :: g_of_r(0:HALF_N)


    if ( t >= rep*REPORT_INTERVAL) then
      rep = rep + 1
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

      ! the actual output happens here
      call save_snapshot(rep, g_of_r)
    end if

  end subroutine consider_progress_report

  subroutine disperse(o_listA, oN)
    type(org), intent(inout) :: o_listA(MAX_POP)
    integer, intent(in) :: oN
    integer :: dx
    real(DP) :: rndx

    dx = nint(rnor()*RANGE_DISPERSION*RESOLUTION)

    if (WELL_MIXED) then
      call random_number(rndx)
      o_listA(oN)%x = ceiling(rndx*N_REAL)
    else
      o_listA(oN)%x = mm(o_listA(oN)%x + dx)
    end if

  end subroutine

  subroutine update_occ_p_field(oN)
    integer, intent(in) :: oN

    field(py_N)%occ(o_list(oN, py_N)%x) = &
    & field(py_N)%occ(o_list(oN, py_N)%x) + 1
    field(py_N)%p(o_list(oN, py_N)%x) = &
    & field(py_N)%p(o_list(oN, py_N)%x) + (o_list(oN, py_N)%p)
    field(py_N)%p_e(o_list(oN, py_N)%x) = &
    & field(py_N)%p_e(o_list(oN, py_N)%x) + (o_list(oN, py_N)%p_e)
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
      o_list(oN, py_N)%p = o_list(oN, py_N)%p + delta_p
      if (o_list(oN, py_N)%p > 0d0) then
        o_list(oN, py_N)%p_e = o_list(oN, py_N)%p
      else
        o_list(oN, py_N)%p_e = 0d0
      end if
      o_list(o, py)%tot_Delta_p = o_list(o, py)%tot_Delta_p + delta_p
    end if

  end subroutine consider_mutation

  subroutine save_snapshot(rep, g_of_r)
    integer, intent(in) :: rep
    real(DP), intent(in) :: g_of_r(0:HALF_N)
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
      write(13, "(g0)") "position, grid-point, phenotype, relative fitness"
      do i=1,o_stats(py)%n
        write(13, "(*(g0))") real(o_list(i, py)%x, kind = DP)/RESOLUTION, ', ', &
        & o_list(i, py)%x, ', ', &
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
      open(unit=13, file="occ_field"//fileRep//".txt", form="formatted", &
      & action="write", status="replace")
      write(13, "(*(g0))") "# occupancy at timestep ", t, ", time ", real(t, kind = DP)*DT
      write(13, "(g0)") "position, grid-point, occupancy"
      do i=1,N
        write(13, "(*(g0))") real(i, kind = DP)/RESOLUTION, &
        & ', ', i , ', ', field(py)%occ(i)
      end do
      close(unit=13)
    end if

    ! if desired, outputs the total value of the phenotypes at each position
    if (OUTPUT_P_FIELD) then
      open(unit=13, file="p_field"//fileRep//".txt", &
      & form="formatted", action="write", status="replace")
      write(13, "(*(g0))") "# total phenotype at timestep ", t, ", time ", real(t, kind = DP)*DT
      write(13, "(g0)") "position, grid-point, total phenotype"
      do i=1,N
        write(13, "(*(g0))") real(i, kind = DP)/RESOLUTION, ', ', i, ', ', &
        & real(field(py)%p(i), kind = 4)
      end do
      close(unit=13)
    end if

    ! if desired, outputs the level of altruism experienced at each position
    if (OUTPUT_ALTR_FIELD) then
      open(unit=13, file="altr"//fileRep//".txt", form="formatted", &
      & action="write", status="replace")
      write(13, "(*(g0))") "# level of altruism experienced at timestep ", t, ", time ", real(t, kind = DP)*DT
      write(13, "(g0)") "position, grid-point, level of altruism"
      do i=1,N
        write(13, "(*(g0))") real(i, kind = DP)/RESOLUTION, ', ', i, ', ', &
        & real(field(py)%K_altr(i), kind = 4)
      end do
      close(unit=13)
    end if

    ! if desired, outputs the contribution of each position to local selection
    ! (the LSD multiplied with local density)
    ! based on scale/band width RANGE_CONTR_SB
    if (OUTPUT_CONTR_SB_FIELD) then
      open(unit=13, file="contrSlocal"//fileRep//".txt", form="formatted", action="write", status="replace")
      write(13, "(*(g0))") "# contribution of each position to local selection at timestep ", t, ", time ", real(t, kind = DP)*DT
      write(13, "(g0)") "position, grid-point, contribution to local selection"
      do i=1,N
        write(13, "(*(g0))") real(i, kind = DP)/RESOLUTION, ', ', i, ', ', &
        & real(field(py)%contrib_Slocal(i), kind = 4)
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
        open(12, file="histograms.txt", status="old", &
        & position="append", action="write")
      else
        open(12, file="histograms.txt", status="new", action="write")
      end if
      write(12, *) histo
      close(unit=12)
    end if

    if (OUTPUT_G_OF_R) then
      open(unit=13, file="g_of_r"//fileRep//".txt", &
      & form="formatted", action="write", status="replace")
      write(13, "(*(g0))") "# pair-correlation function at timestep ", t, ", time ", real(t, kind = DP)*DT
      write(13, "(g0)") "position, g(r)"
      do j=0, HALF_N
        if (distances(j) > 0) then
          write(13, "(*(g0))") real(j, kind = DP)/RESOLUTION, ', ', &
          & g_of_r(j)/real(distances(j), kind = DP)
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
      bin = floor(o_list(o, py)%p_e*real(HIST_NR_BINS, kind = DP)/HIST_MAX) + 1
      if (bin > HIST_NR_BINS) then
        histo(HIST_NR_BINS) = histo(HIST_NR_BINS) + 1
      else
        histo(bin) = histo(bin) + 1
      end if
    end do
  end subroutine bin_p

  ! calculates the growth rate of individual oA in the old field
  real(DP) function growthrate(oA, pyA, o_listA, fieldA) result (gr)
    integer, intent(in) :: oA
    integer, intent(in) :: pyA
    type(org), intent(inout) :: o_listA(MAX_POP, 2)
    type(fld), intent(inout) :: fieldA(2)
    real(DP) :: ben

    ben = real(fieldA(pyA)%K_altr(o_listA(oA, pyA)%x), kind = DP)

    gr = P_R * min( max( &
    & (1d0 - COST*o_listA(oA, pyA)%p_e + MAX_BENEFIT*ben / (ben + MAX_BENEFIT/BENEFIT)) * &
    & (1d0 - real(fieldA(pyA)%K_dens(o_listA(oA, pyA)%x), kind = DP)/K_COMPETITION), 0d0), 1d0)
  end function growthrate

  ! installs individuals at uniformly random positions
  ! sets phenotypes to 0
  subroutine initial_conditions(o_listA, o_statsA, fieldA)
    type(org), intent(inout) :: o_listA(MAX_POP, 2)
    type(org_sts), intent(inout) :: o_statsA(2)
    type(fld), intent(inout) :: fieldA(2)

    real(DP) :: rnd = 0d0
    integer :: pyL = 2 ! at time 0, py = 2

    write(*, '(a)', advance = "no") "# Initializing..."
    ! set initial population size
    o_statsA(pyL)%n = INIT_TOTPOP

    do o = 1, o_statsA(pyL)%n
      call random_number(rnd)
      ! choose a random position
      o_listA(o, pyL)%x = ceiling(rnd*N_REAL)
      ! set phenotypes to 0
      o_listA(o, pyL)%p = 0d0
      ! adjust the occupancy field and the total phenotype field
      fieldA(pyL)%occ(o_listA(o, pyL)%x) = &
      & fieldA(pyL)%occ(o_listA(o, pyL)%x) + 1
      fieldA(pyL)%p(o_listA(o, pyL)%x) = &
      & fieldA(pyL)%p(o_listA(o, pyL)%x) + o_listA(o, pyL)%p
    end do

    print *, "Done."
  end subroutine initial_conditions

  ! promote organism o to the new state
  subroutine keep(o, o_listA, o_statsA)
    integer, intent(in) :: o
    type(org), intent(inout) :: o_listA(MAX_POP, 2)
    type(org_sts), intent(inout) :: o_statsA
    type(org) :: offspring

    ! construct the offspring from the parent
    offspring = o_listA(o,  py)
    ! reset the properties that are not inherited
    offspring%offspring = 0
    offspring%tot_Delta_p = 0d0
    offspring%w = 0d0
    offspring%w_MLS1 = 0d0

    ! add the offspring to the new field
    call add_org(offspring, o_listA(:, py_N), o_statsA)

    ! increase the number of offspring for
    ! organism o of the old state
    o_listA(o, py)%offspring = o_list(o, py)%offspring + 1
  end subroutine keep

  ! add organism org to the list o_listA
  ! and adjust the corresponding total population size
  subroutine add_org(orgA, o_listA, o_statsA)
    type(org), intent(in) :: orgA
    type(org), intent(inout) :: o_listA(MAX_POP)
    type(org_sts), intent(inout) :: o_statsA

    ! increase population size of new state
    o_statsA%n = o_statsA%n + 1
    ! copy the organism to the last position of the new state
    o_listA(o_statsA%n) = orgA
  end subroutine add_org

end program altruism_default_1D
