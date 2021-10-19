module parameters
  implicit none

  integer,  parameter ::  DP = selected_real_kind(12, 60)

  !-------------------------------------------------------------------
  ! >>>>>> PARAMETERS
  !
  ! >> DISCRETIZATION
  !
  real(DP), parameter :: RESOLUTION = 1d1 ! number of grid cells in social range
  real(DP), parameter :: DT = 8d-2 ! time step in units if 1/(basal growth rate)

  !-------------------------------------------------------------------
  !
  ! >> CHOICES
  !
  logical, parameter :: REALIZED_FITNESS = .false.
  logical, parameter :: STEP_FUNCTION_KERNEL = .true.

  logical, parameter :: OUTPUT_CONTR_SB_FIELD = .true. ! output local contribution to S_BELOW
  real(DP), parameter :: RANGE_CONTR_SB = 1d0 ! band width of corresponding kernel

  logical, parameter :: OUTPUT_OCC_FIELD = .true.  ! periodically output occupancy?
  logical, parameter :: OUTPUT_P_FIELD = .true. ! periodically output p_field?
  logical, parameter :: OUTPUT_ALTR_FIELD = .true.   ! periodically output altruism field?
  logical, parameter :: OUTPUT_O_LIST = .true. ! periodically output the full o_list?
  logical, parameter :: OUTPUT_G_OF_R = .true. ! periodically output the g(r)
  logical, parameter :: OUTPUT_ORDER = .true. ! periodically output the order parameter

  logical, parameter :: OUTPUT_HISTOGRAMS = .true.   ! periodically output a histogram of phenotypes?
  integer, parameter :: HIST_NR_BINS = 400 ! number of bins of the histogram
  real(DP), parameter :: HIST_MAX = 5d-2              ! right-hand border of last bin

  logical, parameter :: OUTPUT_COUNT_GLOBAL_LIST = .true.
  logical, parameter :: OUTPUT_ALL_S_CURVES = .true. ! write all S-curves to file
  logical, parameter :: OUTPUT_SELECTION_LIST = .true.

  logical, parameter :: WELL_MIXED = .false. ! use well-mixed conditions.


  !-------------------------------------------------------------------
  !
  ! >> ORGANISM AND ECOLOGY
  !
  real(DP), parameter :: RATE_GROWTH = 5d0 ! basal reproduction rate
  real(DP), parameter :: RATE_DEATH = 1d0 ! death rate, equals 1 per definition of the time unit.
  real(DP), parameter :: P_R = RATE_GROWTH*DT ! basal reproduction probability per time step; 1 per definition
  real(DP), parameter :: P_D = RATE_DEATH*DT ! death probability per time step
  real(DP), parameter :: P_MUT = 1d-3 ! mutation probability per reproduction
  ! the effect size of mutations is exponentially distributed
  real(DP), parameter :: MEAN_SIZE_MUT = 5d-3 ! mean size mutation

  real(DP), parameter :: K_DIFF = 4d-2 ! diffusion conctant of organisms

  ! organisms interact both socally and through resource competition.
  ! both interactions have a range
  real(DP), parameter :: RANGE_SOCIAL = 1d0 ! social interaction range
  real(DP), parameter :: RANGE_COMPETITION =  4d0 ! competition interaction range
  real(DP), parameter :: RANGE_DISPERSION = sqrt(2d0 * K_DIFF * DT) ! sd of the dispersion

  real(DP), parameter :: K_COMPETITION = 4d1  !  carrying capacity density of resource competition; length unit is RANGE_SOCIAL

  ! costs and benefits
  real(DP), parameter :: COST = 1d0 ! cost per unit of altruism; set to 1 per definition
  ! the benefits of the altruism saturate
  real(DP), parameter :: MAX_BENEFIT = 5d0 ! upper limit of benefit
  real(DP), parameter :: BENEFIT = 1d0  ! benefit per unit of altruism at the linear (low-altruism) regime

  ! >> ARENA
  integer, parameter :: N = 2**9!2**10! linear size of the simulated frame. preferably power of 2

  ! >> SIMULATION SETTINGS

  integer, parameter :: seed(33) = 1

  ! last timestep
  integer, parameter :: T_MAX = 101
  ! how frequently to report the state of the system, in timesteps
  integer, parameter :: REPORT_INTERVAL = 50
  integer, parameter :: MEAN_INTERVAL = 1000

  ! at the end of the simulation, we calculate S_local and S_interlocal for various scales sigma
  ! we do this for many time steps and average the result
  ! we start NR_AV time steps before the end of the simulation
  ! and calculate the S_local-curve every SC_IN time steps
  ! making sure to always end at time MAX_T - 1.
  ! (We cannot calculate S_local curves at MAX_T because we need
  ! realized fitnesses for that.)
  integer, parameter :: NR_AV = 3 !
  integer, parameter :: SC_IN = 1 ! interval between scurves
  integer, parameter :: NR_SC = (NR_AV-1)/SC_IN + 1  ! number of S-curves to be computed

  integer, parameter :: NR_SCALES = 82     ! for how many scales sigma do we compute S_below?

  ! we use periodic boundary conditions
  ! this means that the kernel used for calculations of S_below
  ! must take a certain number number of periodic images into account
  ! in both dimensions
  ! we look NR_FIELDS up, down, left, and right
  integer, parameter :: NR_FIELDS = 7

  ! initial conditions
  ! initial population size in units of the estimated steady-state density SS
  real(DP), parameter :: INIT_POP = 1d0

contains

end module parameters
