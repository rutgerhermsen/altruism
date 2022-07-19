module parameters
  implicit none

  integer, parameter  ::  DP = selected_real_kind(12, 60)


  !-------------------------------------------------------------------
  ! >>>>>> PARAMETERS
  !
  ! >> DISCRETIZATION
  !
  real(DP), parameter :: RESOLUTION = 8d1 ! number of grid cells in social range
  real(DP), parameter :: DT = 8d-2 ! time step in units if 1/(death rate)

  !-------------------------------------------------------------------
  !
  ! >> CHOICES
  !
  logical, parameter :: REALIZED_FITNESS = .false.
  logical, parameter :: STEP_FUNCTION_KERNEL = .true.

  logical, parameter :: OUTPUT_CONTR_SB_FIELD = .false. ! output local contribution to S_local
  real(DP), parameter :: RANGE_CONTR_SB = 1d0 ! band width of corresponding kernel

  logical, parameter :: ANALYZE_COLLECTIVES = .true. ! identify and analyze the dynamics of the collectives

  logical, parameter :: OUTPUT_OCC_FIELD = .false.  ! periodically output occupancy?
  logical, parameter :: OUTPUT_P_FIELD = .false. ! periodically output p_field?
  logical, parameter :: OUTPUT_ALTR_FIELD = .false.   ! periodically output altruism field?
  logical, parameter :: OUTPUT_O_LIST = .false. ! periodically output the full o_list?
  logical, parameter :: OUTPUT_G_OF_R = .false. ! periodically output the g(r)
  logical, parameter :: OUTPUT_SMOOTHED_FIELD = .true.
  logical, parameter :: OUTPUT_COLL = .true. ! output the properties of collectives
  logical, parameter :: OUTPUT_COLL_DEMOGRAPHICS = .true. ! output the demographic events at the level of collectives
  logical, parameter :: OUTPUT_COLL_LIST_TIME_LINES = .true. ! output line segments that trace the collectives in space-time
  logical, parameter :: OUTPUT_HAMILTON = .true. ! calculate cost, benefit, relatedness according to Hamilton's rule

  logical, parameter :: OUTPUT_HISTOGRAMS = .false.   ! periodically output a histogram of phenotypes?
  integer, parameter :: HIST_NR_BINS = 400 ! number of bins of the histogram
  real(DP), parameter :: HIST_MAX = 5d-2              ! right-hand border of last bin

  logical, parameter :: OUTPUT_COUNT_GLOBAL_LIST = .false.
  logical, parameter :: OUTPUT_ALL_S_CURVES = .false. ! write all S-curves to file
  logical, parameter :: OUTPUT_MEAN_S_CURVE = .true.  ! write mean of S-curves to file

  logical, parameter :: WELL_MIXED = .false. ! use well-mixed conditions.


  !-------------------------------------------------------------------
  !
  ! >> ORGANISM AND ECOLOGY
  !
  real(DP), parameter :: RATE_GROWTH = 5d0 ! basal reproduction rate
  real(DP), parameter :: RATE_DEATH = 1d0 ! death rate, equals 1 per definition of the time unit.
  real(DP), parameter :: P_R = RATE_GROWTH*DT ! basal reproduction probability per time step
  real(DP), parameter :: P_D = RATE_DEATH*DT ! death probability per time step
  real(DP), parameter :: P_MUT = 5d-4 ! mutation probability per reproduction
  ! the effect size of mutations is exponentially distributed
  real(DP), parameter :: MEAN_SIZE_MUT = 5d-3 ! mean size mutation

  real(DP), parameter :: K_DIFF = 3d-2 ! diffusion conctant of organisms

  ! organisms interact both socally and through resource competition.
  ! both interactions have a range
  real(DP), parameter :: RANGE_SOCIAL = 1d0 ! social interaction range
  real(DP), parameter :: RANGE_COMPETITION =  4d0 ! competition interaction range
  real(DP), parameter :: RANGE_DISPERSION = sqrt(2d0 * K_DIFF * DT) ! sd of the dispersion

  !  carrying capacity density of resource competition; length unit is RANGE_SOCIAL
  real(DP), parameter :: K_COMPETITION = 1d2

  ! costs and benefits
  real(DP), parameter :: COST = 1d0 ! cost per unit of altruism; set to 1 per definition
  ! the benefits of the altruism saturate
  real(DP), parameter :: MAX_BENEFIT = 2d0 ! upper limit of benefit
  real(DP), parameter :: BENEFIT = 5d-1  ! benefit per unit of altruism at the linear (low-altruism) regime

  ! >> ARENA
  integer, parameter :: N = 2**16!4!8192!32768!16384 ! linear size of the simulated frame. preferably power of 2

  ! real(DP), parameter :: TYPE_HIGH = 5d-2
  ! real(DP), parameter :: TYPE_LOW = 1d-2
  !-------------------------------------------------------------------
  !
  ! >> SIMULATION SETTINGS
  integer, parameter :: seed(33) = 1 ! random seed

  ! last timestep
  integer, parameter :: T_MAX = 2000001
  ! how frequently to report the state of the system, in timesteps
  integer, parameter :: REPORT_INTERVAL = 5000
  ! for some variables, we keep track of a rolling mean
  ! this parameter sets the size of the window in timesteps
  integer, parameter :: MEAN_INTERVAL = 1000 !

  ! at the end of the simulation, we calculate S_local and S_interlocal for various scales sigma
  ! we do this for many time steps and average the result
  ! we start NR_AV time steps before the end of the simulation
  ! and calculate the S_local-curve every SC_IN time steps
  ! making sure to always end at time MAX_T - 1.
  ! (We cannot calculate S_local curves at MAX_T because we need
  ! realized fitnesses for that.)
  integer, parameter :: NR_AV = 0 !
  integer, parameter :: SC_IN = 20 ! interval between scurves
  integer, parameter :: NR_SC = (NR_AV-1)/SC_IN + 1  ! number of S-curves to be computed

  integer, parameter :: NR_SCALES = 82     ! for how many scales sigma do we compute S_local?
  !real(DP), parameter :: S_FROM = 1d0  ! lowest scale
  !real(DP), parameter :: S_TO = 1d1     ! largest scale


  ! we use periodic boundary conditions
  ! this means that the kernel used for calculations of S_local
  ! must take a certain number number of periodic images into account
  ! we look NR_FIELDS left, and right
  integer, parameter :: NR_FIELDS = 7


  ! the Price equation involves a timesteps
  ! this is the timestep used for the analysis of collectives
  integer, parameter :: ANALYSIS_INTERVAL = 1000

  ! to identify collectives, we smoothen the occupacy matrix
  ! here's the bandwidth
  real(DP), parameter :: WIDTH_SMOOTH = 1d0
  ! optionally, the occupancy can also be smoothed in time
  ! if this parameter is set to 0,
  ! the smoothed density only depends on the occupancy
  ! at the current time
  real(DP), parameter :: TIME_SMOOTHING = 0d0
  ! threshold density for a break point beteen collectives to considered lost
  real(DP), parameter :: K_THRES_WEAK = K_COMPETITION*7d-1
  ! threshold density for a break point beteen collectives to be gained
  real(DP), parameter :: K_THRES_STRICT = K_COMPETITION/5d0
  ! maximum number of collectives counted on.  Don't exceed!
  integer(DP), parameter :: MAX_COLL = 1000

  ! resolution of the smoothed field exported
  integer(DP), parameter :: RES_SMOOTH = 4000

  ! initial conditions
  ! initial population size in units of the estimated steady-state density SS
  real(DP), parameter :: INIT_POP = 1d0

contains

end module parameters
