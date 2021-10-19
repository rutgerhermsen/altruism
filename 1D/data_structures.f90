! Rutger Hermsen, 2017
module data_structures
  use parameters

  implicit none


  !-------------------------------------------------------------------
  ! >>>>>> DATA STRUCTURES
  !
  ! data structure containing info about one organism
  type org
    sequence
    integer :: x = 0 ! x coordinate
    real(DP) :: p = 0d0! phenotype (level of altruism)
    real(DP) :: p_e = 0d0 ! effective phenotype -- equals p unless p is negative
    real(DP) :: gr = 0d0  ! growth rate of each organism
    integer :: offspring = 0 ! realized absolute fitness
    integer :: offspring_MLS1 = 0 ! realized absolute fitness over time step of MLS1
    real(DP) :: W_abs = 0d0 ! absolute fitness (realized or expected, per settings)
    real(DP) :: w = 0d0 ! relative fitness (realized or expected, per settings)
    real(DP) :: w_MLS1 = 0d0 ! relative fitness (realized) over time step of MLS1
    real(DP) :: tot_Delta_p = 0d0 ! sum of changes in phenotype of offspring due to mutation
    integer :: c = 1
    integer :: a = 0 ! ancestor in previous population MLS1
  end type org

  ! data structure containing info about one collective
  type coll
    sequence
    real(DP) :: x = 1d0  ! center of mass of collective
    integer :: l = 0! left boundary
    integer :: r = 0 ! right boundary
    integer :: minr = 1! location of densitiy nimimum associated with the right border
    real(DP) :: p = 0d0 ! mean phenotype
    integer :: offspring_MLS1 = 0 ! realized absolute fitness of type MLS1
    integer :: offspring_MLS2 = 0 ! realized absolute fitness of type MLS2
    real(DP) :: w_MLS1 = 0d0 ! relative realized fitness of type MLS1
    real(DP) :: w_MLS2 = 0d0 ! relativer realized fitness of type MLS2
    integer :: parent = 0
    real(DP) :: tot_Delta_p = 0d0 ! sum of changes in phenotype of offspring
    integer :: n = 0 ! number of organisms in this collective
    integer :: c_id = 0
  end type coll

  type fld
    sequence
    integer :: occ(N) = 0 !  occupancy (number of organisms on each site)
    real(DP) :: p(N) = 0d0 ! total phenotype (altruism)
    real(DP) :: p_e(N) = 0d0 ! total effective phenotype (altruism) per grid cell
    real(DP) :: w(N) = 0d0       ! sum of fitnesses at each site
    real(DP) :: w_p(N) = 0d0   ! sum of the products of w and p
    complex(DP) :: K_dens(N) = 0d0    ! kernel density estimate for resource competition
    complex(DP) :: K_altr(N) = 0d0    ! unnormalized kernel mean of phenotype (altruism)
    real(DP) :: K_smooth(N) = 0d0    ! smoothed (kernel density) estimate for finding collectives
    complex(DP) :: contrib_Slocal(N) = 0d0  ! contribution to selection at some scale
  end type fld

  ! some population statistics
  type org_sts
    integer :: n = 0            ! total population size
    real(DP) :: mean_w = 0d0          ! population mean of relative fitness (1 by definition!)
    real(DP) :: mean_p = 0d0         ! population mean of phenotype
  end type org_sts

  type coll_sts
    integer :: n = 0            ! total population size
    real(DP) :: mean_w_MLS1 = 1d0       ! population mean of relative fitness (1 by definition!)
    real(DP) :: mean_w_MLS2 = 1d0       ! population mean of relative fitness (1 by definition!)
    real(DP) :: mean_p = 0d0        ! population mean of mean phenotype
  end type coll_sts

contains


end module data_structures
