! Rutger Hermsen, 2017
module data_structures
  use parameters

  implicit none

  !-------------------------------------------------------------------
  ! >>>>>> DATA STRUCTURES
  !
  ! define data structure containing info about one organism
  type org
    sequence
    integer :: x  = 0 ! x coordinate
    integer :: y = 0 ! y coordinate
    real(DP) :: p = 0d0 ! phenotype (level of altruism)
    real(DP) :: gr = 0d0  ! growth rate of each organism
    integer :: offspring = 0 ! realized fitness
    real(DP) :: W_abs = 0d0 ! absolute fitness (realized or expected, per settings)
    real(DP) :: w = 0d0 ! relative fitness (realized or expected, per settings)
    real(DP) :: tot_Delta_p = 0d0 ! sum of changes in phenotype of offspring due to mutation
    integer :: m = 0 ! marker that can be used for various things
  end type org

  ! some population statistics
  type org_sts
    integer :: n = 0            ! total population size
    real(DP) :: mean_w = 0d0          ! population mean of relative fitness (1 by definition!)
    real(DP) :: mean_p = 0d0         ! population mean of phenotype
  end type org_sts

  type fld
    sequence
    integer :: occ(N, N) = 0 !  occupancy (number of organisms on each site)
    real(DP) :: p(N, N) = 0d0 ! total phenotype (altruism)
    real(DP) :: w(N, N) = 0d0       ! sum of fitnesses at each site
    real(DP) :: w_p(N, N) = 0d0   ! sum of the products of w and p
    complex(DP) :: K_dens(N, N) = 0d0    ! kernel density estimate for resource competition
    complex(DP) :: K_altr(N, N) = 0d0    ! unnormalized kernel mean of phenotype (altruism)
    complex(DP) :: contrib_Slocal(N, N) = 0d0  ! contribution to selection at some scale
  end type fld

contains


end module data_structures
