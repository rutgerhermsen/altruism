! Rutger Hermsen, 2017
!

module shorthands
  use parameters

  implicit none

  !-------------------------------------------------------------------
  ! >>>>>> SHORTHANDS

  real(DP) :: SS ! estimated steady-state density for non-altruistic population
  integer :: INIT_TOTPOP ! initial population size
  integer :: MAX_POP ! maximal population size; has to be an over-estimate
  integer, parameter :: HALF_N = N/2 !
  real(DP), parameter :: ANALYSIS_INTERVAL_REAL = real(ANALYSIS_INTERVAL, kind = DP)
  real(DP), parameter :: N_REAL = real(N, kind = DP)
  character(len=1), parameter :: TAB = achar(9)
  real(DP), parameter :: PI = 4.D0*DATAN(1.D0)


contains

  subroutine set_shorthands()
    write(*, '(a)', advance = "no") "# Defining shorthands..."
    SS = K_COMPETITION*(P_R - P_D)/P_R                  ! estimated steady-state density in absence of costs or altruism
    MAX_POP = nint(1.5d0*K_COMPETITION*N_REAL/RESOLUTION)       ! upper bound of maximum population size
    INIT_TOTPOP = nint(INIT_POP*SS*N_REAL/RESOLUTION) ! initial population size chosen
    print *, "Done."
  end subroutine set_shorthands

end module shorthands
