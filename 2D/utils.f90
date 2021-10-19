! Rutger Hermsen, 2017
!
module utils
  use parameters
  use pairwise

  implicit none

contains

  integer function mm(i) result (j)
    integer, intent(in) :: i

    j = modulo(i - 1, N) + 1

  end function mm

  real(DP) function simple_covariance(data1, data2, length) result (covar)
    integer, intent(in) :: length
    real(DP), intent(in) :: data1(length)
    real(DP), intent(in) :: data2(length)
    real(DP) :: sumx
    real(DP) :: sumy
    real(DP) :: sumxy

    sumx = pSum(data1)
    sumy = pSum(data2)
    sumxy = pSum(data1*data2)

    covar = (sumxy - sumx*sumy/real(length, kind = DP))/ real(length, kind = DP)

  end function simple_covariance

  real(DP) function random_exponential( ) result(rexp)
    real(DP) :: r

    call random_number(r)
    rexp = -log(r)

  end function random_exponential

end module utils
