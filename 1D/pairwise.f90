module pairwise
  implicit none

  integer,  parameter, private :: DP = selected_real_kind(12, 60)
  integer, parameter, private :: NORMALSUM = 51

contains

  ! One-dimensional arrays of double-precision reals.
  recursive real(DP) function pSum(array) result (s)
    real(DP), intent(in) :: array(:)
    integer :: n, j, m

    n = size(array)

    if (n < NORMALSUM) then
      s = array(1)
      do j = 2, n
        s = s + array(j)
      end do
    else
      m = floor(dble(n)/2d0)
      s = pSum(array(1:m)) + pSum(array((m+1):n))
    end if

  end function pSum

  ! Two-dimensional arrays of double-precision real numbers.
  ! The idea is to sum up each column using the 1D version
  ! and add up the results for each column.
  recursive real(DP) function pSum_2D(array) result (s)
    real(DP), intent(in) :: array(:, :)
    integer :: n2, m

    n2 = size(array, 2)

    if (n2 < 2) then
      s = pSum(array(:, 1))
    else
      m = floor(dble(n2)/2d0)
      s = pSum_2D(array(:, 1:m)) + pSum_2D(array(:, (m+1):n2))
    end if

  end function pSum_2D

  ! One-dimensional arrays of double-precision complex numbers.
  recursive complex(DP) function pSumC(array) result (s)
    complex(DP), intent(in) :: array(:)
    integer :: n, j, m

    n = size(array)

    if (n < NORMALSUM) then
      s = array(1)
      do j = 2, n
        s = s + array(j)
      end do
    else
      m = floor(dble(n)/2d0)
      s = pSumC(array(1:m)) + pSumC(array((m+1):n))
    end if

  end function pSumC

  ! Two-dimensional arrays of double-precision complex numbers.
  ! The idea is to sum up each column using the 1D version
  ! and add up the results
  recursive complex(DP) function pSumC_2D(array) result (s)
    complex(DP), intent(in) :: array(:, :)
    integer :: n2, m

    n2 = size(array, 2)

    if (n2 < 2) then
      s = pSumC(array(:, 1))
    else
      m = floor(dble(n2)/2d0)
      s = pSumC_2D(array(:, 1:m)) + pSumC_2D(array(:, (m+1):n2))
    end if

  end function pSumC_2D


end module pairwise
