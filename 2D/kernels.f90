module kernels
  use parameters
  use shorthands
  use fft
  use pairwise

  implicit none
  complex(DP) :: fft_kernel(N,N) = 0   ! FFT of competition kernel
  complex(DP) :: fft_kernels(N, N, NR_SCALES) = 0d0 ! FFT used for calculation of S_local and S_interlocal
  complex(DP) :: fft_kernelAltr(N,N) = 0 ! FFT of altruism kernel
  complex(DP) :: fft_kernelLSD(N,N) = 0 ! FFT of kernel used for visualizing Kernel Selection Differential
  real(DP) :: S_curve_vals(NR_SC, NR_SCALES + 1, 3)   ! value of S_local curves for various scales

contains

  subroutine assemble_kernels()
    write(*, '(a)') "# Assembling kernels... "
    write(*, '(a)', advance = "no") "    Competition kernel... "
    call compute_kernel(fft_kernel, RANGE_COMPETITION)
    print *, "Done."
    write(*, '(a)', advance = "no") "    Altruism kernel... "
    call compute_kernel(fft_kernelAltr, RANGE_SOCIAL)
    print *, "Done."
    write(*, '(a)', advance = "no") "    Multiscale selection kernels... "
    if (STEP_FUNCTION_KERNEL) then
      call compute_kernel_stepfunction(fft_kernelLSD, RANGE_CONTR_SB)
    else
      call compute_kernel(fft_kernelLSD, RANGE_CONTR_SB)
    end if
    call precompile_all_kernels()
    print *, "Done."
  end subroutine assemble_kernels

  subroutine compute_kernel(fft_kernelA, sigA)
    complex(DP), intent(out) :: fft_kernelA(N,N)
    real(DP), intent(in) :: sigA
    real(DP) :: kernel1D(N) = 0d0
    real(DP) :: pre
    integer :: col, hfield, j
    real(DP) :: tmp(2*NR_FIELDS + 1)


    pre =  1d0/(2d0*(sigA*RESOLUTION)**2d0)

    do col = 1,N
      do hfield = -NR_FIELDS, NR_FIELDS
        tmp(hfield + NR_FIELDS + 1) = exp(&
        & -pre*(real(col - 1 + hfield*N, kind = DP)**2d0)&
        & )
      end do
      kernel1D(col) = pSum(tmp)
    end do

    forall(j=1:N)
      fft_kernelA(:,j) = cmplx(kernel1D(:) * kernel1D(j), kind = DP)
    end forall
    fft_kernelA = (RESOLUTION**2d0)*fft_kernelA/pSumC_2D(fft_kernelA)

    call cfft2f(N, N, N, fft_kernelA, wsave, lensav, work, lenwrk, ier)

  end subroutine compute_kernel

  subroutine precompile_all_kernels()
    integer :: i

    call choose_s_curve_vals()

    do i = 1, NR_SCALES
      if (STEP_FUNCTION_KERNEL) then
        call compute_kernel_stepfunction(fft_kernels( :, :, i), S_curve_vals (1, i + 1, 1))
      else
        call compute_kernel(fft_kernels( :, :, i), S_curve_vals (1, i + 1, 1))
      end if
    end do

  end subroutine precompile_all_kernels

  subroutine compute_kernel_stepfunction(fft_kernelA, sig)
    complex(8), intent(out) :: fft_kernelA(N,N)
    real(8), intent(in) :: sig
    integer :: col, row
    integer :: drow, dcol

    if (sig > real(HALF_N, kind = DP)) then
      print *, "ERROR!!"
      print *, "compute_kernel_stepfunction() is incorrect if sig > N/2."
      call exit(0)
    end if


    do col = 1,N
      do row = 1, N
        drow = min(row -1, N - row +1)
        dcol = min(col -1, N - col +1)
        if (drow**2 + dcol**2 < floor((sig*RESOLUTION)**2d0) + 1) then
          fft_kernelA(row, col) = 1d0
        else
          fft_kernelA(row, col) = 0d0
        end if
      end do
    end do

    fft_kernelA = (RESOLUTION**2d0)*fft_kernelA/pSumC_2D(fft_kernelA)

    call cfft2f(N, N, N, fft_kernelA, wsave, lensav, work, lenwrk, ier)

  end subroutine compute_kernel_stepfunction

  subroutine choose_s_curve_vals()
    integer :: i

    !routine that sets the scales goes here.

    S_curve_vals(:, 1, 1) = 0d0

    do i = 1, 9
      S_curve_vals (:, i + 1, 1) = real(i, kind = DP)/1d1
    end do

    do i = 10, 21
      S_curve_vals (:, i + 1, 1) = 1d0 + real(i - 10, kind = DP)/3d0
    end do

    do i = 22, 31
      S_curve_vals (:, i + 1, 1) = 5d0 + real(i - 22, kind = DP)/2d0
    end do

    do i = 32, 42
      S_curve_vals (:, i + 1, 1) = 1d1 + real(i - 32, kind = DP)*1d0
    end do

  end subroutine

end module kernels
