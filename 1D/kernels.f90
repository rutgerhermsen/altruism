module kernels
  use parameters
  use shorthands
  use fft
  use pairwise

  implicit none

  complex(DP) :: fft_kernel(N) = 0d0              ! FFT of competition kernel
  complex(DP) :: fft_kernelAltr(N) = 0d0          ! FFT of altruism kernel
  complex(DP) :: fft_kernelLSD(N) = 0d0           ! FFT of kernel used for visualizing Kernel Selection Differential
  complex(DP) :: fft_smooth_kernel(N) = 0d0       ! FFT of kernel used for visualizing smoothed density
  complex(DP) :: fft_kernels(N, NR_SCALES) = 0d0  ! FFT used for calculation of S_local and S_interlocal
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
    write(*, '(a)', advance = "no") "    Smoothing kernel... "
    call compute_kernel(fft_smooth_kernel, WIDTH_SMOOTH)
    print *, "Done."
    call compute_kernel(fft_kernelAltr, RANGE_SOCIAL)
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
    complex(DP), intent(out) :: fft_kernelA(N)
    real(DP), intent(in) :: sigA
    real(DP) :: pre
    integer :: col, hfield
    real(DP) :: tmp(2*NR_FIELDS + 1)


    pre =  1d0/(2d0*(sigA*RESOLUTION)**2d0)

    do col = 1,N
      do hfield = -NR_FIELDS, NR_FIELDS
        tmp(hfield + NR_FIELDS + 1) = exp(-pre*((real(col - 1 + hfield*N, kind = DP))**2d0))
      end do
      fft_kernelA(col) = pSum(tmp)
    end do

    fft_kernelA = RESOLUTION*fft_kernelA/pSumC(fft_kernelA)

    call cfft1f(N, 1, fft_kernelA, N, wsave, lensav, work, lenwrk, ier)

  end subroutine compute_kernel

  subroutine precompile_all_kernels()
    integer :: i

    call choose_s_curve_vals()

    do i = 1, NR_SCALES
      if (STEP_FUNCTION_KERNEL) then
        call compute_kernel_stepfunction(fft_kernels( :, i), S_curve_vals(1, i + 1, 1))
      else
        call compute_kernel(fft_kernels( :, i), S_curve_vals(1, i + 1, 1))
      end if
    end do

  end subroutine precompile_all_kernels

  subroutine compute_kernel_stepfunction(fft_kernelA, sig)
    complex(8), intent(out) :: fft_kernelA(N)
    real(8), intent(in) :: sig
    integer :: i
    integer :: di

    if (sig > real(HALF_N, kind = DP)) then
      print *, "ERROR!!"
      print *, "compute_kernel_stepfunction() is incorrect if sig > N/2."
      call exit(0)
    end if


    do i = 1,N
      di = min(i -1, N - i +1)
      if (di < sig*RESOLUTION + 1) then
        fft_kernelA(i) = 1d0
      else
        fft_kernelA(i) = 0d0
      end if
    end do

    fft_kernelA = (RESOLUTION)*fft_kernelA/pSumC(fft_kernelA)

    call cfft1f(N, 1, fft_kernelA, N, wsave, lensav, work, lenwrk, ier)

  end subroutine compute_kernel_stepfunction

  subroutine choose_s_curve_vals()
    integer :: i
    !routine that sets the scales goes here.

    S_curve_vals(:, 1, 1) = 0d0

    S_curve_vals(:, 2, 1) = (1d0/3d0)*1d-1
    S_curve_vals(:, 3, 1) = (2d0/3d0)*1d-1

    do i = 3, 9
      S_curve_vals (:, i + 1, 1) = 1d-1 + (2d0/3d0)*real(i - 3, kind = DP)*1d-1
    end do

    do i = 10, 14
      S_curve_vals (:, i + 1, 1) = 5d-1 + 1d0*real(i - 9, kind = DP)*1d-1
    end do

    do i = 15, 32
      S_curve_vals (:, i + 1, 1) = 1d0 + (5d0/3d0)*real(i - 14, kind = DP)*1d-1
    end do

    do i = 33, 52
      S_curve_vals (:, i + 1, 1) = 4d0 + (6d0/3d0)*real(i - 32, kind = DP)*1d-1
    end do

    do i = 53, 72
      S_curve_vals (:, i + 1, 1) = 8d0 + (9d0/3d0)*real(i - 52, kind = DP)*1d-1
    end do

    do i = 73, 82
      S_curve_vals (:, i + 1, 1) = 1.4d1 + (1.2d1/3d0)*real(i - 72, kind = DP)*1d-1
    end do

    S_curve_vals(:, :, 1) = 2d0*S_curve_vals(:, :, 1)

  end subroutine

end module kernels
