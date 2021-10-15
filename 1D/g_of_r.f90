module g_of_r
  use parameters
  use shorthands
  use fft

  implicit none

  integer(4) :: distances(0:HALF_N)

contains


  subroutine prepare_distances()
    integer :: col
    integer :: dcol
    integer :: dis

    write(*, '(a)', advance = "no") "# Counting distance frequencies for radial distribution function..."
    distances = 0
    do col = 1,N
      dcol = min(col -1, N - col +1)
      dis = abs(dcol)
      if (dis < HALF_N + 1) then
        distances(dis) = distances(dis) + 1
      end if
    end do
    print *, "Done."

  end subroutine


  subroutine compute_g_of_r(autocorr, g_of_r)
    complex(DP), intent(in) :: autocorr(N)
    real(DP), intent(out):: g_of_r(0:HALF_N)
    integer :: col
    integer :: dcol
    integer :: dis


    g_of_r = 0d0
    do col = 1,N
      dcol = min(col -1, N - col +1)
      dis = abs(dcol)
      if (dis < HALF_N + 1) then
        g_of_r(dis) = g_of_r(dis) + real(autocorr(col), kind = DP)
      end if
    end do

  end subroutine

  subroutine compute_autocorrelation(autocorr, occ, nA)
    complex(DP), intent(out) :: autocorr(N)
    integer, intent(in) :: occ(N)
    integer, intent(in) :: nA
    complex(DP) :: fft_occ_field(N)

    fft_occ_field = cmplx(occ, kind = DP)

    call cfft1f(N, 1, fft_occ_field, N, wsave, lensav, work, lenwrk, ier)

    autocorr = fft_occ_field * conjg(fft_occ_field)
    call cfft1b(N, 1, autocorr, N, wsave, lensav, work, lenwrk, ier)
    autocorr = autocorr*(cmplx(N, kind = DP)/cmplx(nA, kind = DP))**2d0

  end subroutine

end module g_of_r
