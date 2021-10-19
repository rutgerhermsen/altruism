module g_of_r
  use parameters
  use shorthands
  use fft

  implicit none

  integer(4) :: distances(0:(HALF_N**2))
  logical :: distance_mask(0:(HALF_N**2))

contains

  subroutine prepare_distances()
    integer :: col, row
    integer :: dcol, drow
    integer :: dis_sq

    write(*, '(a)', advance = "no") "# Counting distance frequencies for radial distribution function..."

    distances = 0
    do col = 1,N
      do row = 1, N
        drow = min(row -1, N - row +1)
        dcol = min(col -1, N - col +1)
        dis_sq = drow**2 + dcol**2
        if (dis_sq < HALF_N**2 + 1) then
          distances(dis_sq) = distances(dis_sq) + 1
        end if
      end do
    end do

    where (distances == 0)
      distance_mask = .false.
    elsewhere
      distance_mask = .true.
    end where
    print *, "Done."

  end subroutine


  subroutine compute_g_of_r(autocorr, g_of_r)
    complex(DP), intent(in) :: autocorr(N,N)
    real(DP), intent(out):: g_of_r(0:HALF_N**2)
    integer :: col, row
    integer :: dcol, drow
    integer :: dis_sq


    g_of_r = 0d0
    do col = 1,N
      do row = 1, N
        drow = min(row -1, N - row +1)
        dcol = min(col -1, N - col +1)
        dis_sq = drow**2 + dcol**2
        if (dis_sq < HALF_N**2 + 1) then
          g_of_r(dis_sq) = g_of_r(dis_sq) + real(autocorr(row, col), kind = DP)
        end if
      end do
    end do

  end subroutine

  subroutine compute_autocorrelation(autocorr, occ, nA)
    complex(DP), intent(out) :: autocorr(N,N)
    integer, intent(in) :: occ(N, N)
    integer, intent(in) :: nA
    complex(DP) :: fft_occ_field(N,N)

    fft_occ_field = cmplx(occ, kind = DP)

    call cfft2f(N, N, N, fft_occ_field, wsave, lensav, work, lenwrk, ier)

    autocorr = fft_occ_field * conjg(fft_occ_field)
    call cfft2b(N, N, N, autocorr, wsave, lensav, work, lenwrk, ier)
    autocorr = autocorr*(cmplx(N**2, kind = DP)/cmplx(nA, kind = DP))**2d0

  end subroutine

end module g_of_r
