module multiscale
  use parameters
  use shorthands
  use pairwise
  use kernels
  use data_structures
  use fft
  use price

  implicit none

  real(DP) :: S_curves_mean(NR_SCALES + 1, 3) = 0d0   ! mean of the S_local curves for various scales (1st dimension)
  real(DP) :: S_curves_sd(NR_SCALES + 1) = 0d0        ! standard deviation of S curves
  integer :: nr_SC_C = 0                              ! number of S-curves already computed

contains

  subroutine compute_S_curve(field, o_stats)
    type(fld), intent(in) :: field
    type(org_sts), intent(in) :: o_stats
    complex(DP) :: fft_occ_field(N,N)
    complex(DP) :: fft_p_field(N,N)
    complex(DP) :: fft_w_field(N,N)
    integer :: i

    nr_SC_C = nr_SC_C + 1
    print *, "This will be S-curve number ", nr_SC_C, " out of ", NR_SC

    fft_occ_field = cmplx(field%occ, kind = DP)
    call cfft2f(N, N, N, fft_occ_field, wsave, lensav, work, lenwrk, ier)

    fft_p_field = cmplx(field%p, kind = DP)
    call cfft2f(N, N, N, fft_p_field, wsave, lensav, work, lenwrk, ier)

    fft_w_field = cmplx(field%w, kind = DP)
    call cfft2f(N, N, N, fft_w_field, wsave, lensav, work, lenwrk, ier)

    ! calculations for sigma = 0
    S_curve_vals(nr_SC_C, 1, 1) = 0d0
    S_curve_vals(nr_SC_C, 1, 2) = compute_S_interlocal_0(field, o_stats)
    S_curve_vals(nr_SC_C, 1, 3) = S - S_curve_vals(nr_SC_C, 1, 2)
    S_vals(nr_SC_C) = S

    ! calculations for all other values of sigma
    do i = 2, (NR_SCALES + 1)
      call compute_S_split(i - 1, fft_occ_field, fft_p_field, fft_w_field, S_curve_vals(nr_SC_C, i, 2:3), o_stats)
    end do

  end subroutine compute_S_curve

  real(DP) function compute_S_interlocal_0(field, o_stats) result (s0)
    type(fld), intent(in) :: field
    type(org_sts), intent(in) :: o_stats
    real(DP) :: tmp(N, N)
    real(DP) :: occR(N, N)

    occR = real(field%occ, kind = DP)

    tmp = 0d0

    where (field%occ > 0)
      tmp = (field%p - occR*o_stats%mean_p) * (field%w - occR*o_stats%mean_w) / occR
    end where

    s0 = pSum_2D(tmp)/real(o_stats%n, kind = DP)
  end function compute_S_interlocal_0

  subroutine compute_S_split(kn, fft_occ_field, fft_p_field, fft_w_field, ssplit, o_stats)
    integer, intent(in) :: kn ! number of the kernel
    complex(DP), intent(in)  :: fft_occ_field(N,N)
    complex(DP), intent(in) :: fft_p_field(N,N)
    complex(DP), intent(in) :: fft_w_field(N,N)
    real(DP), intent(out) :: ssplit(2)
    type(org_sts), intent(in) :: o_stats
    complex(DP) :: d_est(N,N)
    complex(DP) :: p_est_unnorm(N,N)
    complex(DP) :: w_est_unnorm(N,N)
    logical :: mymask(N, N)

    ! use kernel to compute the corresponding kernel density estimate
    d_est = fft_kernels(:, :, kn)*fft_occ_field*cmplx(N**2, kind = DP)
    call cfft2b(N, N, N, d_est, wsave, lensav, work, lenwrk, ier)

    ! compute unnormalized mean phenotype estimate
    p_est_unnorm = fft_kernels(:, :, kn)*fft_p_field*cmplx(N**2, kind = DP)
    call cfft2b(N, N, N, p_est_unnorm, wsave, lensav, work, lenwrk, ier)

    ! compute unnormalized mean fitness estimate
    w_est_unnorm = fft_kernels( :, :, kn)*fft_w_field*cmplx(N**2, kind = DP)
    call cfft2b(N, N, N, w_est_unnorm, wsave, lensav, work, lenwrk, ier)

    where (real(d_est, kind = DP) > 1d-10)
      mymask = .true.
    elsewhere
      mymask = .false.
    end where

    ! s interlocal
    ssplit(1) = pSum(&
    & pack( &
    & real(w_est_unnorm- d_est*o_stats%mean_w, kind = DP),&
    & mymask) * &
    & pack( &
    & real((p_est_unnorm-d_est*o_stats%mean_p)/d_est, kind = DP), &
    & mymask) &
    & )/real(o_stats%n, kind = DP)
    ! s local
    ssplit(2) = S - ssplit(1)

  end subroutine compute_S_split

  subroutine compute_contr_Slocal(field)
    type(fld), intent(inout) :: field
    complex(DP) :: fft_occ_field(N,N)
    complex(DP) :: fft_p_field(N,N)
    complex(DP) :: fft_w_field(N,N)
    complex(DP) :: fft_w_p_field(N,N)
    complex(DP) :: d_est(N,N)
    complex(DP) :: p_est_unnorm(N,N)
    complex(DP) :: w_est_unnorm(N,N)
    complex(DP) :: w_p_est_unnorm(N,N)

    fft_occ_field = cmplx(field%occ, kind = DP)
    call cfft2f(N, N, N, fft_occ_field, wsave, lensav, work, lenwrk, ier)

    fft_p_field = cmplx(field%p, kind = DP)
    call cfft2f(N, N, N, fft_p_field, wsave, lensav, work, lenwrk, ier)

    fft_w_field = cmplx(field%w, kind = DP)
    call cfft2f(N, N, N, fft_w_field, wsave, lensav, work, lenwrk, ier)

    fft_w_p_field = cmplx(field%w_p, kind = DP)
    call cfft2f(N, N, N, fft_w_p_field, wsave, lensav, work, lenwrk, ier)

    ! use kernel to compute the corresponding kernel density estimate
    d_est = fft_kernelLSD*fft_occ_field*cmplx(N**2, kind = DP)
    call cfft2b(N, N, N, d_est, wsave, lensav, work, lenwrk, ier)

    ! compute unnormalized mean phenotype estimate
    p_est_unnorm = fft_kernelLSD*fft_p_field*cmplx(N**2, kind = DP)
    call cfft2b(N, N, N, p_est_unnorm, wsave, lensav, work, lenwrk, ier)

    ! compute unnormalized mean fitness estimate
    w_est_unnorm = fft_kernelLSD*fft_w_field*cmplx(N**2, kind = DP)
    call cfft2b(N, N, N, w_est_unnorm, wsave, lensav, work, lenwrk, ier)

    ! compute unnormalized mean product of fitness and phenotype estimate
    w_p_est_unnorm = fft_kernelLSD*fft_w_p_field*cmplx(N**2, kind = DP)
    call cfft2b(N, N, N, w_p_est_unnorm, wsave, lensav, work, lenwrk, ier)

    where (real(d_est, kind = DP) > 1d-10)
      field%contrib_Slocal = real(w_p_est_unnorm - w_est_unnorm*p_est_unnorm/d_est, kind = DP)
    elsewhere
      field%contrib_Slocal = 0d0
    end where

  end subroutine compute_contr_Slocal

  subroutine update_w_field()
    integer :: i

    field(py)%w = 0d0  !wipe
    do i = 1, o_stats(py)%n
      field(py)%w(o_list(i, py)%y, o_list(i, py)%x) = &
      & field(py)%w(o_list(i, py)%y, o_list(i, py)%x) + o_list(i, py)%w
    end do

  end subroutine update_w_field


  subroutine update_w_p_field()
    integer :: i

    field(py)%w_p = 0d0  !wipe
    do i = 1, o_stats(py)%n
      field(py)%w_p(o_list(i, py)%y, o_list(i, py)%x) = &
      & field(py)%w_p(o_list(i, py)%y, o_list(i, py)%x) + o_list(i, py)%w*o_list(i, py)%p
    end do

  end subroutine update_w_p_field

  subroutine compute_means()
    integer :: i

    S_mean = pSum(S_vals(1:nr_SC_C)) / real(nr_SC_C, kind = DP)

    do i = 1, (NR_SCALES + 1)
      S_curves_mean(i, 2) = pSum(S_curve_vals(1:nr_SC_C, i, 2))/ real(nr_SC_C, kind = DP)
      S_curves_mean(i, 3) =  pSum(S_curve_vals(1:nr_SC_C, i, 3)) / real(nr_SC_C, kind = DP)
      S_curves_mean(i, 1) = S_curve_vals(1, i,1)
      S_curves_sd(i) = sqrt(simple_covariance(S_curve_vals(1:nr_SC_C, i, 3), S_curve_vals(1:nr_SC_C, i, 3), nr_SC_C))
    end do
  end subroutine compute_means


  subroutine output_S_curves()
    integer :: i, j

    open(unit=13, file="S_curves.txt", form="formatted", action="write", status="replace")
    write(13,*) "# writes the S_curves obtained at different time points as columns"
    write(13, '(a)', advance = "no") "scale, "
    do i = 1,nr_SC_C
      write(13, '(a)', advance = "no") "S_local"
      write(13, '(g0)', advance = "no") i
      write(13, '(a)', advance = "no") ", "
    end do
    do i = 1,(nr_SC_C -1)
      write(13, '(a)', advance = "no") "S_interlocal"
      write(13, '(g0)', advance = "no") i
      write(13, '(a)', advance = "no") ", "
    end do
    write(13, '(a)', advance = "no") "S_interlocal"
    write(13, '(g0)') nr_SC_C
    do i=1,(NR_SCALES + 1)
      write(13, "(*(g0, :, ', '))") S_curve_vals(1, i, 1), &
      & (S_curve_vals(j, i, 3), j= 1,nr_SC_C), (S_curve_vals(j, i, 2), j= 1,nr_SC_C)
    end do
    close(unit=13)

    open(unit=13, file="S_curves_mean.txt", form="formatted", action="write", status="replace")
    write(13, "(g0)") "scale, mean S_local, mean S_interlocal, sd S_(inter)local, mean S_interlocal + mean S_local, S_mean"
    do i=1,(NR_SCALES + 1)
      write(13, "(*(g0, :, ', '))") S_curves_mean(i, 1), S_curves_mean(i, 3), S_curves_mean(i, 2),&
      & S_curves_sd(i), S_curves_mean(i, 2) + S_curves_mean(i, 3), S_mean
    end do
    close(unit=13)

  end subroutine output_S_curves

end module multiscale
