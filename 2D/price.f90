module price
  use parameters
  use shorthands
  use data_structures
  use state
  use utils

  implicit none

  real(DP) :: S = 0d0                ! directional selection
  real(DP) :: S_pur, S_pur_alt = 0d0 ! purifying selection
  real(DP) :: S_mean = 0d0           ! mean of selection
  real(DP) :: S_vals(NR_SC) = 0d0    ! values of S at various times
  real(DP) :: S_r_mean(MEAN_INTERVAL) = 0d0      ! rolling mean of S
  real(DP) :: S_pur_r_mean(MEAN_INTERVAL) = 0d0  ! rolling mean purifying selection
  real(DP) :: S_cum = 0d0, S_pur_cum = 0d0       ! cumulative effect of S and purif. S.
  real(DP) :: transmission = 0d0     ! transmission
  real(DP) :: transmission_r_mean(MEAN_INTERVAL) = 0d0 ! rolling mean of transmission
  real(DP) :: transmission_cum = 0d0 ! cumulative effect of transmission
  real(DP) :: drift = 0d0
  real(DP) :: drift_cum = 0d0        ! cumulative effect of drift

contains


  subroutine calculate_price(t)
    integer, intent(in) :: t

    call update_selection(t)
    call calculate_transmission(t)
    call calculate_drift()

  end subroutine calculate_price

  subroutine update_selection(t)
    integer, intent(in) :: t
    real(DP) :: var_p

    if (REALIZED_FITNESS) then
      call update_realized_w_list()
      o_stats(py)%mean_w = 1d0
    else
      call update_w_list()
      o_stats(py)%mean_w = pSum(o_list(1:o_stats(py)%n, py)%w) &
      & / real(o_stats(py)%n, kind = DP)
    end if
    S = simple_covariance( &
    & o_list(1:o_stats(py)%n, py)%w - o_stats(py)%mean_w, & ! subtract the mean for stability (should be ~1d0)
    & o_list(1:o_stats(py)%n, py)%p - o_stats(py)%mean_p, & ! subtract the mean for stability
    & o_stats(py)%n &
    & )

    S_r_mean(modulo(t - 1, MEAN_INTERVAL) + 1) = S
    S_cum = S_cum + S

    var_p = simple_covariance(o_list(1:o_stats(py)%n, py)%p, o_list(1:o_stats(py)%n, py)%p, o_stats(py)%n)

    S_pur = simple_covariance( &
    & o_list(1:o_stats(py)%n, py)%w - o_stats(py)%mean_w, & ! subtract the mean for stability (should be 1d0)
    & (o_list(1:o_stats(py)%n, py)%p - o_stats(py)%mean_p)**2d0 - var_p, & ! subtract the mean for stability
    & o_stats(py)%n &
    & )


    S_pur_alt = S_pur - S**2d0

    S_pur_r_mean(modulo(t - 1, MEAN_INTERVAL) + 1) = S_pur
    S_pur_cum = S_pur_cum + S_pur

  end subroutine update_selection

  subroutine update_realized_w_list()

    o_list(1:o_stats(py)%n, py)%W_abs = real(o_list(1:o_stats(py)%n, py)%offspring, kind = DP)
    o_list(1:o_stats(py)%n, py)%w = o_list(1:o_stats(py)%n, py)%W_abs*real(o_stats(py)%n, kind = DP)/ &
    & pSum(o_list(1:o_stats(py)%n, py)%W_abs)

  end subroutine update_realized_w_list

  subroutine update_w_list()

    o_list(1:o_stats(py)%n, py)%W_abs = (1d0 - P_D)*(o_list(1:o_stats(py)%n, py)%gr + 1d0)
    o_list(1:o_stats(py)%n, py)%w = o_list(1:o_stats(py)%n, py)%W_abs*&
    &real(o_stats(py)%n, kind = DP)/pSum(real(o_list(1:o_stats(py)%n, py)%offspring, kind = DP))
  end subroutine

  subroutine calculate_transmission(t)
    integer, intent(in) :: t

    transmission = pSum(o_list(1:o_stats(py)%n, py)%tot_Delta_p) / real(o_stats(py_N)%n, kind = DP)
    transmission_r_mean(modulo(t - 1, MEAN_INTERVAL) + 1) = transmission
    transmission_cum = transmission_cum + transmission

  end subroutine calculate_transmission

  subroutine calculate_drift()
    real(DP) :: mismatch(MAX_POP)
    real(DP) :: mean_mismatch

    mismatch(1:o_stats(py)%n) = real(o_list(1:o_stats(py)%n, py)%offspring, kind = DP) &
    & - (1d0 - P_D)*(o_list(1:o_stats(py)%n, py)%gr)
    mean_mismatch = pSum(mismatch)/real(o_stats(py)%n, kind = DP)
    ! note that we have omitted a constant, which below would be subtracted anyway.

    drift = simple_covariance( &
    & mismatch(1:o_stats(py)%n) - mean_mismatch, & ! subtract the mean for stability
    & o_list(1:o_stats(py)%n, py)%p - o_stats(py)%mean_p, & ! subtract the mean for stability
    & o_stats(py)%n &
    & ) * real(o_stats(py)%n, kind = DP)/pSum(real(o_list(1:o_stats(py)%n, py)%offspring, kind = DP))

    drift_cum = drift_cum + drift

  end subroutine calculate_drift

end module price
