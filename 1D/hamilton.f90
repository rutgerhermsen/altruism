module hamilton
  use parameters
  use shorthands
  use pairwise
  use data_structures
  use state
  use utils

  implicit none

  ! This module calculates for a given time step the three terms
  ! of Hamilton's Rule: cost, benefit, relatedness.
  ! To do so, we use the formalism due to Queller
  ! using multiple linear regression analysis.
  ! Apart from the trait value the "experienced level of altruism" is used
  ! as the second regression variable.
  ! However, we subtract the contribution by the focal organism itself.

  real(DP) :: coef_cost = 0d0  ! partial regression coefficient measuring cost
  real(DP) :: coef_cost_cum = 0d0  ! cumulative
  real(DP) :: coef_cost_r_mean(MEAN_INTERVAL) = 0d0  ! rolling mean

  real(DP) :: coef_ben = 0d0   ! partial regression coefficient measuring benefit
  real(DP) :: coef_ben_cum = 0d0   ! cumulative
  real(DP) :: coef_ben_r_mean(MEAN_INTERVAL) = 0d0   ! rolling mean

  real(DP) :: relatedness = 0d0 ! measure of relatedness
  real(DP) :: relatedness_cum = 0d0 ! cumulative
  real(DP) :: relatedness_r_mean(MEAN_INTERVAL) = 0d0 ! rolling mean


  real(DP) :: effect_cost = 0d0 ! effect of cost on mean: cost * var
  real(DP) :: effect_cost_cum = 0d0 ! cumulative
  real(DP) :: effect_cost_r_mean(MEAN_INTERVAL) = 0d0 ! rolling mean

  real(DP) :: effect_rel_ben = 0d0 ! relatedness * benefit * var
  real(DP) :: effect_rel_ben_cum = 0d0 ! cumulative
  real(DP) :: effect_rel_ben_r_mean(MEAN_INTERVAL) = 0d0 ! rolling mean


contains

  subroutine calculate_hamilton(t)
    integer, intent(in) :: t

    real(DP) :: cov_p_p = 0d0 ! variance of trait value p
    real(DP) :: cov_A_A = 0d0 ! variance of experienced level of altruism A
    real(DP) :: cov_p_A = 0d0 ! covariance between p and A
    real(DP) :: cov_w_p = 0d0 ! covariance between w and p
    real(DP) :: cov_w_A = 0d0 ! covariance between w and A
    integer :: o

    real(DP) :: context(o_stats(py)%n)  ! level of altruism experienced due to others
    real(DP) :: mean_context = 0d0

    do o = 1,o_stats(py)%n
      context(o) = real( &
        & field(py)%K_altr(o_list(o,py)%x), kind = DP &
        & ) - &
        & o_list(o, py)%p_e / sqrt( 2d0 * PI )
    end do
    mean_context = pSum(context) &
      & / real(o_stats(py)%n, kind = DP)

    cov_p_p = simple_covariance( &
      & o_list(1:o_stats(py)%n, py)%p - o_stats(py)%mean_p, &
      & o_list(1:o_stats(py)%n, py)%p - o_stats(py)%mean_p, &
      & o_stats(py)%n &
      & )
    cov_A_A = simple_covariance( &
      & context - mean_context, &
      & context - mean_context, &
      & o_stats(py)%n &
      & )
    cov_p_A = simple_covariance( &
      & o_list(1:o_stats(py)%n, py)%p - o_stats(py)%mean_p, &
      & context - mean_context, &
      & o_stats(py)%n &
      & )
    cov_w_p = simple_covariance( &
      & o_list(1:o_stats(py)%n, py)%w - o_stats(py)%mean_w, &
      & o_list(1:o_stats(py)%n, py)%p - o_stats(py)%mean_p, &
      & o_stats(py)%n &
      & )
    cov_w_A = simple_covariance( &
      & o_list(1:o_stats(py)%n, py)%w - o_stats(py)%mean_w, &
      & context - mean_context, &
      & o_stats(py)%n &
      & )

    coef_cost = (cov_w_p*cov_A_A - cov_w_A*cov_p_A) / (cov_p_p*cov_A_A - cov_p_A*cov_p_A)
    coef_cost_cum = coef_cost_cum + coef_cost
    coef_cost_r_mean(modulo(t - 1, MEAN_INTERVAL) + 1) = coef_cost

    coef_ben = (cov_w_A*cov_p_p - cov_w_p*cov_p_A) / (cov_p_p*cov_A_A - cov_p_A*cov_p_A)
    coef_ben_cum = coef_ben_cum + coef_ben
    coef_ben_r_mean(modulo(t - 1, MEAN_INTERVAL) + 1) = coef_ben

    relatedness = cov_p_A/cov_p_p
    relatedness_cum = relatedness_cum + relatedness
    relatedness_r_mean(modulo(t - 1, MEAN_INTERVAL) + 1) = relatedness

    effect_cost = coef_cost*cov_p_p
    effect_cost_cum = effect_cost_cum + effect_cost
    effect_cost_r_mean(modulo(t - 1, MEAN_INTERVAL) + 1) = effect_cost

    effect_rel_ben = coef_ben*cov_p_A
    effect_rel_ben_cum = effect_rel_ben_cum + effect_rel_ben
    effect_rel_ben_r_mean(modulo(t - 1, MEAN_INTERVAL) + 1) = effect_rel_ben

  end subroutine calculate_hamilton

end module hamilton
