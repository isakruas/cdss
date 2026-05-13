!> @file cdss_calibration.f90
!! @brief Generalized mutual information utilities for LLR calibration.
!!
!! @details
!! After a successful CRC check, decoded bits can be treated as known labels for
!! assessing LLR calibration. This module computes generalized mutual information
!! (GMI) and searches for a positive multiplicative temperature that improves the
!! match between LLR magnitude and empirical correctness.
!!
!! LLR sign convention: positive LLR favors bit 0. For a correctly signed LLR,
!! the signed product used by the GMI expression is positive.
module cdss_calibration
  use cdss_kinds, only: dp
  use cdss_constants, only: ln2, eps
  implicit none
  private
  public :: gmi_of_llrs, optimal_temperature

contains

  !> Estimate generalized mutual information for labeled LLRs.
  !!
  !! @param[in] llrs      Channel LLRs with positive values favoring bit 0.
  !! @param[in] true_bits Reference bits represented as 0 or 1.
  !! @return GMI estimate in bits per channel use.
  pure real(dp) function gmi_of_llrs(llrs, true_bits) result(gmi)
    real(dp), intent(in) :: llrs(:)
    integer,  intent(in) :: true_bits(:)
    integer :: i, n
    real(dp) :: acc, sgn, x
    n = min(size(llrs), size(true_bits))
    if (n == 0) then
      gmi = 0.0_dp
      return
    end if
    acc = 0.0_dp
    do i = 1, n
      ! Sign convention: LLR > 0 favors bit 0; sgn * L should be positive when correct.
      sgn = merge(1.0_dp, -1.0_dp, true_bits(i) == 0)
      x = sgn * llrs(i)
      acc = acc + log1p_exp_neg(x) / ln2
    end do
    gmi = 1.0_dp - acc / real(n, dp)
  end function gmi_of_llrs

  !> Find the positive LLR scale that maximizes GMI by golden-section search.
  !!
  !! @param[in] llrs      Channel LLRs with positive values favoring bit 0.
  !! @param[in] true_bits Reference bits represented as 0 or 1.
  !! @return Temperature beta that maximizes GMI(beta * llrs, true_bits).
  !!
  !! The one-dimensional search is performed on log(beta) in [-3, +3], which
  !! corresponds approximately to beta in [0.05, 20].
  function optimal_temperature(llrs, true_bits) result(beta_star)
    real(dp), intent(in) :: llrs(:)
    integer,  intent(in) :: true_bits(:)
    real(dp) :: beta_star
    real(dp) :: a, b, c, d, fc, fd, phi, tmp
    real(dp), allocatable :: scaled(:)
    integer :: i, iter

    a = -3.0_dp
    b =  3.0_dp
    phi = (sqrt(5.0_dp) - 1.0_dp) * 0.5_dp     ! 1/golden ratio
    c = b - phi * (b - a)
    d = a + phi * (b - a)

    allocate(scaled(size(llrs)))
    scaled = exp(c) * llrs
    fc = -gmi_of_llrs(scaled, true_bits)
    scaled = exp(d) * llrs
    fd = -gmi_of_llrs(scaled, true_bits)

    do iter = 1, 30
      if (fc < fd) then
        b = d
        d = c
        fd = fc
        c = b - phi * (b - a)
        scaled = exp(c) * llrs
        fc = -gmi_of_llrs(scaled, true_bits)
      else
        a = c
        c = d
        fc = fd
        d = a + phi * (b - a)
        scaled = exp(d) * llrs
        fd = -gmi_of_llrs(scaled, true_bits)
      end if
    end do
    deallocate(scaled)
    beta_star = exp(0.5_dp * (a + b))
  end function optimal_temperature

  !> Evaluate log(1 + exp(-x)) without overflow.
  !!
  !! @param[in] x Real-valued argument.
  !! @return Stable logistic-loss term.
  pure real(dp) function log1p_exp_neg(x) result(y)
    real(dp), intent(in) :: x
    if (x > 30.0_dp) then
      y = exp(-x)
    else if (x < -30.0_dp) then
      y = -x
    else
      y = log(1.0_dp + exp(-x))
    end if
  end function log1p_exp_neg
end module cdss_calibration
