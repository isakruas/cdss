!> @file test_calibration.f90
!! @brief Unit test for GMI ordering and temperature-search sanity.
program test_calibration
  use cdss_kinds, only: dp
  use cdss_calibration, only: gmi_of_llrs, optimal_temperature
  implicit none
  integer, allocatable :: bits(:)
  real(dp), allocatable :: llrs(:), inverted(:)
  real(dp) :: gmi_good, gmi_bad, beta
  integer :: i, n

  n = 256
  print '(a)', "test_calibration: starting GMI and temperature checks"
  print '(a,i0)', "  vector length=", n
  allocate(bits(n), llrs(n), inverted(n))
  do i = 1, n
    bits(i) = mod(i, 2)
    if (bits(i) == 0) then
      llrs(i) =  3.0_dp                ! Confident bit 0 under the LLR convention.
    else
      llrs(i) = -3.0_dp
    end if
  end do
  inverted = -llrs

  gmi_good = gmi_of_llrs(llrs, bits)
  gmi_bad  = gmi_of_llrs(inverted, bits)
  print '(a,es12.4)', "  GMI with correctly signed LLRs=", gmi_good
  print '(a,es12.4)', "  GMI with inverted LLRs=", gmi_bad
  if (gmi_good <= gmi_bad) then
    print '(a)', "test_calibration: failed GMI ordering check"
    print '(a)', "  expected correctly signed LLRs to have higher GMI"
    error stop 1
  end if

  beta = optimal_temperature(llrs, bits)
  print '(a,es12.4)', "  optimal temperature beta=", beta
  if (beta <= 0.0_dp) then
    print '(a)', "test_calibration: beta must be positive"
    error stop 2
  end if
  if (beta > 100.0_dp) then
    print '(a)', "test_calibration: beta is unexpectedly large"
    print '(a,es12.4)', "  beta=", beta
    error stop 3
  end if

  print '(a)', "test_calibration: passed"
  deallocate(bits, llrs, inverted)
end program test_calibration
