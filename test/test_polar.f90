!> @file test_polar.f90
!! @brief Unit test for noiseless Polar encode/decode round-trip behavior.
program test_polar
  use cdss_kinds, only: dp
  use cdss_polar, only: polar_code, polar_encode_stream, polar_decode_stream
  implicit none

  type(polar_code) :: code
  integer :: i, mismatch
  integer, allocatable :: info(:), coded(:), info_back(:)
  real(dp), allocatable :: llrs(:)

  code%n = 128
  code%k = 64
  code%list_size = 4
  print '(a)', "test_polar: starting noiseless Polar round-trip"
  print '(a,i0,a,i0,a,i0)', "  code: N=", code%n, ", K=", code%k, ", list=", code%list_size

  allocate(info(code%k))
  do i = 1, code%k
    info(i) = mod(i, 2)
  end do
  print '(a,i0,a,i0)', "  input bits: ", size(info), ", ones=", count(info == 1)

  call polar_encode_stream(code, info, coded)
  print '(a,i0,a,i0)', "  encoded bits: ", size(coded), ", ones=", count(coded == 1)
  if (size(coded) /= code%n) then
    print '(a)', "test_polar: failed size check after encoding"
    print '(a,i0,a,i0)', "  expected coded size=", code%n, ", actual=", size(coded)
    error stop 1
  end if

  ! Convert coded bits to ideal high-magnitude LLRs using the decoder sign
  ! convention: positive values favor bit 0 and negative values favor bit 1.
  allocate(llrs(code%n))
  do i = 1, code%n
    if (coded(i) == 0) then
      llrs(i) =  10.0_dp
    else
      llrs(i) = -10.0_dp
    end if
  end do
  print '(a,es12.4,a,es12.4)', "  LLR range: min=", minval(llrs), ", max=", maxval(llrs)

  call polar_decode_stream(code, llrs, info_back)
  print '(a,i0,a,i0)', "  decoded bits: ", size(info_back), ", ones=", count(info_back == 1)
  if (size(info_back) /= code%k) then
    print '(a)', "test_polar: failed decoded-size check"
    print '(a,i0,a,i0)', "  expected decoded size=", code%k, ", actual=", size(info_back)
    error stop 2
  end if
  if (any(info_back /= info)) then
    mismatch = first_mismatch(info, info_back)
    print '(a)', "test_polar: decoded bits differ from input bits"
    print '(a,i0)', "  first mismatch index=", mismatch
    if (mismatch > 0) then
      print '(a,i0,a,i0)', "  expected=", info(mismatch), ", actual=", info_back(mismatch)
    end if
    print '(a,i0)', "  total mismatches=", count(info_back /= info)
    error stop 3
  end if

  print '(a)', "test_polar: passed (encode/decode roundtrip on noiseless LLRs)"
  deallocate(info, coded, info_back, llrs)

contains

  integer function first_mismatch(expected, actual) result(idx)
    integer, intent(in) :: expected(:), actual(:)
    integer :: j, n
    idx = 0
    n = min(size(expected), size(actual))
    do j = 1, n
      if (expected(j) /= actual(j)) then
        idx = j
        return
      end if
    end do
  end function first_mismatch
end program test_polar
