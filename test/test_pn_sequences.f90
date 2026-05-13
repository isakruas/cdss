!> @file test_pn_sequences.f90
!! @brief Unit test for m-sequence balance and Gold-code diversity.
program test_pn_sequences
  use cdss_kinds, only: dp
  use cdss_pn,    only: pn_gold, pn_msequence
  implicit none
  integer, allocatable :: a(:), b(:)
  integer :: ones, zeros, period, length, distance

  ! A full-period n=6 m-sequence should be balanced: 32 ones and 31 zeros.
  period = 2**6 - 1
  length = period
  print '(a)', "test_pn_sequences: starting PN sequence checks"
  print '(a,i0,a,i0)', "  m-sequence stages=6, period=", period
  call pn_msequence(6, [6, 5], 1, length, a)
  ones  = count(a == 1)
  zeros = count(a == 0)
  print '(a,i0,a,i0,a,i0)', &
       "  m-sequence counts: ones=", ones, ", zeros=", zeros, ", length=", length
  if (abs(ones - 32) > 1) then
    print '(a)', "test_pn_sequences: m-sequence balance check failed"
    print '(a,i0,a,i0)', "  expected approximately 32 ones, actual ones=", ones
    error stop 1
  end if
  if (zeros + ones /= length) then
    print '(a)', "test_pn_sequences: sequence contains values outside {0,1}"
    print '(a,i0,a,i0)', "  counted binary chips=", zeros + ones, ", length=", length
    error stop 2
  end if

  ! Distinct Gold-code offsets should not collapse to nearly identical chips.
  call pn_gold(6, length, 0, a)
  call pn_gold(6, length, 7, b)
  distance = count(a /= b)
  print '(a,i0,a,i0)', "  Gold-code Hamming distance=", distance, ", threshold=", length / 4
  if (distance < length / 4) then
    print '(a)', "test_pn_sequences: Gold-code diversity check failed"
    print '(a,i0,a,i0)', "  distance=", distance, ", minimum required=", length / 4
    error stop 3
  end if

  print '(a)', "test_pn_sequences: passed"
  deallocate(a, b)
end program test_pn_sequences
