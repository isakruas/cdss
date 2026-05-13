!> @file test_loopback.f90
!! @brief End-to-end loopback test over a noiseless channel.
program test_loopback
  use cdss_kinds, only: dp
  use cdss_api,   only: cdss_modem
  use cdss_types, only: frame_recovery_result
  implicit none

  type(cdss_modem) :: modem
  type(frame_recovery_result) :: frame
  integer, allocatable :: payload(:)
  complex(dp), allocatable :: tx(:)
  integer :: i, mismatch

  print '(a)', "test_loopback: starting noiseless end-to-end modem check"
  allocate(payload(64))
  do i = 1, 64
    payload(i) = mod(i, 2)
  end do
  print '(a,i0,a,i0)', "  payload bits=", size(payload), ", ones=", count(payload == 1)
  print '(a,f0.1,a,i0)', "  sample rate=", modem%waveform%sample_rate, &
       " Hz, spreading factor=", modem%waveform%spreading_factor

  call modem%build_frame(payload, tx)
  print '(a,i0)', "  synthesized samples=", size(tx)

  call modem%recover_frame(tx, 256, frame)
  print '(a,l1)', "  crc_ok=", frame%crc_ok
  print '(a,i0,a,i0)', "  CRC expected=", frame%payload_crc_expected, &
       ", observed=", frame%payload_crc_observed
  print '(a,i0,a,i0)', "  decoded payload bits=", size(frame%payload_bits), &
       ", reported count=", frame%payload_bit_count
  print '(a,i0)', "  refinement iterations used=", frame%refinement_iters_used

  if (.not. frame%crc_ok) then
    print '(a)', "test_loopback: CRC failed (noiseless decode should always pass)"
    error stop 1
  end if
  if (size(frame%payload_bits) /= size(payload)) then
    print '(a)', "test_loopback: decoded payload length mismatch"
    print '(a,i0,a,i0)', "  expected=", size(payload), ", actual=", size(frame%payload_bits)
    error stop 2
  end if
  if (any(frame%payload_bits /= payload)) then
    mismatch = first_mismatch(payload, frame%payload_bits)
    print '(a)', "test_loopback: decoded payload differs from transmitted payload"
    print '(a,i0)', "  first mismatch index=", mismatch
    if (mismatch > 0) then
      print '(a,i0,a,i0)', "  expected=", payload(mismatch), &
           ", actual=", frame%payload_bits(mismatch)
    end if
    print '(a,i0)', "  total mismatches=", count(frame%payload_bits /= payload)
    error stop 3
  end if

  print '(a)', "test_loopback: passed (end-to-end TX/RX over noiseless channel)"
  deallocate(payload, tx)

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
end program test_loopback
