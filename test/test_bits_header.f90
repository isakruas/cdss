!> @file test_bits_header.f90
!! @brief Unit tests for bit conversions, CRCs, headers, and frame sizing.
program test_bits_header
  use cdss_constants, only: frame_header_bits
  use cdss_bits, only: int_to_bits, bits_to_int, bits_to_string, bits_from_string, &
                         crc16_ccitt, payload_crc16
  use cdss_receiver, only: build_header, decode_header, frame_n_bits, frame_n_samples
  use cdss_types, only: waveform_config, coding_config
  implicit none

  integer :: byte_bits(8), value, crc, expected_crc, length, decoded_crc
  integer, allocatable :: parsed(:), payload(:), payload_packed(:)
  integer :: header(frame_header_bits)
  type(waveform_config) :: waveform
  type(coding_config) :: coding
  character(len=:), allocatable :: text

  print '(a)', "test_bits_header: starting bit/CRC/header checks"

  call int_to_bits(int(z'A5'), 8, byte_bits)
  value = bits_to_int(byte_bits)
  print '(a,i0,a,z2.2)', "  int_to_bits/bits_to_int roundtrip: value=", value, ", hex=", value
  if (value /= int(z'A5')) then
    print '(a)', "test_bits_header: integer bit roundtrip failed"
    error stop 1
  end if

  call bits_from_string("101x001", parsed)
  text = bits_to_string(parsed)
  print '(a,a)', "  parsed string=", text
  if (text /= "1010001") then
    print '(a)', "test_bits_header: string parsing/serialization failed"
    print '(a,a)', "  expected=1010001, actual=", text
    error stop 2
  end if

  crc = crc16_ccitt([ichar('1'), ichar('2'), ichar('3'), ichar('4'), ichar('5'), &
                     ichar('6'), ichar('7'), ichar('8'), ichar('9')])
  expected_crc = int(z'29B1')
  print '(a,z4.4,a,z4.4)', "  CRC-16/CCITT test vector: actual=0x", crc, &
       ", expected=0x", expected_crc
  if (crc /= expected_crc) then
    print '(a)', "test_bits_header: CRC-16/CCITT known-vector check failed"
    error stop 3
  end if

  allocate(payload_packed(16))
  call int_to_bits(int(z'ABCD'), 16, payload_packed)
  crc = payload_crc16(payload_packed)
  expected_crc = crc16_ccitt([int(z'AB'), int(z'CD')])
  print '(a,z4.4,a,z4.4)', "  payload CRC: actual=0x", crc, ", expected=0x", expected_crc
  if (crc /= expected_crc) then
    print '(a)', "test_bits_header: payload CRC packing check failed"
    error stop 4
  end if

  allocate(payload(70))
  payload = [(mod(value, 2), value = 1, 70)]
  call build_header(payload, header)
  call decode_header(header, length, decoded_crc)
  print '(a,i0,a,z4.4)', "  decoded header: length=", length, ", crc=0x", decoded_crc
  if (length /= size(payload)) then
    print '(a)', "test_bits_header: header payload length mismatch"
    error stop 5
  end if
  if (decoded_crc /= payload_crc16(payload)) then
    print '(a)', "test_bits_header: header CRC mismatch"
    error stop 6
  end if

  print '(a,i0)', "  frame coded bits for 70 payload bits=", frame_n_bits(coding, 70)
  if (frame_n_bits(coding, 70) /= coding%header_n + 2 * coding%payload_n) then
    print '(a)', "test_bits_header: frame_n_bits should allocate two payload blocks"
    error stop 7
  end if
  if (frame_n_samples(waveform, coding, 70) /= &
      (waveform%preamble_bits + frame_n_bits(coding, 70)) * waveform%samples_per_bit()) then
    print '(a)', "test_bits_header: frame_n_samples geometry check failed"
    error stop 8
  end if

  print '(a)', "test_bits_header: passed"
  deallocate(parsed, payload, payload_packed)
end program test_bits_header
