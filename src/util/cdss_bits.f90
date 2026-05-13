!> @file cdss_bits.f90
!! @brief Bit-vector conversion, text serialization, and CRC utilities.
!!
!! @details
!! The modem represents all encoded and uncoded binary data as integer arrays
!! containing only 0 and 1. This module provides deterministic conversions
!! between integers, strings, packed bytes, and BPSK sign conventions.
module cdss_bits
  use cdss_kinds, only: dp, i32
  implicit none
  private
  public :: int_to_bits, bits_to_int
  public :: payload_crc16, crc16_ccitt
  public :: bits_to_string, bits_from_string
  public :: bits_to_npolar_sign

contains

  !> Convert an integer to a fixed-width big-endian bit vector.
  !!
  !! @param[in]  value Integer value to unpack.
  !! @param[in]  width Number of output bits.
  !! @param[out] bits  Output array; element 1 receives the most significant bit.
  subroutine int_to_bits(value, width, bits)
    integer, intent(in) :: value, width
    integer, intent(out) :: bits(:)
    integer :: i
    do i = 1, width
      bits(i) = iand(ishft(value, -(width - i)), 1)
    end do
  end subroutine int_to_bits

  !> Convert a big-endian bit vector to an integer.
  !!
  !! @param[in] bits Input bits with the most significant bit at element 1.
  !! @return Integer represented by the input vector.
  pure integer function bits_to_int(bits) result(value)
    integer, intent(in) :: bits(:)
    integer :: i
    value = 0
    do i = 1, size(bits)
      value = ishft(value, 1) + iand(bits(i), 1)
    end do
  end function bits_to_int

  !> Compute the CRC-16/CCITT-FALSE checksum over packed bytes.
  !!
  !! @param[in] bytes Integer byte values; only the least significant 8 bits are used.
  !! @return 16-bit CRC value with polynomial 0x1021 and initial state 0xFFFF.
  pure integer function crc16_ccitt(bytes) result(crc)
    integer, intent(in) :: bytes(:)
    integer :: i, j, c
    crc = int(z'FFFF')
    do i = 1, size(bytes)
      c = iand(bytes(i), int(z'FF'))
      crc = ieor(crc, ishft(c, 8))
      do j = 1, 8
        if (iand(crc, int(z'8000')) /= 0) then
          crc = iand(ieor(ishft(crc, 1), int(z'1021')), int(z'FFFF'))
        else
          crc = iand(ishft(crc, 1), int(z'FFFF'))
        end if
      end do
    end do
  end function crc16_ccitt

  !> Compute the payload CRC over an unpacked bit vector.
  !!
  !! Bits are packed most-significant-bit first into bytes. The final byte is
  !! zero-padded when the bit count is not a multiple of eight.
  !!
  !! @param[in] bits Payload bits represented as 0 or 1.
  !! @return CRC-16/CCITT-FALSE checksum.
  pure integer function payload_crc16(bits) result(crc)
    integer, intent(in) :: bits(:)
    integer :: n_bytes, i, j, k, acc
    integer, allocatable :: bytes(:)
    n_bytes = (size(bits) + 7) / 8
    allocate(bytes(n_bytes))
    bytes = 0
    do i = 1, size(bits)
      j = (i - 1) / 8 + 1
      k = mod(i - 1, 8)
      acc = ior(bytes(j), ishft(iand(bits(i), 1), 7 - k))
      bytes(j) = acc
    end do
    crc = crc16_ccitt(bytes)
    deallocate(bytes)
  end function payload_crc16

  !> Serialize a bit vector as a character string of '0' and '1'.
  !!
  !! @param[in] bits Input bit vector.
  !! @return Allocated character string with one character per bit.
  function bits_to_string(bits) result(out)
    integer, intent(in) :: bits(:)
    character(len=:), allocatable :: out
    integer :: i
    allocate(character(len=size(bits)) :: out)
    do i = 1, size(bits)
      if (iand(bits(i), 1) == 1) then
        out(i:i) = '1'
      else
        out(i:i) = '0'
      end if
    end do
  end function bits_to_string

  !> Parse a character string into a bit vector.
  !!
  !! @param[in]  text Input text. Characters equal to '1' become one; all others
  !!                  become zero.
  !! @param[out] bits Allocated output bit vector.
  subroutine bits_from_string(text, bits)
    character(len=*), intent(in) :: text
    integer, allocatable, intent(out) :: bits(:)
    integer :: i
    allocate(bits(len_trim(text)))
    do i = 1, len_trim(text)
      bits(i) = merge(1, 0, text(i:i) == '1')
    end do
  end subroutine bits_from_string

  !> Map binary bits to the non-polar BPSK sign convention.
  !!
  !! @param[in]  bits  Input bits represented as 0 or 1.
  !! @param[out] signs Output signs where bit 1 maps to +1 and bit 0 maps to -1.
  pure subroutine bits_to_npolar_sign(bits, signs)
    integer,  intent(in)  :: bits(:)
    real(dp), intent(out) :: signs(:)
    integer :: i
    do i = 1, size(bits)
      if (iand(bits(i), 1) == 1) then
        signs(i) = 1.0_dp
      else
        signs(i) = -1.0_dp
      end if
    end do
  end subroutine bits_to_npolar_sign
end module cdss_bits
