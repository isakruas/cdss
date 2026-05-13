!> @file cdss_kinds.f90
!! @brief Numeric kind and C-interoperability definitions for the CDSS modem.
!! @details
!! Centralizing kind parameters prevents precision drift across numerical
!! kernels, file I/O, FFTW bindings, and public APIs. The exported aliases are
!! intentionally short because they appear throughout mathematical kernels.

module cdss_kinds
  use, intrinsic :: iso_fortran_env, only: real64, real32, int32, int64, int8
  use, intrinsic :: iso_c_binding,  only: c_int, c_double, c_double_complex, c_size_t
  implicit none
  private
  public :: dp, sp, i32, i64, i8
  public :: c_int, c_double, c_double_complex, c_size_t

  integer, parameter :: dp  = real64  !< IEEE binary64 real kind.
  integer, parameter :: sp  = real32  !< IEEE binary32 real kind.
  integer, parameter :: i32 = int32   !< 32-bit signed integer kind.
  integer, parameter :: i64 = int64   !< 64-bit signed integer kind.
  integer, parameter :: i8  = int8    !< 8-bit signed integer kind.
end module cdss_kinds
