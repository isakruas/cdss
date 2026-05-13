!> @file cdss_constants.f90
!! @brief Mathematical constants and default modem parameters.
!! @details
!! The values in this module define the default numerical assumptions shared by
!! modulation, channel simulation, soft demodulation, and channel coding. Keeping
!! them in one module makes simulation results reproducible and simplifies peer
!! review of the experimental configuration.

module cdss_constants
  use cdss_kinds, only: dp
  implicit none
  private
  public :: pi, two_pi, eps, sqrt2, ln2
  public :: frame_header_bits, max_payload_bits, default_polar_list
  public :: default_polar_n, default_polar_k
  public :: header_interleaver_domain, payload_interleaver_domain, default_interleaver_seed

  real(dp), parameter :: pi      = acos(-1.0_dp)  !< Archimedes' constant.
  real(dp), parameter :: two_pi  = 2.0_dp * pi    !< Radians in one full cycle.
  real(dp), parameter :: eps     = 1.0e-12_dp     !< Lower bound for divisions and variances.
  real(dp), parameter :: sqrt2   = sqrt(2.0_dp)   !< Square root of two.
  real(dp), parameter :: ln2     = log(2.0_dp)    !< Natural logarithm of two.

  !> Number of uncoded header bits: 32-bit payload length plus 16-bit CRC.
  integer, parameter :: frame_header_bits = 48
  !> Absolute payload-size guard in bits used to reject malformed headers.
  integer, parameter :: max_payload_bits  = 1000000

  !> Default short rate-1/2 Polar code used by the prototype.
  integer, parameter :: default_polar_n    = 128  !< Polar codeword length.
  integer, parameter :: default_polar_k    = 64   !< Polar information length.
  integer, parameter :: default_polar_list = 8    !< SCL decoder list size.

  !> Domain-separated interleaver seeds for independent header and payload use.
  integer, parameter :: default_interleaver_seed   = int(z'4344535A')   !< Base seed, "CDSZ".
  integer, parameter :: header_interleaver_domain  = int(z'48445231')   !< Header domain, "HDR1".
  integer, parameter :: payload_interleaver_domain = int(z'504C4431')   !< Payload domain, "PLD1".
end module cdss_constants
