!> @file cdss_chirp_pulse.f90
!! @brief Chirp carrier and chip-envelope synthesis.
!!
!! @details
!! Each coded bit is transmitted as a PN-spread chip sequence over a slow linear
!! chirp. The complex carrier for one bit is
!!
!!     s(t) = exp(j*2*pi*(fc*t + 0.5*mu*t**2)),
!!
!! where t is centered on the bit interval, fc is the audio carrier, and mu is
!! the chirp slope. Chip edges are optionally tapered by a raised-cosine envelope
!! to reduce spectral splatter.
module cdss_chirp_pulse
  use cdss_kinds, only: dp
  use cdss_constants, only: pi, two_pi
  use cdss_types, only: waveform_config
  implicit none
  private
  public :: build_chirp_carrier, build_chip_envelope

contains

  !> Generate the complex chirp carrier for one coded bit.
  !!
  !! @param[in]  cfg     Waveform configuration.
  !! @param[out] carrier Allocated complex carrier with cfg%samples_per_bit() samples.
  subroutine build_chirp_carrier(cfg, carrier)
    type(waveform_config), intent(in) :: cfg
    complex(dp), allocatable, intent(out) :: carrier(:)
    integer :: n_samples, i
    real(dp) :: t, phase, dt, half_dur
    n_samples = cfg%samples_per_bit()
    allocate(carrier(n_samples))
    dt = 1.0_dp / cfg%sample_rate
    half_dur = 0.5_dp * real(n_samples, dp) * dt
    do i = 1, n_samples
      t = (real(i - 1, dp) * dt) - half_dur
      phase = two_pi * (cfg%carrier_hz * t + 0.5_dp * cfg%chirp_slope_hz_s * t * t)
      carrier(i) = cmplx(cos(phase), sin(phase), dp)
    end do
  end subroutine build_chirp_carrier

  !> Generate a raised-cosine amplitude envelope for one PN chip.
  !!
  !! @param[in]  samples_per_chip Number of samples in one chip.
  !! @param[in]  taper_fraction   Fraction of the chip used for each edge ramp.
  !! @param[out] envelope         Allocated real envelope samples.
  subroutine build_chip_envelope(samples_per_chip, taper_fraction, envelope)
    integer, intent(in)  :: samples_per_chip
    real(dp), intent(in) :: taper_fraction
    real(dp), allocatable, intent(out) :: envelope(:)
    integer :: edge, i
    allocate(envelope(samples_per_chip))
    envelope = 1.0_dp
    if (taper_fraction <= 0.0_dp) return
    edge = max(1, nint(taper_fraction * real(samples_per_chip, dp)))
    do i = 1, edge
      envelope(i) = 0.5_dp - 0.5_dp * cos(pi * real(i - 1, dp) / real(edge, dp))
      envelope(samples_per_chip - i + 1) = envelope(i)
    end do
  end subroutine build_chip_envelope
end module cdss_chirp_pulse
