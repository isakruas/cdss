!> @file cdss_waveform.f90
!! @brief Transmit waveform assembly for coded CDSS bit streams.
!!
!! @details
!! Given coded bits in {0,1}, this module produces the complex samples presented
!! to the channel model. Each coded bit is mapped to a BPSK sign, multiplied by a
!! deterministic PN chip sequence, shaped by a chip envelope, and placed on the
!! per-bit chirp carrier. The receiver rebuilds the same reference components for
!! matched dechirping and despreading.
module cdss_waveform
  use cdss_kinds, only: dp
  use cdss_constants, only: pi, two_pi
  use cdss_types, only: waveform_config
  use cdss_pn, only: pn_gold, pn_to_npolar
  use cdss_chirp_pulse, only: build_chirp_carrier, build_chip_envelope
  implicit none
  private
  public :: synthesize_waveform, build_preamble_pattern
  public :: get_pn_signs, get_carrier, get_envelope

contains

  !> Synthesize a contiguous complex baseband waveform from coded bits.
  !!
  !! @param[in]  cfg  Waveform configuration.
  !! @param[in]  bits Coded bits represented as 0 or 1.
  !! @param[out] iq   Allocated complex waveform samples.
  subroutine synthesize_waveform(cfg, bits, iq)
    type(waveform_config), intent(in) :: cfg
    integer, intent(in) :: bits(:)
    complex(dp), allocatable, intent(out) :: iq(:)

    real(dp), allocatable :: pn_signs(:)
    real(dp), allocatable :: envelope(:)
    complex(dp), allocatable :: carrier(:)
    integer :: n_bits, samples_per_bit, samples_per_chip, total_samples
    integer :: b_idx, c_idx, k, base
    real(dp) :: bit_sign, chip_sign

    n_bits = size(bits)
    samples_per_chip = cfg%samples_per_chip()
    samples_per_bit  = cfg%samples_per_bit()
    total_samples    = n_bits * samples_per_bit

    call get_pn_signs(cfg, pn_signs)
    call build_chirp_carrier(cfg, carrier)
    call build_chip_envelope(samples_per_chip, cfg%pulse_taper, envelope)

    allocate(iq(total_samples))
    iq = (0.0_dp, 0.0_dp)

    do b_idx = 1, n_bits
      bit_sign = merge(1.0_dp, -1.0_dp, iand(bits(b_idx), 1) == 1)
      base = (b_idx - 1) * samples_per_bit
      do c_idx = 1, cfg%spreading_factor
        chip_sign = bit_sign * pn_signs(c_idx)
        do k = 1, samples_per_chip
          iq(base + (c_idx - 1) * samples_per_chip + k) = &
               cmplx(chip_sign * envelope(k), 0.0_dp, dp) * &
               carrier((c_idx - 1) * samples_per_chip + k)
        end do
      end do
    end do
    deallocate(pn_signs, carrier, envelope)
  end subroutine synthesize_waveform

  !> Build a deterministic alternating-bit preamble pattern.
  !!
  !! @param[in]  cfg     Waveform configuration.
  !! @param[out] pattern Allocated preamble bits.
  !!
    !! @note The frame encoder uses a configurable m-sequence preamble. This
    !! routine is retained as a deterministic helper for experiments that require
    !! an alternating preamble.
  subroutine build_preamble_pattern(cfg, pattern)
    type(waveform_config), intent(in) :: cfg
    integer, allocatable, intent(out) :: pattern(:)
    integer :: i
    allocate(pattern(cfg%preamble_bits))
    do i = 1, cfg%preamble_bits
      pattern(i) = mod(i + 1, 2)               ! 1,0,1,0,... for deterministic tests.
    end do
  end subroutine build_preamble_pattern

  !> Return the PN spreading signs for the configured waveform.
  !!
  !! @param[in]  cfg   Waveform configuration.
  !! @param[out] signs Allocated signs with length cfg%spreading_factor.
  subroutine get_pn_signs(cfg, signs)
    type(waveform_config), intent(in) :: cfg
    real(dp), allocatable, intent(out) :: signs(:)
    integer, allocatable :: pn_bits(:)
    call pn_gold(cfg%pn_length_bits, cfg%spreading_factor, cfg%pn_seed, pn_bits)
    allocate(signs(size(pn_bits)))
    call pn_to_npolar(pn_bits, signs)
    deallocate(pn_bits)
  end subroutine get_pn_signs

  !> Return the matched chirp carrier for one coded bit.
  !!
  !! @param[in]  cfg     Waveform configuration.
  !! @param[out] carrier Allocated complex carrier samples.
  subroutine get_carrier(cfg, carrier)
    type(waveform_config), intent(in) :: cfg
    complex(dp), allocatable, intent(out) :: carrier(:)
    call build_chirp_carrier(cfg, carrier)
  end subroutine get_carrier

  !> Return the matched chip envelope for the configured waveform.
  !!
  !! @param[in]  cfg      Waveform configuration.
  !! @param[out] envelope Allocated chip-envelope samples.
  subroutine get_envelope(cfg, envelope)
    type(waveform_config), intent(in) :: cfg
    real(dp), allocatable, intent(out) :: envelope(:)
    call build_chip_envelope(cfg%samples_per_chip(), cfg%pulse_taper, envelope)
  end subroutine get_envelope
end module cdss_waveform
