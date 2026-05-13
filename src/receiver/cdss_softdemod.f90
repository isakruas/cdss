!> @file cdss_softdemod.f90
!! @brief Soft-decision demodulator for chirp-DSSS coded bits.
!!
!! @details
!! The demodulator converts aligned complex baseband samples into one LLR per
!! coded bit. For each bit it dechirps the waveform, correlates each chip with
!! the PN sign and chip envelope, weights chips by the STFT-derived exposure
!! estimate, and tracks residual carrier phase with a decision-directed Costas
!! loop. Positive output LLRs favor coded bit 0 to match the Polar decoder
!! convention used by cdss_polar.
module cdss_softdemod
  use cdss_kinds, only: dp
  use cdss_constants, only: eps
  use cdss_types, only: waveform_config, soft_metric_config, soft_demod_result
  use cdss_waveform, only: get_pn_signs, get_carrier, get_envelope
  use cdss_uncertainty, only: compute_chip_uncertainty
  implicit none
  private
  public :: soft_demod_bits, soft_demod_noise_estimate

contains

  !> Estimate baseline complex-sample noise variance from adjacent differences.
  !!
  !! @param[in] rx Complex baseband samples.
  !! @return Estimated noise variance per complex sample.
  !!
  !! The first-difference estimator is appropriate for slowly varying,
  !! oversampled chip waveforms because adjacent signal samples are highly
  !! correlated while independent white-noise samples are not.
  pure function soft_demod_noise_estimate(rx) result(var)
    complex(dp), intent(in) :: rx(:)
    real(dp) :: var
    integer :: i
    real(dp) :: acc
    if (size(rx) < 2) then
      var = 1.0_dp
      return
    end if
    acc = 0.0_dp
    do i = 1, size(rx) - 1
      acc = acc + abs(rx(i + 1) - rx(i))**2
    end do
    var = max(acc / real(2 * (size(rx) - 1), dp), eps)
  end function soft_demod_noise_estimate

  !> Demodulate aligned samples into soft bit decisions.
  !!
  !! @param[in]  cfg        Waveform configuration.
  !! @param[in]  metric     Soft-metric and interference-weighting parameters.
  !! @param[in]  rx         Received samples aligned so rx(1) starts bit 1.
  !! @param[in]  n_bits     Number of coded bits to demodulate.
  !! @param[out] result     LLRs, hard decisions, chip exposure, and weights.
  !! @param[in]  bg_residue Optional residue used only for STFT exposure mapping.
  subroutine soft_demod_bits(cfg, metric, rx, n_bits, result, bg_residue)
    type(waveform_config),    intent(in)  :: cfg
    type(soft_metric_config), intent(in)  :: metric
    complex(dp),              intent(in)  :: rx(:)
    integer,                  intent(in)  :: n_bits
    type(soft_demod_result),  intent(out) :: result
    complex(dp), intent(in), optional     :: bg_residue(:)

    real(dp), allocatable :: pn_signs(:), envelope(:)
    complex(dp), allocatable :: carrier(:)
    real(dp), allocatable :: chip_J(:), chip_w(:)
    integer :: samples_per_chip, samples_per_bit, b_idx, c_idx, k, base, abs_chip
    real(dp) :: bit_sum, env_norm, env_energy
    real(dp) :: noise_var, chip_noise_var, llr
    integer :: total_samples, total_chips

    ! Costas-loop state variables.
    complex(dp) :: chip_corr_cplx, bit_sum_cplx
    real(dp) :: phase_est, freq_est, phase_err, alpha, beta

    samples_per_chip = cfg%samples_per_chip()
    samples_per_bit  = cfg%samples_per_bit()
    total_samples    = n_bits * samples_per_bit
    total_chips      = n_bits * cfg%spreading_factor

    call get_pn_signs(cfg, pn_signs)
    call get_carrier(cfg, carrier)
    call get_envelope(cfg, envelope)
    ! Envelope normalization defines the coherent chip estimator. The matched
    ! real-component noise variance follows from the complex-sample variance.
    env_norm = sum(envelope)
    if (env_norm < eps) env_norm = eps
    env_energy = sum(envelope**2)

    ! Slice rx to the captured payload region (clip to availability).
    if (size(rx) < total_samples) then
      error stop "soft_demod_bits: rx shorter than n_bits * samples_per_bit"
    end if

    if (metric%has_noise_var) then
      noise_var = metric%noise_var
    else
      noise_var = soft_demod_noise_estimate(rx(1:total_samples))
    end if

    ! Estimate per-chip exposure J_chip over the whole captured run.
    ! If a residue is provided, use it to map the true interference floor.
    if (present(bg_residue)) then
      call compute_chip_uncertainty(cfg, metric, bg_residue(1:total_samples), chip_J)
    else
      call compute_chip_uncertainty(cfg, metric, rx(1:total_samples), chip_J)
    end if

    chip_noise_var = max(noise_var * env_energy / (2.0_dp * env_norm * env_norm), eps)
    allocate(chip_w(total_chips))
    do k = 1, total_chips
      chip_w(k) = 1.0_dp / (chip_noise_var * (1.0_dp + metric%exposure_weight * chip_J(k)))
    end do

    allocate(result%bit_llrs(n_bits), result%hard_bits(n_bits))
    allocate(result%chip_J(total_chips), result%chip_weight(total_chips))
    result%chip_J      = chip_J
    result%chip_weight = chip_w
    result%noise_var   = noise_var

    ! The second-order Costas loop tracks constant phase error and slow residual
    ! frequency error without pilot symbols. Squaring in the discriminator below
    ! removes the 180-degree BPSK ambiguity.
    alpha = 0.05_dp      ! Proportional gain: dictates immediate phase correction speed
    beta  = 0.001_dp     ! Integral gain: tracks long-term frequency drift (CFO)
    phase_est = 0.0_dp
    freq_est  = 0.0_dp

    do b_idx = 1, n_bits
      base = (b_idx - 1) * samples_per_bit
      bit_sum_cplx = (0.0_dp, 0.0_dp)
      do c_idx = 1, cfg%spreading_factor
        abs_chip = (b_idx - 1) * cfg%spreading_factor + c_idx
        chip_corr_cplx = (0.0_dp, 0.0_dp)
        do k = 1, samples_per_chip
          ! Matched dechirping converts the chirp-spread chip to a BPSK sample.
          chip_corr_cplx = chip_corr_cplx + envelope(k) * &
               ( rx(base + (c_idx - 1) * samples_per_chip + k) * &
                 conjg(carrier((c_idx - 1) * samples_per_chip + k)) )
        end do
        chip_corr_cplx = chip_corr_cplx / env_norm
        bit_sum_cplx = bit_sum_cplx + chip_w(abs_chip) * pn_signs(c_idx) * chip_corr_cplx
      end do

      ! Derotate the accumulated bit using the current loop phase estimate.
      bit_sum_cplx = bit_sum_cplx * exp(cmplx(0.0_dp, -phase_est, dp))

      ! Squaring collapses the BPSK constellation to one phase state, allowing
      ! residual phase error estimation independent of the unknown data bit.
      phase_err = 0.5_dp * atan2(aimag(bit_sum_cplx**2), real(bit_sum_cplx**2, dp))

      ! Update proportional and integral loop states for the next bit.
      phase_est = phase_est + alpha * phase_err
      freq_est  = freq_est  + beta  * phase_err
      phase_est = phase_est + freq_est

      ! Extract coherent real part for the LLR. The waveform maps coded bit 0 to
      ! -1 and coded bit 1 to +1, whereas the Polar decoder expects positive LLRs
      ! to favor bit 0; therefore the matched BPSK statistic is negated.
      bit_sum = real(bit_sum_cplx, dp)
      llr   = -metric%temperature * 2.0_dp * bit_sum

      llr = max(-metric%llr_clip, min(metric%llr_clip, llr))
      result%bit_llrs(b_idx) = llr
      result%hard_bits(b_idx) = merge(0, 1, llr >= 0.0_dp)  ! Positive LLR favors bit 0.
    end do

    result%temperature_used = metric%temperature

    deallocate(pn_signs, envelope, carrier, chip_J, chip_w)
  end subroutine soft_demod_bits
end module cdss_softdemod
