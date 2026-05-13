!> @file cdss_channel.f90
!! @brief Channel-impairment models for controlled CDSS simulations.
!!
!! @details
!! The routines in this module apply additive white Gaussian noise, carrier
!! offset, linear frequency drift, continuous-wave interferers, sample-clock
!! error, and guard-padding. They are intentionally deterministic when driven by
!! a seeded rng_state so simulation results can be reproduced for peer review.
module cdss_channel
  use cdss_kinds, only: dp
  use cdss_constants, only: pi, two_pi, eps
  use cdss_types, only: waveform_config, channel_config, tone_interferer
  use cdss_rng, only: rng_state, rng_normal
  implicit none
  private
  public :: apply_awgn_channel, apply_cfo_and_drift, apply_tone_interferers
  public :: apply_sample_clock_error, signal_power, simulate_full_channel

contains

  !> Compute the mean complex-sample power of a signal.
  !!
  !! @param[in] signal Complex baseband samples.
  !! @return Average value of |signal|**2, or zero for an empty array.
  pure function signal_power(signal) result(p)
    complex(dp), intent(in) :: signal(:)
    real(dp) :: p
    if (size(signal) == 0) then
      p = 0.0_dp
    else
      p = sum(abs(signal)**2) / real(size(signal), dp)
    end if
  end function signal_power

  !> Add complex AWGN at a requested signal-to-noise ratio.
  !!
  !! @param[in]    input             Input complex samples.
  !! @param[in]    snr_db            Target SNR in decibels.
  !! @param[inout] rng               Random generator state.
  !! @param[out]   output            Allocated noisy samples.
  !! @param[out]   noise_var         Complex-sample noise variance.
  !! @param[in]    true_signal_power Optional active-signal power used when the
  !!                                 input contains zero-padding.
  !!
  !! @note When guard padding is present, estimating power over the padded array
  !! would reduce the measured signal power and underestimate the injected noise.
  !! Passing true_signal_power preserves the intended SNR reference.
  subroutine apply_awgn_channel(input, snr_db, rng, output, noise_var, true_signal_power)
    complex(dp), intent(in)  :: input(:)
    real(dp),    intent(in)  :: snr_db
    type(rng_state), intent(inout) :: rng
    complex(dp), allocatable, intent(out) :: output(:)
    real(dp), intent(out) :: noise_var
    real(dp), intent(in), optional :: true_signal_power
    integer :: i
    real(dp) :: p_sig, sigma

    if (present(true_signal_power)) then
      p_sig = max(true_signal_power, eps)
    else
      p_sig = max(signal_power(input), eps)
    end if

    noise_var = p_sig / max(10.0_dp**(snr_db / 10.0_dp), eps)
    sigma = sqrt(noise_var * 0.5_dp)
    allocate(output(size(input)))
    do i = 1, size(input)
      output(i) = input(i) + sigma * cmplx(rng_normal(rng), rng_normal(rng), dp)
    end do
  end subroutine apply_awgn_channel

  !> Apply carrier-frequency offset and linear frequency drift.
  !!
  !! @param[in]  cfg        Waveform configuration containing the sample rate.
  !! @param[in]  input      Input complex samples.
  !! @param[in]  cfo_hz     Static carrier-frequency offset in hertz.
  !! @param[in]  drift_hz_s Linear frequency drift in hertz per second.
  !! @param[out] output     Allocated impaired samples.
  subroutine apply_cfo_and_drift(cfg, input, cfo_hz, drift_hz_s, output)
    type(waveform_config), intent(in) :: cfg
    complex(dp), intent(in) :: input(:)
    real(dp), intent(in) :: cfo_hz, drift_hz_s
    complex(dp), allocatable, intent(out) :: output(:)
    integer :: i, n
    real(dp) :: t, phase
    n = size(input)
    allocate(output(n))
    !$omp parallel do default(shared) private(i, t, phase) schedule(static)
    do i = 1, n
      t = real(i - 1, dp) / cfg%sample_rate
      phase = two_pi * (cfo_hz * t + 0.5_dp * drift_hz_s * t * t)
      output(i) = input(i) * cmplx(cos(phase), sin(phase), dp)
    end do
    !$omp end parallel do
  end subroutine apply_cfo_and_drift

  !> Add one or more continuous-wave tone interferers.
  !!
  !! @param[in]  cfg    Waveform configuration containing the sample rate.
  !! @param[in]  input  Input complex samples.
  !! @param[in]  tones  Interferer specifications. An empty array is allowed.
  !! @param[out] output Allocated samples with added tones.
  !!
  !! Each complex tone amplitude is set from its requested
  !! signal-to-interference ratio relative to the measured power of the input
  !! sequence, so the tone power is amp**2.
  subroutine apply_tone_interferers(cfg, input, tones, output)
    type(waveform_config), intent(in) :: cfg
    complex(dp), intent(in) :: input(:)
    type(tone_interferer), intent(in) :: tones(:)
    complex(dp), allocatable, intent(out) :: output(:)
    integer :: i, k, n
    real(dp) :: t, phase, p_sig, p_tone, amp
    n = size(input)
    allocate(output(n))
    output = input
    if (size(tones) == 0) return
    p_sig = max(signal_power(input), eps)
    do k = 1, size(tones)
      p_tone = p_sig / max(10.0_dp**(tones(k)%sir_db / 10.0_dp), eps)
      amp = sqrt(p_tone)
      !$omp parallel do default(shared) firstprivate(k, amp) private(i, t, phase) schedule(static)
      do i = 1, n
        t = real(i - 1, dp) / cfg%sample_rate
        if (t < tones(k)%start_s) cycle
        if (tones(k)%has_duration .and. t > tones(k)%start_s + tones(k)%duration_s) cycle
        phase = two_pi * tones(k)%frequency_hz * t + tones(k)%phase_rad
        output(i) = output(i) + amp * cmplx(cos(phase), sin(phase), dp)
      end do
      !$omp end parallel do
    end do
  end subroutine apply_tone_interferers

  !> Apply sample-clock error by linear-interpolation resampling.
  !!
  !! @param[in]  input  Input complex samples.
  !! @param[in]  ppm    Clock offset in parts per million.
  !! @param[out] output Allocated resampled sequence.
  subroutine apply_sample_clock_error(input, ppm, output)
    complex(dp), intent(in) :: input(:)
    real(dp), intent(in) :: ppm
    complex(dp), allocatable, intent(out) :: output(:)
    integer :: n_in, n_out, i, lo, hi
    real(dp) :: scale, src, frac
    n_in = size(input)
    if (abs(ppm) < eps) then
      allocate(output(n_in))
      output = input
      return
    end if
    scale = 1.0_dp + ppm * 1.0e-6_dp
    n_out = int(real(n_in, dp) / scale)
    allocate(output(n_out))
    do i = 1, n_out
      src = real(i - 1, dp) * scale + 1.0_dp
      lo = int(src)
      hi = min(n_in, lo + 1)
      lo = max(1, lo)
      frac = src - real(lo, dp)
      output(i) = (1.0_dp - frac) * input(lo) + frac * input(hi)
    end do
  end subroutine apply_sample_clock_error

  !> Apply the complete simulation channel pipeline.
  !!
  !! @param[in]    cfg       Waveform configuration.
  !! @param[in]    ch        Channel impairment parameters.
  !! @param[in]    tones     Continuous-wave interferer specifications.
  !! @param[in]    input     Transmitted complex samples.
  !! @param[inout] rng       Random generator state for AWGN generation.
  !! @param[out]   output    Allocated received samples.
  !! @param[out]   noise_var AWGN variance used in the final stage.
  !!
  !! The processing order is guard padding, carrier offset/drift, tone
  !! interference, sample-clock error, and AWGN. The active-signal power is
  !! measured before padding so SNR remains referenced to the transmitted frame.
  subroutine simulate_full_channel(cfg, ch, tones, input, rng, output, noise_var)
    type(waveform_config), intent(in) :: cfg
    type(channel_config),  intent(in) :: ch
    type(tone_interferer), intent(in) :: tones(:)
    complex(dp), intent(in) :: input(:)
    type(rng_state), intent(inout) :: rng
    complex(dp), allocatable, intent(out) :: output(:)
    real(dp), intent(out) :: noise_var

    complex(dp), allocatable :: padded(:), shifted(:), with_tones(:), clocked(:)
    integer :: pre, post, i
    real(dp) :: true_p_sig

    ! Measure active-frame power before zero-padding can dilute the SNR reference.
    true_p_sig = max(signal_power(input), eps)

    pre  = nint(ch%prefix_pad_s * cfg%sample_rate)
    post = nint(ch%suffix_pad_s * cfg%sample_rate)
    allocate(padded(pre + size(input) + post))
    padded = (0.0_dp, 0.0_dp)
    do i = 1, size(input)
      padded(pre + i) = input(i)
    end do

    call apply_cfo_and_drift(cfg, padded, ch%cfo_hz, ch%drift_hz_s, shifted)
    deallocate(padded)
    call apply_tone_interferers(cfg, shifted, tones, with_tones)
    deallocate(shifted)
    call apply_sample_clock_error(with_tones, ch%clock_ppm, clocked)
    deallocate(with_tones)
    call apply_awgn_channel(clocked, ch%snr_db, rng, output, noise_var, true_p_sig)
    deallocate(clocked)
  end subroutine simulate_full_channel
end module cdss_channel
