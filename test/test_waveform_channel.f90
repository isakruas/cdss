!> @file test_waveform_channel.f90
!! @brief Unit tests for waveform synthesis and channel impairments.
program test_waveform_channel
  use cdss_kinds, only: dp
  use cdss_constants, only: pi
  use cdss_types, only: waveform_config, channel_config, tone_interferer
  use cdss_waveform, only: synthesize_waveform, get_pn_signs, get_carrier, get_envelope
  use cdss_channel, only: signal_power, apply_awgn_channel, apply_cfo_and_drift, &
                            apply_sample_clock_error, apply_tone_interferers, &
                            simulate_full_channel
  use cdss_rng, only: rng_state, rng_seed
  implicit none

  type(waveform_config) :: cfg
  type(channel_config) :: ch
  type(tone_interferer), allocatable :: tones(:), empty_tones(:)
  type(rng_state) :: rng
  integer, allocatable :: bits(:)
  complex(dp), allocatable :: iq1(:), iq2(:), noisy(:), shifted(:), clocked(:)
  complex(dp), allocatable :: with_tone(:), channel_out(:), constant(:)
  real(dp), allocatable :: pn_signs(:), envelope(:)
  complex(dp), allocatable :: carrier(:)
  real(dp) :: p_sig, noise_var, measured_noise, expected_noise, measured_tone_power, max_err
  integer :: n, expected_samples

  print '(a)', "test_waveform_channel: starting waveform/channel checks"

  cfg%sample_rate = 2000.0_dp
  cfg%carrier_hz = 250.0_dp
  cfg%chip_rate = 250.0_dp
  cfg%spreading_factor = 8
  cfg%chirp_slope_hz_s = 0.0_dp
  cfg%pulse_taper = 0.0_dp
  allocate(bits(16))
  bits = [(mod(n, 2), n = 1, size(bits))]

  call synthesize_waveform(cfg, bits, iq1)
  call synthesize_waveform(cfg, bits, iq2)
  expected_samples = size(bits) * cfg%samples_per_bit()
  print '(a,i0,a,i0)', "  synthesized samples=", size(iq1), ", expected=", expected_samples
  if (size(iq1) /= expected_samples) then
    print '(a)', "test_waveform_channel: synthesized sample count mismatch"
    error stop 1
  end if
  max_err = maxval(abs(iq1 - iq2))
  print '(a,es12.4)', "  waveform determinism max error=", max_err
  if (max_err > 1.0e-12_dp) then
    print '(a)', "test_waveform_channel: waveform synthesis is not deterministic"
    error stop 2
  end if

  call get_pn_signs(cfg, pn_signs)
  call get_carrier(cfg, carrier)
  call get_envelope(cfg, envelope)
  print '(a,i0,a,i0,a,i0)', "  PN/carrier/envelope sizes: ", size(pn_signs), ", ", &
       size(carrier), ", ", size(envelope)
  if (size(pn_signs) /= cfg%spreading_factor) error stop 3
  if (size(carrier) /= cfg%samples_per_bit()) error stop 4
  if (size(envelope) /= cfg%samples_per_chip()) error stop 5
  if (abs(signal_power(iq1) - 1.0_dp) > 1.0e-12_dp) then
    print '(a,es12.4)', "test_waveform_channel: unexpected waveform power=", signal_power(iq1)
    error stop 6
  end if

  allocate(constant(20000))
  constant = (1.0_dp, 0.0_dp)
  call rng_seed(rng, 1234)
  call apply_awgn_channel(constant, 10.0_dp, rng, noisy, noise_var)
  expected_noise = 0.1_dp
  measured_noise = signal_power(noisy - constant)
  print '(a,es12.4,a,es12.4)', "  AWGN variance returned=", noise_var, &
       ", measured=", measured_noise
  if (abs(noise_var - expected_noise) > 1.0e-12_dp) then
    print '(a)', "test_waveform_channel: AWGN returned variance mismatch"
    error stop 7
  end if
  if (abs(measured_noise - expected_noise) > 1.5e-2_dp) then
    print '(a)', "test_waveform_channel: measured AWGN power outside tolerance"
    error stop 8
  end if

  cfg%sample_rate = 8.0_dp
  call apply_cfo_and_drift(cfg, constant(1:8), 2.0_dp, 0.0_dp, shifted)
  max_err = abs(shifted(2) - cmplx(0.0_dp, 1.0_dp, dp)) + &
            abs(shifted(3) - cmplx(-1.0_dp, 0.0_dp, dp))
  print '(a,es12.4)', "  CFO quarter-cycle check error=", max_err
  if (max_err > 1.0e-12_dp) then
    print '(a)', "test_waveform_channel: CFO phase rotation check failed"
    error stop 9
  end if

  call apply_sample_clock_error(constant(1:1000), 0.0_dp, clocked)
  if (size(clocked) /= 1000 .or. maxval(abs(clocked - constant(1:1000))) > 0.0_dp) then
    print '(a)', "test_waveform_channel: zero-ppm clock path should be identity"
    error stop 10
  end if
  call apply_sample_clock_error(constant(1:1000), 1000.0_dp, clocked)
  print '(a,i0)', "  +1000 ppm output samples=", size(clocked)
  if (size(clocked) /= int(1000.0_dp / 1.001_dp)) then
    print '(a)', "test_waveform_channel: sample-clock output length mismatch"
    error stop 11
  end if

  cfg%sample_rate = 2000.0_dp
  allocate(tones(1))
  allocate(empty_tones(0))
  tones(1)%frequency_hz = 0.0_dp
  tones(1)%sir_db = 0.0_dp
  tones(1)%phase_rad = 0.0_dp
  call apply_tone_interferers(cfg, constant(1:128), tones, with_tone)
  measured_tone_power = signal_power(with_tone - constant(1:128))
  print '(a,es12.4,a,es12.4)', "  tone-interfered first sample magnitude=", &
       abs(with_tone(1)), ", measured tone power=", measured_tone_power
  if (abs(with_tone(1) - constant(1)) < 0.5_dp) then
    print '(a)', "test_waveform_channel: tone interferer did not modify the signal"
    error stop 12
  end if
  if (abs(measured_tone_power - signal_power(constant(1:128))) > 1.0e-12_dp) then
    print '(a)', "test_waveform_channel: 0 dB complex tone SIR should match signal power"
    error stop 13
  end if

  ch%snr_db = 120.0_dp
  ch%prefix_pad_s = 0.002_dp
  ch%suffix_pad_s = 0.003_dp
  call rng_seed(rng, 22)
  call simulate_full_channel(cfg, ch, empty_tones, iq1, rng, channel_out, noise_var)
  p_sig = max(signal_power(iq1), 1.0e-12_dp)
  expected_samples = size(iq1) + nint(ch%prefix_pad_s * cfg%sample_rate) + &
       nint(ch%suffix_pad_s * cfg%sample_rate)
  print '(a,i0,a,i0,a,es12.4)', "  full channel samples=", size(channel_out), &
       ", expected=", expected_samples, ", noise_var=", noise_var
  if (size(channel_out) /= expected_samples) then
    print '(a)', "test_waveform_channel: full channel padding length mismatch"
    error stop 14
  end if
  if (abs(noise_var - p_sig / 10.0_dp**12) > 1.0e-18_dp) then
    print '(a)', "test_waveform_channel: full channel noise variance mismatch"
    error stop 15
  end if

  print '(a)', "test_waveform_channel: passed"
  deallocate(bits, iq1, iq2, noisy, shifted, clocked, with_tone, channel_out, constant)
  deallocate(pn_signs, carrier, envelope, tones, empty_tones)
end program test_waveform_channel
