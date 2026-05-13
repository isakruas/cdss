!> @file test_channel_regression.f90
!! @brief Deterministic channel-regression tests for end-to-end frame recovery.
program test_channel_regression
  use cdss_kinds, only: dp
  use cdss_types, only: waveform_config, soft_metric_config, coding_config, &
                          channel_config, tone_interferer, frame_recovery_result
  use cdss_receiver, only: encode_frame, recover_frame
  use cdss_waveform, only: synthesize_waveform
  use cdss_channel, only: simulate_full_channel
  use cdss_rng, only: rng_state, rng_seed
  implicit none

  type(waveform_config) :: cfg
  type(soft_metric_config) :: metric
  type(coding_config) :: coding
  integer, allocatable :: payload(:), coded_bits(:)
  complex(dp), allocatable :: tx(:)

  print '(a)', "test_channel_regression: starting deterministic channel regressions"
  cfg%sample_rate = 2000.0_dp
  cfg%carrier_hz = 250.0_dp
  cfg%chip_rate = 250.0_dp
  cfg%spreading_factor = 8
  cfg%chirp_slope_hz_s = 0.0_dp
  cfg%pulse_taper = 0.0_dp
  metric%stft_window = 128
  metric%stft_hop = 32
  metric%refinement_iters = 2

  allocate(payload(64))
  call fill_payload(payload)
  call encode_frame(cfg, coding, payload, coded_bits)
  call synthesize_waveform(cfg, coded_bits, tx)
  print '(a,i0,a,i0)', "  payload bits=", size(payload), ", TX samples=", size(tx)

  call run_case("snr_35_seed_1", 35.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, .false., 1)
  call run_case("snr_25_seed_2", 25.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, .false., 2)
  call run_case("small_cfo_drift", 35.0_dp, 0.005_dp, 0.0002_dp, 0.0_dp, .false., 3)
  call run_case("tone_interference", 35.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, .true., 4)

  print '(a)', "test_channel_regression: passed"
  deallocate(payload, coded_bits, tx)

contains

  subroutine run_case(label, snr_db, cfo_hz, drift_hz_s, clock_ppm, add_tone, seed)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: snr_db, cfo_hz, drift_hz_s, clock_ppm
    logical, intent(in) :: add_tone
    integer, intent(in) :: seed
    type(channel_config) :: ch
    type(tone_interferer), allocatable :: tones(:)
    type(rng_state) :: rng
    type(frame_recovery_result) :: frame
    complex(dp), allocatable :: rx(:)
    real(dp) :: noise_var
    integer :: errors, i

    ch%snr_db = snr_db
    ch%cfo_hz = cfo_hz
    ch%drift_hz_s = drift_hz_s
    ch%clock_ppm = clock_ppm
    if (add_tone) then
      allocate(tones(1))
      tones(1)%frequency_hz = cfg%carrier_hz + 25.0_dp
      tones(1)%sir_db = 18.0_dp
      tones(1)%phase_rad = 0.25_dp
    else
      allocate(tones(0))
    end if

    call rng_seed(rng, seed)
    call simulate_full_channel(cfg, ch, tones, tx, rng, rx, noise_var)
    call recover_frame(cfg, metric, coding, rx, size(payload) + coding%payload_k, frame)
    errors = size(payload)
    if (frame%payload_bit_count == size(payload)) then
      errors = 0
      do i = 1, size(payload)
        if (payload(i) /= frame%payload_bits(i)) errors = errors + 1
      end do
    end if

    print '(a,a,a,l1,a,i0,a,es12.4)', "  case ", trim(label), ": crc_ok=", &
         frame%crc_ok, ", errors=", errors, ", noise_var=", noise_var
    if (.not. frame%crc_ok .or. errors /= 0) then
      print '(a,a)', "test_channel_regression: failed case ", trim(label)
      print '(a,i0,a,i0)', "  decoded bits=", frame%payload_bit_count, &
           ", expected=", size(payload)
      print '(a,i0,a,i0)', "  CRC expected=", frame%payload_crc_expected, &
           ", observed=", frame%payload_crc_observed
      error stop 1
    end if
    deallocate(tones, rx)
  end subroutine run_case

  subroutine fill_payload(bits)
    integer, intent(out) :: bits(:)
    integer :: i
    do i = 1, size(bits)
      bits(i) = merge(1, 0, mod(i * 7 + i / 3, 11) < 5)
    end do
  end subroutine fill_payload
end program test_channel_regression
