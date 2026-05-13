!> @file test_receiver_pipeline.f90
!! @brief Receiver tests for soft demodulation, uncertainty, sync, and guards.
program test_receiver_pipeline
  use cdss_kinds, only: dp
  use cdss_constants, only: two_pi
  use cdss_types, only: waveform_config, soft_metric_config, coding_config, &
                          frame_recovery_result, soft_demod_result
  use cdss_waveform, only: synthesize_waveform
  use cdss_softdemod, only: soft_demod_bits, soft_demod_noise_estimate
  use cdss_uncertainty, only: compute_chip_uncertainty
  use cdss_iterative, only: refine_soft_demod
  use cdss_receiver, only: sync_preamble, recover_frame, frame_n_samples
  implicit none

  type(waveform_config) :: cfg
  type(soft_metric_config) :: metric
  type(coding_config) :: coding
  type(soft_demod_result) :: demod, refined
  type(frame_recovery_result) :: frame
  integer, allocatable :: bits(:), decoded(:), payload(:)
  complex(dp), allocatable :: iq(:), padded(:), aligned(:), tone(:), short_rx(:)
  real(dp), allocatable :: chip_j(:)
  integer :: i, prefix, best_idx, n_frame
  real(dp) :: best_score, phase
  complex(dp) :: best_peak, phase_rot

  print '(a)', "test_receiver_pipeline: starting receiver signal-path checks"

  cfg%sample_rate = 2000.0_dp
  cfg%carrier_hz = 250.0_dp
  cfg%chip_rate = 250.0_dp
  cfg%spreading_factor = 8
  cfg%chirp_slope_hz_s = 0.0_dp
  cfg%pulse_taper = 0.0_dp
  metric%has_noise_var = .true.
  metric%noise_var = 1.0e-3_dp
  metric%stft_window = 128
  metric%stft_hop = 32
  metric%llr_clip = 100.0_dp

  call test_soft_demodulation()
  call test_uncertainty_map()
  call test_sync_and_guards()
  call test_variable_preamble()

  print '(a)', "test_receiver_pipeline: passed"

contains

  subroutine test_soft_demodulation()
    integer :: mismatches
    type(soft_metric_config) :: metric_low_noise, metric_high_noise
    type(soft_demod_result) :: demod_low_noise, demod_high_noise
    real(dp) :: llr_low_noise, llr_high_noise
    allocate(bits(16))
    bits = [0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0]
    call synthesize_waveform(cfg, bits, iq)
    print '(a,i0,a,i0)', "  soft-demod bits=", size(bits), ", samples=", size(iq)

    if (abs(soft_demod_noise_estimate(iq(1:1)) - 1.0_dp) > 1.0e-12_dp) then
      print '(a)', "test_receiver_pipeline: one-sample noise estimate should default to 1"
      error stop 1
    end if

    call soft_demod_bits(cfg, metric, iq, size(bits), demod)
    decoded = demod%hard_bits
    mismatches = count(decoded /= bits)
    print '(a,i0,a,es12.4,a,es12.4)', "  soft-demod mismatches=", mismatches, &
         ", min LLR=", minval(demod%bit_llrs), ", max LLR=", maxval(demod%bit_llrs)
    if (mismatches /= 0) then
      print '(a)', "test_receiver_pipeline: noiseless soft demod produced wrong hard bits"
      error stop 2
    end if
    if (size(demod%chip_J) /= size(bits) * cfg%spreading_factor) then
      print '(a)', "test_receiver_pipeline: chip_J size mismatch"
      error stop 3
    end if

    metric_low_noise = metric
    metric_high_noise = metric
    metric_low_noise%llr_clip = 1.0e9_dp
    metric_high_noise%llr_clip = 1.0e9_dp
    metric_low_noise%noise_var = 1.0e-2_dp
    metric_high_noise%noise_var = 1.0_dp
    call soft_demod_bits(cfg, metric_low_noise, iq, size(bits), demod_low_noise)
    call soft_demod_bits(cfg, metric_high_noise, iq, size(bits), demod_high_noise)
    llr_low_noise = sum(abs(demod_low_noise%bit_llrs)) / real(size(bits), dp)
    llr_high_noise = sum(abs(demod_high_noise%bit_llrs)) / real(size(bits), dp)
    print '(a,es12.4,a,es12.4)', "  LLR noise scaling: low-noise mean=", &
         llr_low_noise, ", high-noise mean=", llr_high_noise
    if (llr_low_noise <= 10.0_dp * llr_high_noise) then
      print '(a)', "test_receiver_pipeline: LLR magnitude should scale with supplied noise variance"
      error stop 4
    end if

    call refine_soft_demod(cfg, metric, iq, size(bits), bits, refined)
    decoded = refined%hard_bits
    mismatches = count(decoded /= bits)
    print '(a,i0)', "  refinement mismatches=", mismatches
    if (mismatches /= 0) then
      print '(a)', "test_receiver_pipeline: refinement changed a noiseless correct decode"
      error stop 5
    end if

    deallocate(bits, decoded, iq)
  end subroutine test_soft_demodulation

  subroutine test_uncertainty_map()
    integer :: n_samples
    n_samples = 32 * cfg%samples_per_chip()
    allocate(tone(n_samples))
    tone = (0.0_dp, 0.0_dp)
    call compute_chip_uncertainty(cfg, metric, tone, chip_j)
    print '(a,es12.4)', "  zero-signal uncertainty max=", maxval(chip_j)
    if (maxval(chip_j) > 1.0e-12_dp) then
      print '(a)', "test_receiver_pipeline: zero signal should have zero exposure"
      error stop 6
    end if

    do i = 1, n_samples
      tone(i) = cmplx(cos(two_pi * cfg%carrier_hz * real(i - 1, dp) / cfg%sample_rate), &
                      sin(two_pi * cfg%carrier_hz * real(i - 1, dp) / cfg%sample_rate), dp)
    end do
    call compute_chip_uncertainty(cfg, metric, tone, chip_j)
    print '(a,es12.4,a,es12.4)', "  tone uncertainty max=", maxval(chip_j), &
         ", mean=", sum(chip_j) / real(size(chip_j), dp)
    if (maxval(chip_j) < 1.0_dp) then
      print '(a)', "test_receiver_pipeline: carrier tone should raise chip exposure"
      error stop 7
    end if
    deallocate(tone, chip_j)
  end subroutine test_uncertainty_map

  subroutine test_sync_and_guards()
    allocate(payload(16))
    payload = [(mod(i, 2), i = 1, size(payload))]
    call encode_and_pad_with_phase(payload, 53, 0.4_dp, padded, n_frame)

    call sync_preamble(cfg, padded, best_idx, best_score, best_peak)
    print '(a,i0,a,es12.4)', "  sync index=", best_idx, ", expected=53, score=", best_score
    if (best_idx /= 53) then
      print '(a)', "test_receiver_pipeline: sync did not recover the known prefix length"
      error stop 8
    end if

    phase_rot = conjg(best_peak / abs(best_peak))
    allocate(aligned(n_frame))
    aligned = padded(best_idx + 1:best_idx + n_frame) * phase_rot
    call recover_frame(cfg, metric, coding, aligned, 128, frame)
    print '(a,l1,a,i0)', "  synced recovery crc_ok=", frame%crc_ok, &
         ", payload bits=", frame%payload_bit_count
    if (.not. frame%crc_ok .or. any(frame%payload_bits /= payload)) then
      print '(a)', "test_receiver_pipeline: synced noiseless recovery failed"
      error stop 9
    end if

    allocate(short_rx(cfg%samples_per_bit()))
    short_rx = aligned(1:cfg%samples_per_bit())
    call recover_frame(cfg, metric, coding, short_rx, 128, frame)
    print '(a,i0,a,l1)', "  truncated recovery payload bits=", frame%payload_bit_count, &
         ", crc_ok=", frame%crc_ok
    if (frame%payload_bit_count /= 0 .or. frame%crc_ok) then
      print '(a)', "test_receiver_pipeline: truncated input should produce an empty failed result"
      error stop 10
    end if

    call recover_frame(cfg, metric, coding, aligned, 8, frame)
    print '(a,i0,a,l1)', "  max-hint rejection payload bits=", frame%payload_bit_count, &
         ", crc_ok=", frame%crc_ok
    if (frame%payload_bit_count /= 0 .or. frame%crc_ok) then
      print '(a)', "test_receiver_pipeline: max_payload_hint guard did not reject the frame"
      error stop 11
    end if

    deallocate(payload, padded, aligned, short_rx)
  end subroutine test_sync_and_guards

  subroutine test_variable_preamble()
    integer :: original_preamble
    original_preamble = cfg%preamble_bits
    cfg%preamble_bits = 17
    allocate(payload(12))
    payload = [(merge(1, 0, mod(5 * i + 2, 7) < 3), i = 1, size(payload))]
    call encode_and_pad_with_phase(payload, 37, -0.25_dp, padded, n_frame)

    call sync_preamble(cfg, padded, best_idx, best_score, best_peak)
    print '(a,i0,a,es12.4)', "  variable preamble sync index=", &
         best_idx, ", expected=37, score=", best_score
    if (best_idx /= 37) then
      print '(a)', "test_receiver_pipeline: sync failed for non-default preamble length"
      error stop 12
    end if

    phase_rot = conjg(best_peak / abs(best_peak))
    allocate(aligned(n_frame))
    aligned = padded(best_idx + 1:best_idx + n_frame) * phase_rot
    call recover_frame(cfg, metric, coding, aligned, 64, frame)
    print '(a,l1,a,i0)', "  variable preamble recovery crc_ok=", frame%crc_ok, &
         ", payload bits=", frame%payload_bit_count
    if (.not. frame%crc_ok .or. any(frame%payload_bits /= payload)) then
      print '(a)', "test_receiver_pipeline: non-default preamble recovery failed"
      error stop 13
    end if

    deallocate(payload, padded, aligned)
    cfg%preamble_bits = original_preamble
  end subroutine test_variable_preamble

  subroutine encode_and_pad_with_phase(payload_bits, pad, phase_rad, out, frame_samples)
    integer, intent(in) :: payload_bits(:), pad
    real(dp), intent(in) :: phase_rad
    complex(dp), allocatable, intent(out) :: out(:)
    integer, intent(out) :: frame_samples
    complex(dp), allocatable :: frame_iq(:)
    complex(dp) :: rot
    call build_frame_samples(payload_bits, frame_iq)
    frame_samples = frame_n_samples(cfg, coding, size(payload_bits))
    if (frame_samples /= size(frame_iq)) then
      print '(a)', "test_receiver_pipeline: frame_n_samples disagrees with synthesized frame"
      error stop 14
    end if
    rot = cmplx(cos(phase_rad), sin(phase_rad), dp)
    allocate(out(pad + size(frame_iq)))
    out(1:pad) = (0.0_dp, 0.0_dp)
    out(pad + 1:) = frame_iq * rot
    deallocate(frame_iq)
  end subroutine encode_and_pad_with_phase

  subroutine build_frame_samples(payload_bits, frame_iq)
    use cdss_receiver, only: encode_frame
    integer, intent(in) :: payload_bits(:)
    complex(dp), allocatable, intent(out) :: frame_iq(:)
    integer, allocatable :: coded_bits(:)
    call encode_frame(cfg, coding, payload_bits, coded_bits)
    call synthesize_waveform(cfg, coded_bits, frame_iq)
    deallocate(coded_bits)
  end subroutine build_frame_samples

end program test_receiver_pipeline
