!> @file cdss_receiver.f90
!! @brief Frame encoding, header handling, synchronization, and receiver orchestration.
!!
!! @details
!! CDSS frames are arranged as
!!
!!     preamble | Polar-coded header | Polar-coded payload.
!!
!! The header carries the payload bit length and a CRC-16 checksum. Receiver
!! processing demodulates and decodes the header first, validates the payload
!! length against caller-provided bounds, then demodulates and decodes the
!! payload. Optional feedback iterations can improve interference exposure
!! estimates when the first payload decode is not CRC-clean.
module cdss_receiver
  use cdss_kinds, only: dp
  use cdss_constants, only: eps, frame_header_bits, max_payload_bits, two_pi
  use cdss_types, only: waveform_config, soft_metric_config, coding_config, &
                          soft_demod_result, frame_recovery_result
  use cdss_softdemod, only: soft_demod_bits
  use cdss_iterative, only: refine_soft_demod
  use cdss_polar, only: polar_code, polar_encode_stream, polar_decode_stream
  use cdss_bits, only: payload_crc16, int_to_bits, bits_to_int
  use cdss_chirp_pulse, only: build_chirp_carrier, build_chip_envelope
  use cdss_pn, only: pn_msequence
  use cdss_fft, only: fft_forward, fft_inverse
  implicit none
  private
  public :: encode_frame, recover_frame, sync_preamble
  public :: build_header, decode_header, frame_n_bits, frame_n_samples

contains

  !> Build the uncoded 48-bit frame header.
  !!
  !! @param[in]  payload_bits Payload bits used for length and CRC fields.
  !! @param[out] header       Header bits: 32-bit length followed by 16-bit CRC.
  subroutine build_header(payload_bits, header)
    integer, intent(in)  :: payload_bits(:)
    integer, intent(out) :: header(frame_header_bits)
    integer :: crc, length_bits(32), crc_bits(16)
    call int_to_bits(size(payload_bits), 32, length_bits)
    crc = payload_crc16(payload_bits)
    call int_to_bits(crc, 16, crc_bits)
    header(1:32)  = length_bits
    header(33:48) = crc_bits
  end subroutine build_header

  !> Decode the uncoded 48-bit frame header.
  !!
  !! @param[in]  header_bits    Header bit vector.
  !! @param[out] payload_length Payload length in bits.
  !! @param[out] crc_expected   CRC value carried by the header.
  subroutine decode_header(header_bits, payload_length, crc_expected)
    integer, intent(in)  :: header_bits(:)
    integer, intent(out) :: payload_length, crc_expected
    payload_length = bits_to_int(header_bits(1:32))
    crc_expected   = bits_to_int(header_bits(33:48))
  end subroutine decode_header

  !> Compute the number of coded bits after header and payload coding.
  !!
  !! @param[in] coding            Coding configuration.
  !! @param[in] payload_bit_count Uncoded payload length in bits.
  !! @return Number of coded bits excluding the preamble.
  pure integer function frame_n_bits(coding, payload_bit_count) result(total)
    type(coding_config), intent(in) :: coding
    integer, intent(in) :: payload_bit_count
    integer :: header_n, payload_blocks, payload_n
    header_n = coding%header_n
    payload_blocks = (payload_bit_count + coding%payload_k - 1) / coding%payload_k
    payload_n = payload_blocks * coding%payload_n
    total = header_n + payload_n
  end function frame_n_bits

  !> Compute the number of IQ samples in a complete frame.
  !!
  !! @param[in] cfg               Waveform configuration.
  !! @param[in] coding            Coding configuration.
  !! @param[in] payload_bit_count Uncoded payload length in bits.
  !! @return Sample count including preamble, header, and payload.
  pure integer function frame_n_samples(cfg, coding, payload_bit_count) result(samples)
    type(waveform_config), intent(in) :: cfg
    type(coding_config),   intent(in) :: coding
    integer, intent(in) :: payload_bit_count
    integer :: total_bits
    total_bits = cfg%preamble_bits + frame_n_bits(coding, payload_bit_count)
    samples = total_bits * cfg%samples_per_bit()
  end function frame_n_samples

  !> Encode payload bits into the transmitted coded-bit frame.
  !!
  !! @param[in]  cfg            Waveform configuration, including preamble length.
  !! @param[in]  coding         Polar coding configuration.
  !! @param[in]  payload_bits   Uncoded payload bits.
  !! @param[out] all_coded_bits Allocated bit stream: preamble, header, payload.
  subroutine encode_frame(cfg, coding, payload_bits, all_coded_bits)
    type(waveform_config), intent(in) :: cfg
    type(coding_config),   intent(in) :: coding
    integer, intent(in)  :: payload_bits(:)
    integer, allocatable, intent(out) :: all_coded_bits(:)

    integer :: header(frame_header_bits)
    integer, allocatable :: header_info(:), header_coded(:), payload_coded(:), pre_bits(:)
    integer :: total
    type(polar_code) :: hcode, pcode

    hcode%n = coding%header_n
    hcode%k = coding%header_k
    hcode%list_size = coding%list_size
    pcode%n = coding%payload_n
    pcode%k = coding%payload_k
    pcode%list_size = coding%list_size

    call build_header(payload_bits, header)
    allocate(header_info(coding%header_k))
    header_info = 0
    header_info(1:frame_header_bits) = header

    call polar_encode_stream(hcode, header_info, header_coded)
    call polar_encode_stream(pcode, payload_bits, payload_coded)

    total = cfg%preamble_bits + size(header_coded) + size(payload_coded)
    allocate(all_coded_bits(total))

    ! Transmit a deterministic m-sequence for acquisition.
    call pn_msequence(5, [5, 3], 1, cfg%preamble_bits, pre_bits)
    all_coded_bits(1:cfg%preamble_bits) = pre_bits

    all_coded_bits(cfg%preamble_bits + 1:cfg%preamble_bits + size(header_coded)) = header_coded
    all_coded_bits(cfg%preamble_bits + size(header_coded) + 1:) = payload_coded

    deallocate(header_info, header_coded, payload_coded, pre_bits)
  end subroutine encode_frame

  !> Recover a payload from an aligned received frame.
  !!
  !! @param[in]  cfg              Waveform configuration.
  !! @param[in]  metric           Soft-demodulator configuration.
  !! @param[in]  coding           Polar coding configuration.
  !! @param[in]  rx               Received samples. The first sample is assumed
  !!                              to be the start of the transmitted frame.
  !! @param[in]  max_payload_hint Caller-provided payload length guard in bits.
  !! @param[out] result           Decoded payload, CRC status, and diagnostics.
  !! @param[in]  acq_peak         Optional complex correlation peak from preamble sync.
  subroutine recover_frame(cfg, metric, coding, rx, max_payload_hint, result, acq_peak)
    use cdss_types, only: tracking_state
    type(waveform_config),    intent(in)  :: cfg
    type(soft_metric_config), intent(in)  :: metric
    type(coding_config),      intent(in)  :: coding
    complex(dp),              intent(in)  :: rx(:)
    integer,                  intent(in)  :: max_payload_hint
    type(frame_recovery_result), intent(out) :: result
    complex(dp), intent(in), optional     :: acq_peak

    type(polar_code) :: hcode, pcode
    type(soft_demod_result) :: header_demod, payload_demod, refined
    real(dp), allocatable :: header_llrs_polar(:), payload_llrs_polar(:)
    integer, allocatable :: header_info(:), payload_bits(:)
    integer :: header_n_bits, payload_n_bits, payload_blocks
    integer :: payload_length, crc_expected, observed_crc
    integer :: total_samples_needed, iter, max_iter
    integer :: header_offset, i, preamble_idx
    type(tracking_state) :: state
    complex(dp), allocatable :: rx_used(:)
    complex(dp) :: preamble_peak, phase_rot
    real(dp) :: preamble_score, preamble_cfo_hz, t, phase

    hcode%n = coding%header_n
    hcode%k = coding%header_k
    hcode%list_size = coding%list_size
    pcode%n = coding%payload_n
    pcode%k = coding%payload_k
    pcode%list_size = coding%list_size

    allocate(rx_used(size(rx)))
    rx_used = rx

    ! If the caller did not already provide synchronized/CFO-corrected samples,
    ! use the aligned preamble to remove commercial-radio-scale residual CFO
    ! before the header decoder sees it.
    if (.not. present(acq_peak)) then
      if (size(rx) >= cfg%preamble_bits * cfg%samples_per_bit()) then
        call sync_preamble(cfg, rx, preamble_idx, preamble_score, preamble_peak, preamble_cfo_hz)
        if (abs(preamble_idx) <= 4) then
          if (abs(preamble_peak) > eps) then
            phase_rot = conjg(preamble_peak / abs(preamble_peak))
          else
            phase_rot = (1.0_dp, 0.0_dp)
          end if
          do i = 1, size(rx)
            t = real(i - 1, dp) / cfg%sample_rate
            phase = -two_pi * preamble_cfo_hz * t
            rx_used(i) = rx(i) * cmplx(cos(phase), sin(phase), dp) * phase_rot
          end do
        end if
      end if
    end if

    ! Initialize tracking state. If acq_peak is provided, use its phase.
    state%phase_est = 0.0_dp
    state%freq_est  = 0.0_dp
    state%timing_error = 0.0_dp
    if (present(acq_peak)) then
      if (abs(acq_peak) > eps) then
        state%phase_est = atan2(aimag(acq_peak), real(acq_peak, dp))
      end if
    end if

    header_offset = cfg%preamble_bits

    header_n_bits = coding%header_n
    if (size(rx_used) < (header_offset + header_n_bits) * cfg%samples_per_bit()) then
      call empty_result(result)
      return
    end if

    ! Demodulate header using persistent state.
    call soft_demod_bits(cfg, metric, &
         rx_used(header_offset * cfg%samples_per_bit() + 1:), header_n_bits, header_demod, state=state)
    
    allocate(header_llrs_polar(header_n_bits))
    header_llrs_polar = header_demod%bit_llrs
    call polar_decode_stream(hcode, header_llrs_polar, header_info)

    call decode_header(header_info(1:frame_header_bits), payload_length, crc_expected)
    if (payload_length < 0 .or. payload_length > min(max_payload_bits, max_payload_hint)) then
      call empty_result(result, header_demod)
      return
    end if

    payload_blocks = (payload_length + coding%payload_k - 1) / coding%payload_k
    payload_n_bits = payload_blocks * coding%payload_n
    total_samples_needed = (header_offset + header_n_bits + payload_n_bits) * cfg%samples_per_bit()
    if (size(rx_used) < total_samples_needed) then
      call empty_result(result, header_demod)
      return
    end if

    ! Demodulate payload using state carried over from header.
    call soft_demod_bits(cfg, metric, &
         rx_used((header_offset + header_n_bits) * cfg%samples_per_bit() + 1:), &
         payload_n_bits, payload_demod, state=state)
    
    allocate(payload_llrs_polar(payload_n_bits))
    payload_llrs_polar = payload_demod%bit_llrs
    call polar_decode_stream(pcode, payload_llrs_polar, payload_bits)

    if (size(payload_bits) > payload_length) then
      payload_bits = payload_bits(1:payload_length)
    end if

    observed_crc = payload_crc16(payload_bits)

    max_iter = max(0, metric%refinement_iters)
    iter = 0
    do while (iter < max_iter .and. observed_crc /= crc_expected)
      block
        type(tracking_state) :: state_refine
        state_refine = state 
        call refine_loop_once(cfg, metric, coding, rx_used, header_offset, header_n_bits, &
             payload_bits, payload_n_bits, refined, state_refine)
      end block
      payload_demod = refined
      payload_llrs_polar = payload_demod%bit_llrs
      call polar_decode_stream(pcode, payload_llrs_polar, payload_bits)
      if (size(payload_bits) > payload_length) payload_bits = payload_bits(1:payload_length)
      observed_crc = payload_crc16(payload_bits)
      iter = iter + 1
    end do

    allocate(result%payload_bits(payload_length))
    result%payload_bits = payload_bits(1:payload_length)
    allocate(result%polar_payload_bits(payload_length))
    result%polar_payload_bits = result%payload_bits
    allocate(result%payload_llrs(payload_n_bits))
    result%payload_llrs = payload_llrs_polar
    result%soft_detection = payload_demod
    result%payload_bit_count       = payload_length
    result%payload_crc_expected    = crc_expected
    result%payload_crc_observed    = observed_crc
    result%refinement_iters_used   = iter
    result%crc_ok                  = (observed_crc == crc_expected)
    result%acquisition_sample      = 0
    result%acquisition_score       = 1.0_dp

    deallocate(header_llrs_polar, payload_llrs_polar, header_info)
  end subroutine recover_frame

  !> Locate the frame preamble by matched complex correlation.
  !!
  !! @param[in]  cfg        Waveform configuration.
  !! @param[in]  rx         Received samples to search.
  !! @param[out] best_idx   Zero-based start index of the largest correlation.
  !! @param[out] best_score Magnitude of the largest correlation.
  !! @param[out] best_peak  Optional complex correlation peak, including phase.
  !! @param[out] best_cfo_hz Optional CFO estimate from the preamble in hertz.
  !!
  !! The reference waveform is synthesized from the same m-sequence,
  !! spreading code, chirp carrier, and chip envelope used by the transmitter.
  !! The returned peak phase can be used for bulk phase rotation before demodulation.
  subroutine sync_preamble(cfg, rx, best_idx, best_score, best_peak, best_cfo_hz)
    type(waveform_config), intent(in) :: cfg
    complex(dp), intent(in) :: rx(:)
    integer, intent(out) :: best_idx
    real(dp), intent(out) :: best_score
    complex(dp), intent(out), optional :: best_peak
    real(dp), intent(out), optional :: best_cfo_hz

    complex(dp), allocatable :: ref(:), bit_ref(:), carrier(:), fft_rx(:), fft_ref(:)
    complex(dp), allocatable :: corr_freq(:), corr_time(:), rx_pad(:), ref_pad(:)
    real(dp), allocatable :: envelope(:), pn_signs(:)
    integer, allocatable :: pre_bits(:)
    integer :: i, j, k, n_samples, n_search, samples_per_bit, samples_per_chip, nfft
    real(dp) :: bit_sign, chip_sign
    complex(dp) :: corr, peak
    real(dp) :: cfo_hz

    samples_per_chip = cfg%samples_per_chip()
    samples_per_bit  = cfg%samples_per_bit()
    n_samples = cfg%preamble_bits * samples_per_bit

    allocate(ref(n_samples), bit_ref(samples_per_bit), carrier(samples_per_bit), envelope(samples_per_chip), &
             pn_signs(cfg%spreading_factor))

    call build_chirp_carrier(cfg, carrier)
    call build_chip_envelope(samples_per_chip, cfg%pulse_taper, envelope)

    block
      use cdss_waveform, only: get_pn_signs
      call get_pn_signs(cfg, pn_signs)
    end block

    ref = (0.0_dp, 0.0_dp)
    bit_ref = (0.0_dp, 0.0_dp)
    do j = 1, cfg%spreading_factor
      chip_sign = pn_signs(j)
      do k = 1, samples_per_chip
        bit_ref((j-1)*samples_per_chip + k) = &
             chip_sign * envelope(k) * carrier((j-1)*samples_per_chip + k)
      end do
    end do

    call pn_msequence(5, [5, 3], 1, cfg%preamble_bits, pre_bits)
    do i = 1, cfg%preamble_bits
      ! Correct BPSK mapping: bit 1 maps to +1, bit 0 maps to -1.
      bit_sign = merge(1.0_dp, -1.0_dp, pre_bits(i) == 1)
      do j = 1, cfg%spreading_factor
        chip_sign = bit_sign * pn_signs(j)
        do k = 1, samples_per_chip
          ref((i-1)*samples_per_bit + (j-1)*samples_per_chip + k) = &
               chip_sign * envelope(k) * carrier((j-1)*samples_per_chip + k)
        end do
      end do
    end do

    best_idx = 0
    best_score = -1.0_dp
    peak = (0.0_dp, 0.0_dp)
    cfo_hz = 0.0_dp
    n_search = size(rx) - n_samples
    if (n_search < 0) then
      if (present(best_peak)) best_peak = peak
      if (present(best_cfo_hz)) best_cfo_hz = cfo_hz
      deallocate(ref, bit_ref, carrier, envelope, pn_signs, pre_bits)
      return
    end if

    ! First acquire timing with a CFO-tolerant differential preamble metric.  A
    ! full coherent preamble matched filter assumes zero CFO; over the default
    ! 7.9 s preamble, even 0.5 Hz rotates through several cycles and can cancel
    ! the true peak.  The one-bit matched-filter stream preserves the known
    ! preamble bit signs while allowing a constant inter-bit phase increment.
    nfft = 1
    do while (nfft < size(rx) + samples_per_bit - 1)
      nfft = nfft * 2
    end do
    allocate(rx_pad(nfft), ref_pad(nfft), corr_freq(nfft))
    rx_pad = (0.0_dp, 0.0_dp)
    ref_pad = (0.0_dp, 0.0_dp)
    rx_pad(1:size(rx)) = rx
    do j = 1, samples_per_bit
      ref_pad(j) = conjg(bit_ref(samples_per_bit - j + 1))
    end do
    call fft_forward(rx_pad, fft_rx)
    call fft_forward(ref_pad, fft_ref)
    corr_freq = fft_rx * fft_ref
    call fft_inverse(corr_freq, corr_time)

    do i = 0, n_search
      block
        complex(dp) :: prev, curr, phase_sum
        real(dp) :: local_score
        phase_sum = (0.0_dp, 0.0_dp)
        bit_sign = merge(1.0_dp, -1.0_dp, pre_bits(1) == 1)
        prev = bit_sign * corr_time(i + samples_per_bit)
        do j = 2, cfg%preamble_bits
          bit_sign = merge(1.0_dp, -1.0_dp, pre_bits(j) == 1)
          curr = bit_sign * corr_time(i + (j - 1) * samples_per_bit + samples_per_bit)
          phase_sum = phase_sum + curr * conjg(prev)
          prev = curr
        end do
        local_score = abs(phase_sum)
        if (local_score > best_score) then
          best_score = local_score
          best_idx = i
        end if
      end block
    end do

    deallocate(fft_rx, fft_ref, corr_freq, corr_time, rx_pad, ref_pad)

    call refine_preamble_cfo(cfg, rx, ref, best_idx, cfo_hz, peak, best_score)
    call refine_preamble_timing_with_cfo(cfg, rx, ref, cfo_hz, best_idx, peak, best_score)
    call refine_preamble_cfo(cfg, rx, ref, best_idx, cfo_hz, peak, best_score)

    if (present(best_peak)) best_peak = peak
    if (present(best_cfo_hz)) best_cfo_hz = cfo_hz

    deallocate(ref, bit_ref, carrier, envelope, pn_signs, pre_bits)
  end subroutine sync_preamble

  !> Refine CFO at a known frame candidate by direct preamble correlation.
  subroutine refine_preamble_cfo(cfg, rx, ref, start_idx, cfo_hz, peak, score)
    type(waveform_config), intent(in) :: cfg
    complex(dp), intent(in) :: rx(:), ref(:)
    integer, intent(in) :: start_idx
    real(dp), intent(out) :: cfo_hz, score
    complex(dp), intent(out) :: peak
    integer :: fbin, s
    real(dp), parameter :: cfo_limit_hz = 8.0_dp
    real(dp), parameter :: cfo_step_hz = 0.05_dp
    real(dp) :: f_hz, t, local_score
    complex(dp) :: corr, rot

    cfo_hz = 0.0_dp
    score = -1.0_dp
    peak = (0.0_dp, 0.0_dp)
    if (start_idx < 0 .or. start_idx + size(ref) > size(rx)) return

    do fbin = nint(-cfo_limit_hz / cfo_step_hz), nint(cfo_limit_hz / cfo_step_hz)
      f_hz = real(fbin, dp) * cfo_step_hz
      corr = (0.0_dp, 0.0_dp)
      do s = 1, size(ref)
        t = real(start_idx + s - 1, dp) / cfg%sample_rate
        rot = cmplx(cos(-two_pi * f_hz * t), sin(-two_pi * f_hz * t), dp)
        corr = corr + rx(start_idx + s) * rot * conjg(ref(s))
      end do
      local_score = abs(corr)
      if (local_score > score) then
        score = local_score
        cfo_hz = f_hz
        peak = corr
      end if
    end do
  end subroutine refine_preamble_cfo

  !> Refine frame timing after derotating by a coarse CFO estimate.
  subroutine refine_preamble_timing_with_cfo(cfg, rx, ref, cfo_hz, best_idx, peak, score)
    type(waveform_config), intent(in) :: cfg
    complex(dp), intent(in) :: rx(:), ref(:)
    real(dp), intent(in) :: cfo_hz
    integer, intent(inout) :: best_idx
    complex(dp), intent(out) :: peak
    real(dp), intent(out) :: score

    complex(dp), allocatable :: fft_rx(:), fft_ref(:), corr_freq(:), corr_time(:)
    complex(dp), allocatable :: rx_pad(:), ref_pad(:)
    integer :: i, j, n_search, nfft
    real(dp) :: t, phase, local_score
    complex(dp) :: corr

    n_search = size(rx) - size(ref)
    if (n_search < 0) return

    nfft = 1
    do while (nfft < size(rx) + size(ref) - 1)
      nfft = nfft * 2
    end do
    allocate(rx_pad(nfft), ref_pad(nfft), corr_freq(nfft))
    rx_pad = (0.0_dp, 0.0_dp)
    ref_pad = (0.0_dp, 0.0_dp)

    do i = 1, size(rx)
      t = real(i - 1, dp) / cfg%sample_rate
      phase = -two_pi * cfo_hz * t
      rx_pad(i) = rx(i) * cmplx(cos(phase), sin(phase), dp)
    end do
    do j = 1, size(ref)
      ref_pad(j) = conjg(ref(size(ref) - j + 1))
    end do

    call fft_forward(rx_pad, fft_rx)
    call fft_forward(ref_pad, fft_ref)
    corr_freq = fft_rx * fft_ref
    call fft_inverse(corr_freq, corr_time)

    score = -1.0_dp
    peak = (0.0_dp, 0.0_dp)
    do i = 0, n_search
      corr = corr_time(i + size(ref))
      local_score = abs(corr)
      if (local_score > score) then
        score = local_score
        best_idx = i
        peak = corr
      end if
    end do

    deallocate(fft_rx, fft_ref, corr_freq, corr_time, rx_pad, ref_pad)
  end subroutine refine_preamble_timing_with_cfo

  !> Perform one payload feedback-refinement iteration.
  !!
  !! @param[in]  cfg             Waveform configuration.
  !! @param[in]  metric          Soft-demodulator configuration.
  !! @param[in]  coding          Polar coding configuration.
  !! @param[in]  rx              Full received frame samples.
  !! @param[in]  header_offset   Header start offset in coded bits.
  !! @param[in]  header_n_bits   Coded header length in bits.
  !! @param[in]  payload_bits    Current decoded payload estimate.
  !! @param[in]  payload_n_bits  Coded payload length in bits.
  !! @param[out] refined         Updated soft-demodulation result.
  !! @param[inout] state         Optional tracking state.
  subroutine refine_loop_once(cfg, metric, coding, rx, header_offset, header_n_bits, &
                                payload_bits, payload_n_bits, refined, state)
    use cdss_types, only: tracking_state
    type(waveform_config),    intent(in) :: cfg
    type(soft_metric_config), intent(in) :: metric
    type(coding_config),      intent(in) :: coding
    complex(dp),              intent(in) :: rx(:)
    integer, intent(in) :: header_offset, header_n_bits, payload_n_bits
    integer, intent(in) :: payload_bits(:)
    type(soft_demod_result), intent(out) :: refined
    type(tracking_state), intent(inout), optional :: state
    integer, allocatable :: payload_recoded(:)
    type(polar_code) :: pcode
    pcode%n = coding%payload_n
    pcode%k = coding%payload_k
    pcode%list_size = coding%list_size
    call polar_encode_stream(pcode, payload_bits, payload_recoded)
    if (size(payload_recoded) /= payload_n_bits) then
      error stop "refine_loop_once: re-encoded payload length mismatch"
    end if
    call refine_soft_demod(cfg, metric, &
         rx((header_offset + header_n_bits) * cfg%samples_per_bit() + 1:), &
         payload_n_bits, payload_recoded, refined, state=state)
    deallocate(payload_recoded)
  end subroutine refine_loop_once

  !> Initialize an empty failed-recovery result.
  !!
  !! @param[out] result Result object to initialize.
  !! @param[in]  soft   Optional soft-demodulator diagnostics to preserve.
  subroutine empty_result(result, soft)
    type(frame_recovery_result), intent(out) :: result
    type(soft_demod_result), intent(in), optional :: soft
    if (present(soft)) result%soft_detection = soft
    allocate(result%payload_bits(0))
    allocate(result%polar_payload_bits(0))
    allocate(result%payload_llrs(0))
    result%payload_bit_count = 0
    result%crc_ok = .false.
  end subroutine empty_result
end module cdss_receiver
