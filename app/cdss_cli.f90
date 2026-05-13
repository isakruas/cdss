!> @file cdss_cli.f90
!! @brief Command-line interface for CDSS waveform generation and recovery.
!!
!! @details
!! The executable exposes a small JSON-emitting interface for reproducible
!! experiments: modem geometry reporting, modulation to waveform files,
!! synchronization, demodulation, and end-to-end channel simulation.
program cdss_cli
  use, intrinsic :: iso_fortran_env, only: error_unit
  use cdss_kinds, only: dp
  use cdss_types, only: waveform_config, soft_metric_config, coding_config, &
                          channel_config, tone_interferer, frame_recovery_result
  use cdss_api,   only: cdss_modem
  use cdss_io,    only: save_waveform, load_waveform, infer_format
  use cdss_channel, only: simulate_full_channel
  use cdss_rng,   only: rng_state, rng_seed, rng_bits
  use cdss_bits,  only: bits_to_string, bits_from_string
  implicit none

  character(len=32) :: command
  integer :: argc

  argc = command_argument_count()
  if (argc < 1) then
    call usage()
    stop 2, quiet=.true.
  end if
  call get_command_argument(1, command)
  select case (trim(command))
  case ("help", "--help", "-h")
    call usage()
  case ("info")
    call cmd_info()
  case ("mod")
    call cmd_mod()
  case ("sync")
    call cmd_sync()
  case ("demod")
    call cmd_demod()
  case ("simulate")
    call cmd_simulate()
  case default
    write(error_unit, '(a,a)') "error: unknown command: ", trim(command)
    call usage()
    stop 2, quiet=.true.
  end select

contains

  !> Synchronize an input waveform to the detected CDSS preamble and save it.
  !!
  !! Options:
  !! --input selects the source waveform, --output selects the synchronized
  !! waveform path, and --bits supplies the expected payload size used to compute
  !! the copied frame length.
  subroutine cmd_sync()
    use cdss_receiver, only: sync_preamble, frame_n_samples
    use cdss_constants, only: eps
    type(cdss_modem) :: modem
    complex(dp), allocatable :: iq(:), synced(:)
    character(len=512) :: input, output, fmt
    integer :: acq_idx, nbits, total_samples
    real(dp) :: acq_score
    complex(dp) :: acq_peak, phase_rot
    input  = "frame.cf32"
    output = "sync.cf32"
    fmt    = "auto"
    nbits  = 64 ! Default payload size for synchronization experiments.
    call parse_string_option_default("--input", input)
    call parse_string_option_default("--output", output)
    call parse_int_option_default("--bits", nbits)
    
    fmt = infer_format(trim(input), "auto")
    call load_waveform(trim(input), trim(fmt), modem%waveform, iq)
    
    call sync_preamble(modem%waveform, iq, acq_idx, acq_score, acq_peak)
    
    if (abs(acq_peak) > eps) then
      phase_rot = conjg(acq_peak / abs(acq_peak))
    else
      phase_rot = (1.0_dp, 0.0_dp)
    end if

    ! Calculate the expected frame length after acquisition.
    total_samples = frame_n_samples(modem%waveform, modem%coding, nbits)
    
    if (size(iq) >= acq_idx + total_samples) then
      allocate(synced(total_samples))
      synced = iq(acq_idx + 1 : acq_idx + total_samples) * phase_rot
      fmt = infer_format(trim(output), "auto")
      call save_waveform(trim(output), trim(fmt), modem%waveform, synced)
    else if (size(iq) > acq_idx) then
      ! If the capture is truncated, save the available suffix for diagnostics.
      allocate(synced(size(iq) - acq_idx))
      synced = iq(acq_idx + 1 :) * phase_rot
      fmt = infer_format(trim(output), "auto")
      call save_waveform(trim(output), trim(fmt), modem%waveform, synced)
    end if
    
    block
      character(len=512) :: jbuf
      integer :: jpos
      call json_begin(jbuf, jpos)
      call json_int(jbuf, jpos, "detected_sample", acq_idx)
      call json_real(jbuf, jpos, "acquisition_score", acq_score)
      call json_int(jbuf, jpos, "frame_samples", total_samples)
      call json_str(jbuf, jpos, "output", trim(output))
      call json_end(jbuf, jpos)
      print '(a)', jbuf(:jpos)
    end block
  end subroutine cmd_sync

  !> Print modem geometry and default receiver settings as JSON.
  subroutine cmd_info()
    type(cdss_modem) :: modem
    character(len=512) :: jbuf
    integer :: jpos
    call json_begin(jbuf, jpos)
    call json_real(jbuf, jpos, "sample_rate_hz",       modem%waveform%sample_rate)
    call json_real(jbuf, jpos, "carrier_hz",           modem%waveform%carrier_hz)
    call json_real(jbuf, jpos, "chip_rate_hz",         modem%waveform%chip_rate)
    call json_int (jbuf, jpos, "spreading_factor",     modem%waveform%spreading_factor)
    call json_real(jbuf, jpos, "chirp_slope_hz_s",     modem%waveform%chirp_slope_hz_s)
    call json_int (jbuf, jpos, "samples_per_chip",     modem%waveform%samples_per_chip())
    call json_int (jbuf, jpos, "samples_per_bit",      modem%waveform%samples_per_bit())
    call json_real(jbuf, jpos, "bandwidth_hz",         modem%waveform%bandwidth_hz())
    call json_real(jbuf, jpos, "bit_rate_bps", &
         modem%waveform%chip_rate / real(modem%waveform%spreading_factor, dp))
    call json_int (jbuf, jpos, "polar_n",              modem%coding%payload_n)
    call json_int (jbuf, jpos, "polar_k",              modem%coding%payload_k)
    call json_int (jbuf, jpos, "polar_list",           modem%coding%list_size)
    call json_real(jbuf, jpos, "exposure_weight",      modem%metric%exposure_weight)
    call json_int (jbuf, jpos, "refinement_iters",     modem%metric%refinement_iters)
    call json_end (jbuf, jpos)
    print '(a)', jbuf(:jpos)
  end subroutine cmd_info

  !> Encode a command-line bit string and write the resulting waveform.
  subroutine cmd_mod()
    type(cdss_modem) :: modem
    integer, allocatable :: bits(:)
    complex(dp), allocatable :: iq(:)
    character(len=512) :: bit_string, output, fmt
    logical :: has_bits
    bit_string = ""
    output     = "frame.cf32"
    fmt        = "auto"
    call parse_string_option("--bit-string", bit_string, has_bits)
    call parse_string_option_default("--output", output)
    call parse_string_option_default("--output-format", fmt)
    if (.not. has_bits) error stop "mod requires --bit-string"
    call bits_from_string(trim(bit_string), bits)
    call modem%build_frame(bits, iq)
    fmt = infer_format(trim(output), trim(fmt))
    call save_waveform(trim(output), trim(fmt), modem%waveform, iq)
    block
      character(len=512) :: jbuf
      integer :: jpos
      call json_begin(jbuf, jpos)
      call json_str(jbuf, jpos, "output", trim(output))
      call json_str(jbuf, jpos, "format", trim(fmt))
      call json_int(jbuf, jpos, "input_bits", size(bits))
      call json_int(jbuf, jpos, "samples", size(iq))
      call json_end(jbuf, jpos)
      print '(a)', jbuf(:jpos)
    end block
  end subroutine cmd_mod

  !> Load a waveform, attempt frame recovery, and print the decoded result.
  subroutine cmd_demod()
    type(cdss_modem) :: modem
    type(frame_recovery_result) :: frame
    complex(dp), allocatable :: iq(:)
    character(len=512) :: input, fmt
    integer :: max_payload
    input = "frame.cf32"
    fmt   = "auto"
    max_payload = 512
    call parse_string_option_default("--input", input)
    call parse_string_option_default("--input-format", fmt)
    call parse_int_option_default("--max-payload-bits", max_payload)
    fmt = infer_format(trim(input), trim(fmt))
    call load_waveform(trim(input), trim(fmt), modem%waveform, iq)
    call modem%recover_frame(iq, max_payload, frame)
    call print_frame_json(frame)
  end subroutine cmd_demod

  !> Run an end-to-end randomized transmit/channel/receive simulation.
  !!
  !! The command is intended for receiver regression runs under AWGN, carrier
  !! offset, drift, clock error, optional guard padding, and optional CW tone
  !! interference. It prints bit-error statistics and CRC status as JSON.
  subroutine cmd_simulate()
    use cdss_constants, only: eps
    use cdss_receiver, only: sync_preamble
    type(cdss_modem) :: modem
    type(channel_config) :: ch
    type(tone_interferer), allocatable :: tones(:)
    complex(dp), allocatable :: tx(:), rx(:)
    integer, allocatable :: tx_bits(:)
    type(rng_state) :: rng
    type(frame_recovery_result) :: frame
    integer :: nbits, seed, errors, i, frame_samples
    real(dp) :: noise_var, tone_frequency_hz, tone_sir_db, tone_phase_rad
    real(dp) :: tone_start_s, tone_duration_s, clock_scale
    real(dp) :: acquisition_score
    logical :: has_tone, found, has_tone_duration, sync_search, sync_ok
    integer :: known_prefix_samples, expected_sample, detected_sample, timing_error
    complex(dp) :: acquisition_peak, phase_rot
    complex(dp), allocatable :: rx_aligned(:)

    nbits = 64
    seed  = 1
    ch%snr_db = -27.0_dp
    tone_frequency_hz = 0.0_dp
    tone_sir_db = 20.0_dp
    tone_phase_rad = 0.0_dp
    tone_start_s = 0.0_dp
    tone_duration_s = 0.0_dp
    has_tone = .false.
    has_tone_duration = .false.
    sync_search = has_flag("--sync-search")
    call parse_int_option_default("--bits", nbits)
    call parse_int_option_default("--seed", seed)
    call parse_real_option_default("--snr-db", ch%snr_db)
    call parse_real_option_default("--cfo-hz", ch%cfo_hz)
    call parse_real_option_default("--drift-hz-s", ch%drift_hz_s)
    call parse_real_option_default("--clock-ppm", ch%clock_ppm)
    call parse_real_option_default("--prefix-pad-s", ch%prefix_pad_s)
    call parse_real_option_default("--suffix-pad-s", ch%suffix_pad_s)
    call parse_real_option("--tone-frequency-hz", tone_frequency_hz, found)
    has_tone = has_tone .or. found
    call parse_real_option("--tone-sir-db", tone_sir_db, found)
    has_tone = has_tone .or. found
    call parse_real_option("--tone-phase-rad", tone_phase_rad, found)
    has_tone = has_tone .or. found
    call parse_real_option("--tone-start-s", tone_start_s, found)
    has_tone = has_tone .or. found
    call parse_real_option("--tone-duration-s", tone_duration_s, found)
    has_tone = has_tone .or. found
    has_tone_duration = found

    if (has_tone .and. tone_frequency_hz == 0.0_dp) then
      tone_frequency_hz = modem%waveform%carrier_hz
    end if

    call rng_seed(rng, seed)
    call rng_bits(rng, nbits, tx_bits)

    call modem%build_frame(tx_bits, tx)
    if (has_tone) then
      allocate(tones(1))
      tones(1)%frequency_hz = tone_frequency_hz
      tones(1)%sir_db = tone_sir_db
      tones(1)%phase_rad = tone_phase_rad
      tones(1)%start_s = tone_start_s
      tones(1)%has_duration = has_tone_duration
      tones(1)%duration_s = tone_duration_s
    else
      allocate(tones(0))
    end if
    call simulate_full_channel(modem%waveform, ch, tones, tx, rng, rx, noise_var)
    known_prefix_samples = nint(ch%prefix_pad_s * modem%waveform%sample_rate)
    clock_scale = max(1.0_dp + ch%clock_ppm * 1.0e-6_dp, eps)
    expected_sample = nint(real(known_prefix_samples, dp) / clock_scale)
    detected_sample = expected_sample
    timing_error = 0
    acquisition_score = 0.0_dp
    sync_ok = .true.
    frame_samples = modem%estimate_frame_samples(nbits)

    if (sync_search) then
      call sync_preamble(modem%waveform, rx, detected_sample, acquisition_score, acquisition_peak)
      timing_error = detected_sample - expected_sample
      sync_ok = abs(timing_error) <= 4
      if (abs(acquisition_peak) > eps) then
        phase_rot = conjg(acquisition_peak / abs(acquisition_peak))
      else
        phase_rot = (1.0_dp, 0.0_dp)
      end if
      if (detected_sample >= 0 .and. size(rx) >= detected_sample + frame_samples) then
        allocate(rx_aligned(frame_samples))
        rx_aligned = rx(detected_sample + 1:detected_sample + frame_samples) * phase_rot
        call modem%recover_frame(rx_aligned, nbits + 64, frame, acq_peak=acquisition_peak * phase_rot)
        deallocate(rx_aligned)
      else if (detected_sample >= 0 .and. size(rx) > detected_sample) then
        allocate(rx_aligned(size(rx) - detected_sample))
        rx_aligned = rx(detected_sample + 1:) * phase_rot
        call modem%recover_frame(rx_aligned, nbits + 64, frame, acq_peak=acquisition_peak * phase_rot)
        deallocate(rx_aligned)
      else
        call modem%recover_frame(rx, nbits + 64, frame, acq_peak=acquisition_peak)
      end if
    else if (known_prefix_samples > 0 .and. known_prefix_samples < size(rx)) then
      allocate(rx_aligned(size(rx) - known_prefix_samples))
      rx_aligned = rx(known_prefix_samples + 1:)
      call modem%recover_frame(rx_aligned, nbits + 64, frame)
      deallocate(rx_aligned)
    else
      call modem%recover_frame(rx, nbits + 64, frame)
    end if

    errors = 0
    if (frame%payload_bit_count > 0) then
      do i = 1, min(size(tx_bits), size(frame%payload_bits))
        if (tx_bits(i) /= frame%payload_bits(i)) errors = errors + 1
      end do
    else
      errors = nbits
    end if

    block
      character(len=2048) :: jbuf
      integer :: jpos
      call json_begin(jbuf, jpos)
      call json_int (jbuf, jpos, "bits", nbits)
      call json_real(jbuf, jpos, "snr_db", ch%snr_db)
      call json_real(jbuf, jpos, "cfo_hz", ch%cfo_hz)
      call json_real(jbuf, jpos, "drift_hz_s", ch%drift_hz_s)
      call json_real(jbuf, jpos, "clock_ppm", ch%clock_ppm)
      call json_real(jbuf, jpos, "prefix_pad_s", ch%prefix_pad_s)
      call json_real(jbuf, jpos, "suffix_pad_s", ch%suffix_pad_s)
      call json_bool(jbuf, jpos, "sync_search", sync_search)
      call json_bool(jbuf, jpos, "sync_ok", sync_ok)
      call json_int (jbuf, jpos, "expected_sample", expected_sample)
      call json_int (jbuf, jpos, "detected_sample", detected_sample)
      call json_int (jbuf, jpos, "timing_error_samples", timing_error)
      call json_real(jbuf, jpos, "acquisition_score", acquisition_score)
      call json_int (jbuf, jpos, "tone_count", size(tones))
      if (has_tone) then
        call json_real(jbuf, jpos, "tone_frequency_hz", tones(1)%frequency_hz)
        call json_real(jbuf, jpos, "tone_sir_db", tones(1)%sir_db)
      end if
      call json_real(jbuf, jpos, "noise_var", noise_var)
      call json_int (jbuf, jpos, "bit_errors", errors)
      call json_real(jbuf, jpos, "ber", real(errors, dp) / real(max(nbits, 1), dp))
      call json_bool(jbuf, jpos, "crc_ok", frame%crc_ok)
      call json_int (jbuf, jpos, "refinement_iters_used", frame%refinement_iters_used)
      call json_int (jbuf, jpos, "payload_bit_count", frame%payload_bit_count)
      call json_end (jbuf, jpos)
      print '(a)', jbuf(:jpos)
    end block
  end subroutine cmd_simulate

  !> Print a frame-recovery result as a compact JSON object.
  !!
  !! @param[in] frame Decoded frame result and diagnostics.
  subroutine print_frame_json(frame)
    type(frame_recovery_result), intent(in) :: frame
    character(len=:), allocatable :: jbuf
    integer :: jpos, cap
    cap = 4096 + max(size(frame%payload_bits), 1)
    allocate(character(len=cap) :: jbuf)
    call json_begin(jbuf, jpos)
    call json_int (jbuf, jpos, "payload_bit_count", frame%payload_bit_count)
    call json_int (jbuf, jpos, "crc_expected", frame%payload_crc_expected)
    call json_int (jbuf, jpos, "crc_observed", frame%payload_crc_observed)
    call json_bool(jbuf, jpos, "crc_ok", frame%crc_ok)
    call json_int (jbuf, jpos, "refinement_iters_used", frame%refinement_iters_used)
    if (frame%payload_bit_count > 0) then
      call json_str(jbuf, jpos, "hard_bits", bits_to_string(frame%payload_bits))
    else
      call json_str(jbuf, jpos, "hard_bits", "")
    end if
    call json_end(jbuf, jpos)
    print '(a)', jbuf(:jpos)
    deallocate(jbuf)
  end subroutine print_frame_json

  !> Initialize a JSON object in a fixed-size character buffer.
  subroutine json_begin(buf, used)
    character(len=*), intent(inout) :: buf
    integer, intent(out) :: used
    used = 0
    call append(buf, used, '{')
  end subroutine json_begin

  !> Close a JSON object in a fixed-size character buffer.
  subroutine json_end(buf, used)
    character(len=*), intent(inout) :: buf
    integer, intent(inout) :: used
    call append(buf, used, '}')
  end subroutine json_end

  !> Append a JSON field name and separator.
  subroutine json_field(buf, used, name)
    character(len=*), intent(inout) :: buf
    integer, intent(inout) :: used
    character(len=*), intent(in) :: name
    if (used > 1) call append(buf, used, ',')
    call append(buf, used, '"')
    call append(buf, used, name)
    call append(buf, used, '":')
  end subroutine json_field

  !> Append a string-valued JSON field.
  subroutine json_str(buf, used, name, value)
    character(len=*), intent(inout) :: buf
    integer, intent(inout) :: used
    character(len=*), intent(in) :: name, value
    call json_field(buf, used, name)
    call append(buf, used, '"')
    call append(buf, used, value)
    call append(buf, used, '"')
  end subroutine json_str

  !> Append an integer-valued JSON field.
  subroutine json_int(buf, used, name, value)
    character(len=*), intent(inout) :: buf
    integer, intent(inout) :: used
    character(len=*), intent(in) :: name
    integer, intent(in) :: value
    character(len=32) :: tmp
    call json_field(buf, used, name)
    write(tmp, '(i0)') value
    call append(buf, used, trim(tmp))
  end subroutine json_int

  !> Append a real-valued JSON field using scientific notation.
  subroutine json_real(buf, used, name, value)
    character(len=*), intent(inout) :: buf
    integer, intent(inout) :: used
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value
    character(len=64) :: tmp
    call json_field(buf, used, name)
    write(tmp, '(es16.8)') value
    call append(buf, used, trim(adjustl(tmp)))
  end subroutine json_real

  !> Append a logical-valued JSON field.
  subroutine json_bool(buf, used, name, value)
    character(len=*), intent(inout) :: buf
    integer, intent(inout) :: used
    character(len=*), intent(in) :: name
    logical, intent(in) :: value
    call json_field(buf, used, name)
    if (value) then
      call append(buf, used, "true")
    else
      call append(buf, used, "false")
    end if
  end subroutine json_bool

  !> Append raw text to a fixed-size character buffer.
  !!
  !! The caller is responsible for ensuring that the buffer is large enough.
  subroutine append(buf, used, str)
    character(len=*), intent(inout) :: buf
    integer, intent(inout) :: used
    character(len=*), intent(in) :: str
    integer :: n
    n = len(str)
    buf(used + 1:used + n) = str
    used = used + n
  end subroutine append

  !> Parse an option that expects a following string argument.
  subroutine parse_string_option(flag, value, found)
    character(len=*), intent(in) :: flag
    character(len=*), intent(out) :: value
    logical, intent(out) :: found
    integer :: i
    character(len=512) :: arg
    found = .false.
    do i = 2, command_argument_count() - 1
      call get_command_argument(i, arg)
      if (trim(arg) == trim(flag)) then
        call get_command_argument(i + 1, value)
        found = .true.
      end if
    end do
  end subroutine parse_string_option

  !> Parse a string option and keep the existing value when absent.
  subroutine parse_string_option_default(flag, value)
    character(len=*), intent(in) :: flag
    character(len=*), intent(inout) :: value
    logical :: found
    character(len=len(value)) :: tmp
    call parse_string_option(flag, tmp, found)
    if (found) value = tmp
  end subroutine parse_string_option_default

  !> Parse an integer option and keep the existing value when absent.
  subroutine parse_int_option_default(flag, value)
    character(len=*), intent(in) :: flag
    integer, intent(inout) :: value
    character(len=128) :: text
    logical :: found
    call parse_string_option(flag, text, found)
    if (found) read(text, *) value
  end subroutine parse_int_option_default

  !> Parse a real-valued option and keep the existing value when absent.
  subroutine parse_real_option_default(flag, value)
    character(len=*), intent(in) :: flag
    real(dp), intent(inout) :: value
    character(len=128) :: text
    logical :: found
    call parse_string_option(flag, text, found)
    if (found) read(text, *) value
  end subroutine parse_real_option_default

  !> Parse a real-valued option and report whether it was present.
  subroutine parse_real_option(flag, value, found)
    character(len=*), intent(in) :: flag
    real(dp), intent(inout) :: value
    logical, intent(out) :: found
    character(len=128) :: text
    call parse_string_option(flag, text, found)
    if (found) read(text, *) value
  end subroutine parse_real_option

  !> Return true when a standalone flag appears on the command line.
  logical function has_flag(flag) result(found)
    character(len=*), intent(in) :: flag
    integer :: i
    character(len=512) :: arg
    found = .false.
    do i = 2, command_argument_count()
      call get_command_argument(i, arg)
      if (trim(arg) == trim(flag)) then
        found = .true.
        return
      end if
    end do
  end function has_flag

  !> Print CLI usage information.
  subroutine usage()
    print '(a)', 'usage: cdss <command> [options]'
    print '(a)', ''
    print '(a)', 'Commands:'
    print '(a)', '  info'
    print '(a)', '      Print modem geometry and receiver defaults as JSON.'
    print '(a)', '  mod --bit-string <bits> --output <path>'
    print '(a)', '      [--output-format auto|cf32|wav-audio|wav-iq]'
    print '(a)', '      Encode payload bits to an IQ/audio waveform.'
    print '(a)', '  sync --input <path> --output <path> [--bits N]'
    print '(a)', '      Find the preamble, phase-align the frame, and write synchronized waveform.'
    print '(a)', '  demod --input <path>'
    print '(a)', '      [--input-format auto|cf32|wav-audio|wav-iq] [--max-payload-bits N]'
    print '(a)', '      Decode a synchronized waveform and print CRC/result JSON.'
    print '(a)', '  simulate [--bits N] [--seed S] [channel options]'
    print '(a)', '      Run TX -> channel -> RX and print BER/CRC JSON.'
    print '(a)', ''
    print '(a)', 'Simulation channel options:'
    print '(a)', '  --snr-db X                 AWGN SNR in dB (default: -27)'
    print '(a)', '  --cfo-hz X                 carrier frequency offset in Hz'
    print '(a)', '  --drift-hz-s X             linear frequency drift in Hz/s'
    print '(a)', '  --clock-ppm X              sample clock error in ppm'
    print '(a)', '  --prefix-pad-s X           leading silence before the frame'
    print '(a)', '  --suffix-pad-s X           trailing silence after the frame'
    print '(a)', '  --sync-search              validate preamble acquisition before decoding'
    print '(a)', '  --tone-frequency-hz X      add one CW interferer at this frequency'
    print '(a)', '  --tone-sir-db X            tone signal-to-interference ratio in dB'
    print '(a)', '  --tone-phase-rad X         tone initial phase in radians'
    print '(a)', '  --tone-start-s X           tone start time in seconds'
    print '(a)', '  --tone-duration-s X        finite tone duration in seconds'
  end subroutine usage
end program cdss_cli
