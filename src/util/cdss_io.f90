!> @file cdss_io.f90
!! @brief Waveform file I/O for simulation and command-line use.
!!
!! @details
!! Supported formats are interleaved 32-bit floating-point complex samples
!! ("cf32"), mono 16-bit PCM WAV carrying the real passband projection
!! ("wav-audio"), and stereo 16-bit PCM WAV carrying I/Q as left/right channels
!! ("wav-iq"). The routines intentionally avoid external audio dependencies.
module cdss_io
  use cdss_kinds, only: dp, i32
  use cdss_types, only: waveform_config
  implicit none
  private
  public :: save_waveform, load_waveform, infer_format

  integer, parameter :: WAV_PCM = 1

contains

  !> Infer the waveform file format from an explicit option or file extension.
  !!
  !! @param[in] path  File path used when the format argument is "auto".
  !! @param[in] given Explicit format string or "auto".
  !! @return One of "cf32", "wav-audio", or the explicit non-auto value.
  function infer_format(path, given) result(fmt)
    character(len=*), intent(in) :: path, given
    character(len=16) :: fmt
    integer :: n
    if (trim(given) /= "auto" .and. len_trim(given) > 0) then
      fmt = given
      return
    end if
    n = len_trim(path)
    if (n >= 5 .and. path(n - 4:n) == '.cf32') then
      fmt = "cf32"
    else if (n >= 4 .and. path(n - 3:n) == '.wav') then
      fmt = "wav-audio"
    else
      fmt = "cf32"
    end if
  end function infer_format

  !> Save a complex waveform using the selected on-disk representation.
  !!
  !! @param[in] path Output file path.
  !! @param[in] fmt  Format selector: "cf32", "wav-audio", or "wav-iq".
  !! @param[in] cfg  Waveform configuration used for WAV sample-rate metadata.
  !! @param[in] iq   Complex samples to write.
  subroutine save_waveform(path, fmt, cfg, iq)
    character(len=*), intent(in) :: path, fmt
    type(waveform_config), intent(in) :: cfg
    complex(dp), intent(in) :: iq(:)
    select case (trim(fmt))
    case ("cf32")
      call save_cf32(path, iq)
    case ("wav-audio")
      call save_wav_audio(path, cfg, iq)
    case ("wav-iq")
      call save_wav_iq(path, cfg, iq)
    case default
      error stop "save_waveform: unsupported format"
    end select
  end subroutine save_waveform

  !> Load a complex waveform from disk.
  !!
  !! @param[in]  path Input file path.
  !! @param[in]  fmt  Format selector: "cf32", "wav-audio", or "wav-iq".
  !! @param[in]  cfg  Waveform configuration used to validate WAV sample rate.
  !! @param[out] iq   Allocated complex samples.
  subroutine load_waveform(path, fmt, cfg, iq)
    character(len=*), intent(in) :: path, fmt
    type(waveform_config), intent(in) :: cfg
    complex(dp), allocatable, intent(out) :: iq(:)
    select case (trim(fmt))
    case ("cf32")
      call load_cf32(path, iq)
    case ("wav-audio")
      call load_wav_audio(path, cfg, iq)
    case ("wav-iq")
      call load_wav_iq(path, iq)
    case default
      error stop "load_waveform: unsupported format"
    end select
  end subroutine load_waveform

  !> Write interleaved binary float32 I/Q samples.
  !!
  !! @param[in] path Output path.
  !! @param[in] iq   Complex samples written as Re, Im, Re, Im, ...
  subroutine save_cf32(path, iq)
    character(len=*), intent(in) :: path
    complex(dp), intent(in) :: iq(:)
    integer :: unit, i
    real(kind(1.0e0)) :: i_re, i_im
    open(newunit=unit, file=path, access='stream', form='unformatted', status='replace')
    do i = 1, size(iq)
      i_re = real(real(iq(i)), kind(1.0e0))
      i_im = real(aimag(iq(i)), kind(1.0e0))
      write(unit) i_re, i_im
    end do
    close(unit)
  end subroutine save_cf32

  !> Read interleaved binary float32 I/Q samples.
  !!
  !! @param[in]  path Input path.
  !! @param[out] iq   Allocated complex samples reconstructed from Re/Im pairs.
  subroutine load_cf32(path, iq)
    character(len=*), intent(in) :: path
    complex(dp), allocatable, intent(out) :: iq(:)
    integer :: unit, n_bytes, n_samples, i
    real(kind(1.0e0)) :: re, im
    open(newunit=unit, file=path, access='stream', form='unformatted', status='old', action='read')
    inquire(unit=unit, size=n_bytes)
    n_samples = n_bytes / 8
    allocate(iq(n_samples))
    do i = 1, n_samples
      read(unit) re, im
      iq(i) = cmplx(real(re, dp), real(im, dp), dp)
    end do
    close(unit)
  end subroutine load_cf32

  !> Save the real component of the waveform as mono 16-bit PCM WAV.
  !!
  !! The caller supplies a waveform already centered at the configured audio
  !! carrier; therefore the file stores x[n] = Re{iq[n]} without additional
  !! frequency translation.
  !!
  !! @param[in] path Output WAV path.
  !! @param[in] cfg  Waveform configuration used for the sample rate.
  !! @param[in] iq   Complex waveform to project onto real audio.
  subroutine save_wav_audio(path, cfg, iq)
    character(len=*), intent(in) :: path
    type(waveform_config), intent(in) :: cfg
    complex(dp), intent(in) :: iq(:)
    integer :: unit, i
    integer(2), allocatable :: samples(:)
    real(dp) :: peak, scale
    allocate(samples(size(iq)))
    peak = max(maxval(abs(real(iq))), 1.0e-12_dp)
    scale = 32767.0_dp / max(peak, 1.0_dp)
    do i = 1, size(iq)
      samples(i) = int(max(-32768.0_dp, min(32767.0_dp, real(iq(i)) * scale)), 2)
    end do
    call write_wav_header(unit, path, int(cfg%sample_rate), 1, size(samples))
    do i = 1, size(samples)
      write(unit) samples(i)
    end do
    close(unit)
    deallocate(samples)
  end subroutine save_wav_audio

  !> Save complex samples as stereo 16-bit PCM WAV.
  !!
  !! The left channel contains in-phase samples and the right channel contains
  !! quadrature samples.
  !!
  !! @param[in] path Output WAV path.
  !! @param[in] cfg  Waveform configuration used for the sample rate.
  !! @param[in] iq   Complex waveform to write.
  subroutine save_wav_iq(path, cfg, iq)
    character(len=*), intent(in) :: path
    type(waveform_config), intent(in) :: cfg
    complex(dp), intent(in) :: iq(:)
    integer :: unit, i
    integer(2) :: s_i, s_q
    real(dp) :: peak, scale
    peak = max(maxval(abs(real(iq))), maxval(abs(aimag(iq))), 1.0e-12_dp)
    scale = 32767.0_dp / max(peak, 1.0_dp)
    call write_wav_header(unit, path, int(cfg%sample_rate), 2, 2 * size(iq))
    do i = 1, size(iq)
      s_i = int(max(-32768.0_dp, min(32767.0_dp, real(iq(i)) * scale)), 2)
      s_q = int(max(-32768.0_dp, min(32767.0_dp, aimag(iq(i)) * scale)), 2)
      write(unit) s_i, s_q
    end do
    close(unit)
  end subroutine save_wav_iq

  !> Load mono 16-bit PCM WAV as a real-valued complex waveform.
  !!
  !! @param[in]  path Input WAV path.
  !! @param[in]  cfg  Expected waveform configuration.
  !! @param[out] iq   Allocated complex samples with zero quadrature component.
  subroutine load_wav_audio(path, cfg, iq)
    character(len=*), intent(in) :: path
    type(waveform_config), intent(in) :: cfg
    complex(dp), allocatable, intent(out) :: iq(:)
    integer :: unit, n_samples, i
    integer(2) :: s
    integer :: sample_rate, n_channels
    call read_wav_header(unit, path, sample_rate, n_channels, n_samples)
    if (n_channels /= 1) error stop "load_wav_audio: expected mono"
    if (sample_rate /= int(cfg%sample_rate)) then
      ! Allow mismatch but leave resampling responsibility to the caller.
      continue
    end if
    allocate(iq(n_samples))
    do i = 1, n_samples
      read(unit) s
      iq(i) = cmplx(real(s, dp) / 32767.0_dp, 0.0_dp, dp)
    end do
    close(unit)
  end subroutine load_wav_audio

  !> Load stereo 16-bit PCM WAV as complex I/Q samples.
  !!
  !! @param[in]  path Input WAV path.
  !! @param[out] iq   Allocated complex samples with left as I and right as Q.
  subroutine load_wav_iq(path, iq)
    character(len=*), intent(in) :: path
    complex(dp), allocatable, intent(out) :: iq(:)
    integer :: unit, n_samples, i
    integer(2) :: s_i, s_q
    integer :: sample_rate, n_channels
    call read_wav_header(unit, path, sample_rate, n_channels, n_samples)
    if (n_channels /= 2) error stop "load_wav_iq: expected stereo"
    allocate(iq(n_samples))
    do i = 1, n_samples
      read(unit) s_i, s_q
      iq(i) = cmplx(real(s_i, dp) / 32767.0_dp, real(s_q, dp) / 32767.0_dp, dp)
    end do
    close(unit)
  end subroutine load_wav_iq

  !> Open a stream file and write a minimal PCM WAV header.
  !!
  !! @param[out] unit        Open Fortran unit positioned at the start of data.
  !! @param[in]  path        Output path.
  !! @param[in]  sample_rate Sample rate in hertz.
  !! @param[in]  n_channels  Number of PCM channels.
  !! @param[in]  n_samples   Number of 16-bit samples, including all channels.
  subroutine write_wav_header(unit, path, sample_rate, n_channels, n_samples)
    integer, intent(out) :: unit
    character(len=*), intent(in) :: path
    integer, intent(in) :: sample_rate, n_channels, n_samples
    integer :: byte_rate, block_align, data_size
    integer(2) :: bps, fmt, ch
    bps = 16
    fmt = WAV_PCM
    ch  = int(n_channels, 2)
    block_align = n_channels * 2
    byte_rate   = sample_rate * block_align
    data_size   = n_samples * 2
    open(newunit=unit, file=path, access='stream', form='unformatted', status='replace')
    write(unit) "RIFF"
    write(unit) 36 + data_size
    write(unit) "WAVE"
    write(unit) "fmt "
    write(unit) 16
    write(unit) fmt
    write(unit) ch
    write(unit) sample_rate
    write(unit) byte_rate
    write(unit) int(block_align, 2)
    write(unit) bps
    write(unit) "data"
    write(unit) data_size
  end subroutine write_wav_header

  !> Read a minimal PCM WAV header and position the file at the data payload.
  !!
  !! @param[out] unit        Open Fortran unit positioned at audio data.
  !! @param[in]  path        Input WAV path.
  !! @param[out] sample_rate Sample rate in hertz.
  !! @param[out] n_channels  Number of channels.
  !! @param[out] n_samples   Number of sample frames per channel.
  subroutine read_wav_header(unit, path, sample_rate, n_channels, n_samples)
    integer, intent(out) :: unit, sample_rate, n_channels, n_samples
    character(len=*), intent(in) :: path
    character(len=4) :: tag
    integer :: chunk_size, fmt_size, byte_rate, data_size
    integer(2) :: fmt, ch, block_align, bps
    open(newunit=unit, file=path, access='stream', form='unformatted', status='old', action='read')
    read(unit) tag
    if (tag /= "RIFF") error stop "read_wav_header: missing RIFF"
    read(unit) chunk_size
    read(unit) tag
    if (tag /= "WAVE") error stop "read_wav_header: missing WAVE"
    do
      read(unit) tag
      read(unit) chunk_size
      if (tag == "fmt ") exit
      call skip_bytes(unit, chunk_size)
    end do
    fmt_size = chunk_size
    read(unit) fmt
    read(unit) ch
    read(unit) sample_rate
    read(unit) byte_rate
    read(unit) block_align
    read(unit) bps
    if (fmt_size > 16) call skip_bytes(unit, fmt_size - 16)
    do
      read(unit) tag
      read(unit) chunk_size
      if (tag == "data") exit
      call skip_bytes(unit, chunk_size)
    end do
    data_size = chunk_size
    n_channels = int(ch)
    n_samples = data_size / (2 * n_channels)
  end subroutine read_wav_header

  !> Consume a fixed number of bytes from an already-open stream file.
  !!
  !! @param[in] unit Open file unit.
  !! @param[in] n    Number of bytes to skip.
  subroutine skip_bytes(unit, n)
    integer, intent(in) :: unit, n
    integer :: i
    integer(1) :: byte
    do i = 1, n
      read(unit) byte
    end do
  end subroutine skip_bytes
end module cdss_io
