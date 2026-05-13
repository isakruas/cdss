!> @file test_io_roundtrip.f90
!! @brief Round-trip tests for supported waveform file formats.
program test_io_roundtrip
  use cdss_kinds, only: dp
  use cdss_types, only: waveform_config
  use cdss_io, only: save_waveform, load_waveform, infer_format
  implicit none

  type(waveform_config) :: cfg
  complex(dp), allocatable :: iq(:), loaded(:)
  character(len=*), parameter :: cf32_path = "/tmp/cdss_test_io_roundtrip.cf32"
  character(len=*), parameter :: waviq_path = "/tmp/cdss_test_io_roundtrip_iq.wav"
  character(len=*), parameter :: wavaudio_path = "/tmp/cdss_test_io_roundtrip_audio.wav"
  real(dp) :: err
  integer :: i

  print '(a)', "test_io_roundtrip: starting waveform I/O checks"
  cfg%sample_rate = 8000.0_dp
  allocate(iq(128))
  do i = 1, size(iq)
    iq(i) = cmplx(0.5_dp * sin(0.03_dp * real(i, dp)), &
                  0.25_dp * cos(0.07_dp * real(i, dp)), dp)
  end do

  print '(a,a)', "  infer .cf32 -> ", trim(infer_format("x.cf32", "auto"))
  print '(a,a)', "  infer .wav -> ", trim(infer_format("x.wav", "auto"))
  if (trim(infer_format("x.cf32", "auto")) /= "cf32") error stop 1
  if (trim(infer_format("x.wav", "auto")) /= "wav-audio") error stop 2
  if (trim(infer_format("x.bin", "wav-iq")) /= "wav-iq") error stop 3

  call save_waveform(cf32_path, "cf32", cfg, iq)
  call load_waveform(cf32_path, "cf32", cfg, loaded)
  err = maxval(abs(loaded - iq))
  print '(a,es12.4)', "  cf32 max roundtrip error=", err
  if (size(loaded) /= size(iq) .or. err > 1.0e-6_dp) then
    print '(a)', "test_io_roundtrip: cf32 roundtrip failed"
    error stop 4
  end if

  call save_waveform(waviq_path, "wav-iq", cfg, iq)
  call load_waveform(waviq_path, "wav-iq", cfg, loaded)
  err = maxval(abs(loaded - iq))
  print '(a,es12.4)', "  wav-iq max roundtrip error=", err
  if (size(loaded) /= size(iq) .or. err > 5.0e-5_dp) then
    print '(a)', "test_io_roundtrip: wav-iq roundtrip failed"
    error stop 5
  end if

  call save_waveform(wavaudio_path, "wav-audio", cfg, iq)
  call load_waveform(wavaudio_path, "wav-audio", cfg, loaded)
  err = maxval(abs(real(loaded, dp) - real(iq, dp)))
  print '(a,es12.4,a,es12.4)', "  wav-audio real error=", err, &
       ", max imag=", maxval(abs(aimag(loaded)))
  if (size(loaded) /= size(iq) .or. err > 5.0e-5_dp) then
    print '(a)', "test_io_roundtrip: wav-audio real component roundtrip failed"
    error stop 6
  end if
  if (maxval(abs(aimag(loaded))) > 1.0e-12_dp) then
    print '(a)', "test_io_roundtrip: wav-audio should load with zero quadrature"
    error stop 7
  end if

  call delete_file(cf32_path)
  call delete_file(waviq_path)
  call delete_file(wavaudio_path)
  print '(a)', "test_io_roundtrip: passed"
  deallocate(iq, loaded)

contains

  subroutine delete_file(path)
    character(len=*), intent(in) :: path
    integer :: unit, ios
    open(newunit=unit, file=path, status='old', iostat=ios)
    if (ios == 0) close(unit, status='delete')
  end subroutine delete_file
end program test_io_roundtrip
