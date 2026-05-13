!> @file cdss_api.f90
!! @brief Public facade for constructing, modulating, and demodulating CDSS frames.
!!
!! @details
!! This module provides a compact object-oriented facade over the lower-level
!! waveform, channel-coding, and receiver modules. It is intended for command-line
!! tools, tests, and external applications that need a stable entry point without
!! managing each processing stage directly.
module cdss_api
  use cdss_kinds, only: dp
  use cdss_types, only: waveform_config, soft_metric_config, coding_config, &
                          frame_recovery_result
  use cdss_waveform, only: synthesize_waveform
  use cdss_receiver, only: encode_frame, recover_frame, frame_n_samples
  implicit none
  private
  public :: cdss_modem

  !> High-level modem configuration and processing facade.
  !!
  !! The object stores the waveform geometry, soft-demodulator tuning, and Polar
  !! coding parameters used consistently by both the transmitter and receiver.
  type :: cdss_modem
    type(waveform_config)    :: waveform  !< Physical waveform and spreading parameters.
    type(soft_metric_config) :: metric    !< Soft-demodulation and uncertainty settings.
    type(coding_config)      :: coding    !< Polar coding and interleaver parameters.
  contains
    procedure :: build_frame
    procedure :: recover_frame  => modem_recover
    procedure :: estimate_frame_samples
  end type cdss_modem

contains

  !> Encode payload bits into a complete complex baseband frame.
  !!
  !! @param[in]  self         Modem instance containing waveform and coding settings.
  !! @param[in]  payload_bits Uncoded payload bits, each represented as 0 or 1.
  !! @param[out] iq           Allocated IQ samples for preamble, coded header, and
  !!                          coded payload.
  subroutine build_frame(self, payload_bits, iq)
    class(cdss_modem), intent(in) :: self
    integer, intent(in) :: payload_bits(:)
    complex(dp), allocatable, intent(out) :: iq(:)
    integer, allocatable :: all_coded_bits(:)
    call encode_frame(self%waveform, self%coding, payload_bits, all_coded_bits)
    call synthesize_waveform(self%waveform, all_coded_bits, iq)
    deallocate(all_coded_bits)
  end subroutine build_frame

  !> Recover a payload frame from received complex baseband samples.
  !!
  !! @param[in]  self             Modem instance containing receiver settings.
  !! @param[in]  rx               Received IQ samples. The current receiver path
  !!                              assumes the frame is already sample-aligned.
  !! @param[in]  max_payload_hint Safety bound for decoded payload length in bits.
  !! @param[out] result           Decoding result, including CRC status and metrics.
  !! @param[in]  acq_peak         Optional complex correlation peak from preamble sync.
  subroutine modem_recover(self, rx, max_payload_hint, result, acq_peak)
    class(cdss_modem), intent(in) :: self
    complex(dp), intent(in) :: rx(:)
    integer, intent(in) :: max_payload_hint
    type(frame_recovery_result), intent(out) :: result
    complex(dp), intent(in), optional :: acq_peak
    call recover_frame(self%waveform, self%metric, self%coding, rx, max_payload_hint, result, acq_peak)
  end subroutine modem_recover

  !> Estimate the number of IQ samples occupied by a frame.
  !!
  !! @param[in] self              Modem instance containing waveform and coding settings.
  !! @param[in] payload_bit_count Number of uncoded payload bits.
  !! @return Total frame length in samples, including preamble and coding overhead.
  pure integer function estimate_frame_samples(self, payload_bit_count) result(n)
    class(cdss_modem), intent(in) :: self
    integer, intent(in) :: payload_bit_count
    n = frame_n_samples(self%waveform, self%coding, payload_bit_count)
  end function estimate_frame_samples
end module cdss_api
