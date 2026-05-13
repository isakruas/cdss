!> @file test_payload_edges.f90
!! @brief Edge-case tests for payload lengths and Polar block boundaries.
program test_payload_edges
  use cdss_kinds, only: dp
  use cdss_types, only: waveform_config, soft_metric_config, coding_config, &
                          frame_recovery_result
  use cdss_receiver, only: encode_frame, recover_frame, frame_n_bits, frame_n_samples
  use cdss_waveform, only: synthesize_waveform
  implicit none

  type(waveform_config) :: cfg
  type(soft_metric_config) :: metric
  type(coding_config) :: coding

  print '(a)', "test_payload_edges: starting payload boundary checks"
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
  metric%refinement_iters = 1

  call check_payload_length(0)
  call check_payload_length(1)
  call check_payload_length(coding%payload_k - 1)
  call check_payload_length(coding%payload_k)
  call check_payload_length(coding%payload_k + 1)
  call check_payload_length(2 * coding%payload_k + 2)

  print '(a)', "test_payload_edges: passed"

contains

  subroutine check_payload_length(n_bits)
    integer, intent(in) :: n_bits
    integer, allocatable :: payload(:), coded_bits(:)
    complex(dp), allocatable :: iq(:)
    type(frame_recovery_result) :: frame
    integer :: i, expected_bits, expected_samples

    allocate(payload(n_bits))
    do i = 1, n_bits
      payload(i) = merge(1, 0, mod(3 * i + 1, 5) < 2)
    end do

    call encode_frame(cfg, coding, payload, coded_bits)
    call synthesize_waveform(cfg, coded_bits, iq)
    expected_bits = cfg%preamble_bits + frame_n_bits(coding, n_bits)
    expected_samples = frame_n_samples(cfg, coding, n_bits)
    print '(a,i0,a,i0,a,i0)', "  payload=", n_bits, ", coded bits=", &
         size(coded_bits), ", samples=", size(iq)

    if (size(coded_bits) /= expected_bits) then
      print '(a)', "test_payload_edges: coded-bit count mismatch"
      error stop 1
    end if
    if (size(iq) /= expected_samples) then
      print '(a)', "test_payload_edges: sample count mismatch"
      error stop 2
    end if

    call recover_frame(cfg, metric, coding, iq, max(n_bits + coding%payload_k, 1), frame)
    print '(a,l1,a,i0)', "    recovery crc_ok=", frame%crc_ok, &
         ", decoded bits=", frame%payload_bit_count
    if (.not. frame%crc_ok) then
      print '(a,i0)', "test_payload_edges: CRC failed for payload length ", n_bits
      error stop 3
    end if
    if (frame%payload_bit_count /= n_bits) then
      print '(a)', "test_payload_edges: decoded payload length mismatch"
      error stop 4
    end if
    if (n_bits > 0) then
      if (any(frame%payload_bits /= payload)) then
        print '(a,i0)', "test_payload_edges: payload mismatch for length ", n_bits
        error stop 5
      end if
    end if

    deallocate(payload, coded_bits, iq)
  end subroutine check_payload_length
end program test_payload_edges
