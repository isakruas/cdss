!> @file test_noisy_integration.f90
!! @brief High-SNR integration test through modem, channel, and receiver.
program test_noisy_integration
  use cdss_kinds, only: dp
  use cdss_api, only: cdss_modem
  use cdss_types, only: channel_config, tone_interferer, frame_recovery_result
  use cdss_channel, only: simulate_full_channel
  use cdss_rng, only: rng_state, rng_seed, rng_bits
  implicit none

  type(cdss_modem) :: modem
  type(channel_config) :: ch
  type(tone_interferer), allocatable :: tones(:)
  type(rng_state) :: rng
  type(frame_recovery_result) :: frame
  integer, allocatable :: tx_bits(:)
  complex(dp), allocatable :: tx(:), rx(:)
  real(dp) :: noise_var
  integer :: errors, i

  print '(a)', "test_noisy_integration: starting high-SNR channel integration check"
  allocate(tones(0))
  ch%snr_db = 35.0_dp
  ch%cfo_hz = 0.0_dp
  ch%drift_hz_s = 0.0_dp
  ch%clock_ppm = 0.0_dp

  call rng_seed(rng, 2026)
  call rng_bits(rng, 64, tx_bits)
  print '(a,i0,a,i0)', "  payload bits=", size(tx_bits), ", ones=", count(tx_bits == 1)

  call modem%build_frame(tx_bits, tx)
  print '(a,i0)', "  TX samples=", size(tx)

  call simulate_full_channel(modem%waveform, ch, tones, tx, rng, rx, noise_var)
  print '(a,es12.4,a,i0)', "  channel noise_var=", noise_var, ", RX samples=", size(rx)

  call modem%recover_frame(rx, size(tx_bits) + 64, frame)
  errors = size(tx_bits)
  if (frame%payload_bit_count == size(tx_bits)) then
    errors = 0
    do i = 1, size(tx_bits)
      if (tx_bits(i) /= frame%payload_bits(i)) errors = errors + 1
    end do
  end if

  print '(a,l1,a,i0,a,i0)', "  crc_ok=", frame%crc_ok, ", decoded bits=", &
       frame%payload_bit_count, ", bit errors=", errors
  print '(a,i0,a,i0)', "  CRC expected=", frame%payload_crc_expected, &
       ", observed=", frame%payload_crc_observed
  if (.not. frame%crc_ok) then
    print '(a)', "test_noisy_integration: high-SNR frame failed CRC"
    error stop 1
  end if
  if (errors /= 0) then
    print '(a)', "test_noisy_integration: high-SNR frame has bit errors"
    error stop 2
  end if

  print '(a)', "test_noisy_integration: passed"
  deallocate(tones, tx_bits, tx, rx)
end program test_noisy_integration
