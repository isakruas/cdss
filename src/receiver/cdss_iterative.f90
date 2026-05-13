!> @file cdss_iterative.f90
!! @brief Feedback-aided soft-demodulation refinement.
!!
!! @details
!! The refinement step re-synthesizes the current decoded coded-bit estimate and
!! subtracts it from the received waveform. The resulting residue is used only to
!! refresh the STFT-based interference exposure map, while demodulation remains
!! referenced to the original received samples. This reduces self-contamination
!! of the exposure estimate when the current decode is correct or nearly correct.
module cdss_iterative
  use cdss_kinds, only: dp
  use cdss_constants, only: eps
  use cdss_types, only: waveform_config, soft_metric_config, soft_demod_result
  use cdss_waveform, only: synthesize_waveform
  use cdss_softdemod, only: soft_demod_bits
  implicit none
  private
  public :: refine_soft_demod

contains

  !> Recompute soft metrics using a decoded-codeword residue.
  !!
  !! @param[in]  cfg                Waveform configuration.
  !! @param[in]  metric             Soft-demodulator configuration.
  !! @param[in]  rx                 Received samples aligned to the payload.
  !! @param[in]  n_bits             Number of coded bits represented by rx.
  !! @param[in]  decoded_coded_bits Current best coded-bit estimate.
  !! @param[out] refined            Soft-demodulation result using refreshed
  !!                                chip exposure metrics.
  !! @param[inout] state            Optional tracking state.
  subroutine refine_soft_demod(cfg, metric, rx, n_bits, decoded_coded_bits, refined, state)
    use cdss_types, only: tracking_state
    type(waveform_config),    intent(in)  :: cfg
    type(soft_metric_config), intent(in)  :: metric
    complex(dp),              intent(in)  :: rx(:)
    integer,                  intent(in)  :: n_bits
    integer,                  intent(in)  :: decoded_coded_bits(:)   ! current best estimate
    type(soft_demod_result),  intent(out) :: refined
    type(tracking_state), intent(inout), optional :: state

    complex(dp), allocatable :: tx_estimate(:), residue(:)
    integer :: total_samples, i

    total_samples = n_bits * cfg%samples_per_bit()
    if (size(rx) < total_samples) then
      error stop "refine_soft_demod: rx shorter than expected"
    end if
    if (size(decoded_coded_bits) /= n_bits) then
      error stop "refine_soft_demod: decoded_coded_bits length mismatch"
    end if

    ! Re-synthesize using the decoded coded bits.
    call synthesize_waveform(cfg, decoded_coded_bits, tx_estimate)

    ! Subtraction approximates a noise-plus-interference residue when the decode
    ! is correct.
    allocate(residue(total_samples))
    do i = 1, total_samples
      residue(i) = rx(i) - tx_estimate(i)
    end do

    ! Re-run demod on the SAME rx but passing the residue solely to refresh the
    ! STFT exposure estimate. This enables the receiver to "see past" its own
    ! signal and accurately map the underlying interference floor.
    block
      type(soft_metric_config) :: metric_iter
      real(dp) :: residue_power
      metric_iter = metric
      residue_power = sum(abs(residue)**2) / real(total_samples, dp)
      metric_iter%has_noise_var = .true.
      metric_iter%noise_var = max(residue_power, eps)
      call soft_demod_bits(cfg, metric_iter, rx, n_bits, refined, state=state, bg_residue=residue)
    end block

    deallocate(tx_estimate, residue)
  end subroutine refine_soft_demod
end module cdss_iterative
