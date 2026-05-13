!> @file cdss_uncertainty.f90
!! @brief STFT-based per-chip interference exposure estimation.
!!
!! @details
!! The soft demodulator needs a reliability estimate for each PN chip. This
!! module computes a power spectrogram, estimates a global out-of-band floor, and
!! measures excess energy along the expected chirp trajectory. The resulting
!! normalized exposure values are clipped and converted to chip weights by the
!! demodulator.
module cdss_uncertainty
  use cdss_kinds, only: dp
  use cdss_constants, only: pi, two_pi, eps
  use cdss_types, only: waveform_config, soft_metric_config
  use cdss_fft, only: fft_stft
  implicit none
  private
  public :: compute_chip_uncertainty

contains

  !> Estimate normalized interference exposure for each received chip.
  !!
  !! @param[in]  cfg    Waveform configuration.
  !! @param[in]  metric STFT and exposure-clipping parameters.
  !! @param[in]  rx     Complex samples to analyze.
  !! @param[out] chip_J Allocated exposure vector with one value per full chip.
  !!
  !! Values near zero indicate no detectable excess energy above the estimated
  !! floor. Larger values indicate that the chip's time-frequency location is
  !! contaminated by narrow-band or impulsive interference.
  subroutine compute_chip_uncertainty(cfg, metric, rx, chip_J)
    type(waveform_config),    intent(in) :: cfg
    type(soft_metric_config), intent(in) :: metric
    complex(dp),              intent(in) :: rx(:)
    real(dp), allocatable,    intent(out) :: chip_J(:)

    real(dp), allocatable :: spectrogram(:, :), bin_floor(:)
    real(dp) :: df, fs, dt_per_frame, t_chip, t_in_bit, f_inst, floor_global
    integer  :: nfft, hop, n_frames, total_chips, k, frame_idx, bin_idx, lobe, idx
    real(dp) :: sum_lobe, denom, norm
    real(dp), allocatable :: outband_pool(:)
    integer :: n_outband, total_samples, samples_per_chip, chips_per_bit
    real(dp) :: band_hz, signal_center_hz

    total_samples    = size(rx)
    samples_per_chip = cfg%samples_per_chip()
    chips_per_bit    = cfg%spreading_factor
    total_chips      = total_samples / samples_per_chip

    if (total_chips <= 0) then
      allocate(chip_J(0))
      return
    end if

    nfft = metric%stft_window
    hop  = metric%stft_hop
    fs   = cfg%sample_rate
    df   = fs / real(nfft, dp)

    call fft_stft(rx, nfft, hop, spectrogram, n_frames)
    if (n_frames == 0) then
      allocate(chip_J(total_chips))
      chip_J = 0.0_dp
      return
    end if

    ! Per-bin floor: mean across frames approximates stationary background per bin.
    allocate(bin_floor(nfft))
    do bin_idx = 1, nfft
      bin_floor(bin_idx) = sum(spectrogram(:, bin_idx)) / real(n_frames, dp)
    end do

    ! Out-of-band global floor: bins whose freq is outside the signal band.
    signal_center_hz = cfg%carrier_hz
    band_hz          = cfg%bandwidth_hz() + metric%outband_guard_hz
    allocate(outband_pool(nfft))
    n_outband = 0
    do bin_idx = 1, nfft
      if (bin_idx - 1 <= nfft / 2) then
        f_inst = real(bin_idx - 1, dp) * df
      else
        f_inst = real(bin_idx - 1 - nfft, dp) * df
      end if
      if (abs(f_inst - signal_center_hz) > band_hz .and. &
          abs(f_inst + signal_center_hz) > band_hz .and. &
          abs(f_inst) > band_hz) then
        n_outband = n_outband + 1
        outband_pool(n_outband) = bin_floor(bin_idx)
      end if
    end do
    if (n_outband > 0) then
      floor_global = median_of(outband_pool(1:n_outband))
    else
      floor_global = max(sum(bin_floor) / real(nfft, dp), eps)
    end if
    deallocate(outband_pool)

    dt_per_frame = real(hop, dp) / fs

    allocate(chip_J(total_chips))
    do k = 1, total_chips
      ! Map absolute chip center time to the nearest STFT frame index.
      t_chip = (real(k - 1, dp) + 0.5_dp) * real(samples_per_chip, dp) / fs
      frame_idx = nint(t_chip / dt_per_frame) + 1
      frame_idx = max(1, min(n_frames, frame_idx))

      ! The chirp resets each bit, so frequency must be evaluated relative to the
      ! current bit rather than absolute frame time.
      t_in_bit = (real(mod(k - 1, chips_per_bit), dp) + 0.5_dp) * real(samples_per_chip, dp) / fs
      f_inst = cfg%carrier_hz + cfg%chirp_slope_hz_s * &
               (t_in_bit - real(chips_per_bit, dp) * 0.5_dp / cfg%chip_rate)

      ! Map to a wrapped FFT bin. Negative frequencies wrap above nfft/2.
      bin_idx = nint(modulo(f_inst / df + real(nfft, dp), real(nfft, dp))) + 1
      bin_idx = max(1, min(nfft, bin_idx))

      ! Average the target bin and its neighbors to absorb spectral leakage from
      ! finite-window STFT analysis and raised-cosine chip tapering.
      sum_lobe = 0.0_dp
      do lobe = -1, 1
        idx = modulo(bin_idx - 1 + lobe, nfft) + 1
        sum_lobe = sum_lobe + max(spectrogram(frame_idx, idx) - floor_global, 0.0_dp)
      end do
      sum_lobe = sum_lobe / 3.0_dp

      denom = floor_global + eps
      norm = sum_lobe / denom
      chip_J(k) = max(0.0_dp, min(norm, metric%exposure_clip))
    end do

    deallocate(spectrogram, bin_floor)
  end subroutine compute_chip_uncertainty

  !> Compute the median of a real-valued vector.
  !!
  !! @param[in] values Input values.
  !! @return Median value, or zero for an empty array.
  pure function median_of(values) result(med)
    real(dp), intent(in) :: values(:)
    real(dp) :: med
    real(dp), allocatable :: sorted(:)
    integer :: n
    n = size(values)
    if (n == 0) then
      med = 0.0_dp
      return
    end if
    allocate(sorted(n))
    sorted = values
    call insertion_sort(sorted)
    if (mod(n, 2) == 0) then
      med = 0.5_dp * (sorted(n / 2) + sorted(n / 2 + 1))
    else
      med = sorted(n / 2 + 1)
    end if
    deallocate(sorted)
  end function median_of

  !> Sort a real-valued array in ascending order using insertion sort.
  !!
  !! @param[inout] a Array to sort in place.
  !!
  !! The arrays sorted here are small STFT summary pools, so insertion sort keeps
  !! the dependency footprint minimal without affecting receiver complexity.
  pure subroutine insertion_sort(a)
    real(dp), intent(inout) :: a(:)
    integer :: i, j
    real(dp) :: key
    do i = 2, size(a)
      key = a(i)
      j = i - 1
      do while (j >= 1)
        if (a(j) <= key) exit
        a(j + 1) = a(j)
        j = j - 1
      end do
      a(j + 1) = key
    end do
  end subroutine insertion_sort
end module cdss_uncertainty
