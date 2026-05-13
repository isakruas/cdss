!> @file cdss_types.f90
!! @brief Core data structures and configuration types for the CDSS modem.
!! @details
!! The types declared here are the stable contracts between transmitter,
!! channel model, soft receiver, and application code. Field documentation is
!! intentionally explicit because these values define the experimental geometry
!! reported in simulations and publications.

module cdss_types
  use cdss_kinds, only: dp, i32
  use cdss_constants, only: default_polar_n, default_polar_k, &
                            default_polar_list, default_interleaver_seed
  implicit none
  private
  public :: waveform_config, channel_config, soft_metric_config, coding_config
  public :: tone_interferer
  public :: soft_demod_result, frame_recovery_result, tracking_state

  !> Persistent receiver synchronization and tracking state.
  type :: tracking_state
    real(dp) :: phase_est = 0.0_dp      !< Current carrier phase estimate in radians.
    real(dp) :: freq_est  = 0.0_dp      !< Current carrier frequency estimate in radians per bit.
    real(dp) :: timing_error = 0.0_dp   !< Accumulated timing error in samples.
    real(dp) :: timing_freq_est = 0.0_dp !< Timing frequency estimate (samples per bit).
  end type tracking_state

  !> Physical-layer waveform geometry.
  !!
  !! The default configuration targets a narrow audio-band mode with
  !! direct-sequence spreading over a slow chirp carrier. The transmitter and
  !! matched receiver must use identical values.
  type :: waveform_config
    real(dp) :: sample_rate       = 6000.0_dp     !< Internal sample rate in hertz.
    real(dp) :: carrier_hz        = 1500.0_dp     !< Audio-frequency carrier in hertz.
    real(dp) :: chip_rate         = 250.0_dp      !< PN chip rate in chips per second.
    integer  :: spreading_factor  = 64            !< Number of chips per coded bit.
    real(dp) :: chirp_slope_hz_s  = 50.0_dp       !< Linear chirp slope in hertz per second.
    real(dp) :: pulse_taper       = 0.20_dp       !< Raised-cosine edge fraction per chip.
    integer  :: pn_seed           = 17            !< Gold-code phase offset and seed selector.
    integer  :: pn_length_bits    = 6             !< LFSR stage count for PN generation.
    integer  :: preamble_bits     = 31            !< Known synchronization preamble length.
  contains
    procedure :: samples_per_chip
    procedure :: samples_per_bit
    procedure :: bandwidth_hz
  end type waveform_config

  !> Channel-emulator impairments and receiver-side channel priors.
  type :: channel_config
    real(dp) :: snr_db          = -27.0_dp        !< Signal-to-noise ratio in decibels.
    real(dp) :: cfo_hz          = 0.0_dp          !< Carrier-frequency offset in hertz.
    real(dp) :: drift_hz_s      = 0.0_dp          !< Linear frequency drift in hertz per second.
    real(dp) :: clock_ppm       = 0.0_dp          !< Sample-clock error in parts per million.
    integer  :: n_tones         = 0               !< Number of continuous-wave interferers.
    real(dp) :: prefix_pad_s    = 0.0_dp          !< Leading silence duration in seconds.
    real(dp) :: suffix_pad_s    = 0.0_dp          !< Trailing silence duration in seconds.
  end type channel_config

  !> Continuous-wave interferer specification.
  type :: tone_interferer
    real(dp) :: frequency_hz = 0.0_dp             !< Interferer frequency in hertz.
    real(dp) :: sir_db       = 0.0_dp             !< Signal-to-interference ratio in decibels.
    real(dp) :: phase_rad    = 0.0_dp             !< Initial interferer phase in radians.
    real(dp) :: start_s      = 0.0_dp             !< Interferer start time in seconds.
    logical  :: has_duration = .false.            !< True when the tone is time-gated.
    real(dp) :: duration_s   = 0.0_dp             !< Tone duration in seconds when gated.
  end type tone_interferer

  !> Soft-demodulation and interference-exposure tuning parameters.
  type :: soft_metric_config
    logical  :: has_noise_var       = .false.     !< True when noise_var is externally supplied.
    real(dp) :: noise_var           = 1.0_dp      !< Complex-sample noise variance.
    real(dp) :: exposure_weight     = 0.20_dp     !< Weight applied to per-chip exposure.
    real(dp) :: exposure_clip       = 8.0_dp      !< Upper bound for normalized exposure.
    integer  :: stft_window         = 256         !< STFT window length in samples.
    integer  :: stft_hop            = 64          !< STFT hop size in samples.
    real(dp) :: temperature         = 1.0_dp      !< Multiplicative LLR calibration factor.
    real(dp) :: llr_clip            = 30.0_dp     !< Absolute LLR saturation limit.
    integer  :: refinement_iters    = 2           !< Maximum feedback-refinement iterations.
    logical  :: gmi_calibration     = .true.      !< Enable GMI-based temperature calibration.
    real(dp) :: outband_guard_hz    = 30.0_dp     !< Guard band for out-of-band floor estimation.
  end type soft_metric_config

  !> Polar coding and interleaver parameters.
  type :: coding_config
    integer :: header_n  = default_polar_n        !< Header codeword length.
    integer :: header_k  = default_polar_k        !< Header information length.
    integer :: payload_n = default_polar_n        !< Payload codeword length.
    integer :: payload_k = default_polar_k        !< Payload information length.
    integer :: list_size = default_polar_list     !< Successive-cancellation list size.
    integer :: interleaver_seed = default_interleaver_seed !< Base pseudo-random interleaver seed.
  end type coding_config

  !> Intermediate results produced by soft demodulation.
  type :: soft_demod_result
    integer,  allocatable :: hard_bits(:)        !< Hard coded-bit decisions before channel decoding.
    real(dp), allocatable :: bit_llrs(:)         !< LLRs where positive values favor coded bit 0.
    real(dp), allocatable :: chip_J(:)           !< Per-chip interference exposure estimate.
    real(dp), allocatable :: chip_weight(:)      !< Reliability weight assigned to each chip.
    real(dp)              :: noise_var = 1.0_dp  !< Noise variance used for LLR scaling.
    real(dp)              :: gmi       = 0.0_dp  !< Generalized mutual information estimate.
    real(dp)              :: temperature_used = 1.0_dp !< Applied LLR temperature factor.
  end type soft_demod_result

  !> Complete frame-recovery result and diagnostic metadata.
  type :: frame_recovery_result
    integer,  allocatable :: payload_bits(:)      !< Final decoded payload bits.
    integer,  allocatable :: polar_payload_bits(:) !< Payload bits selected by the SCL decoder.
    real(dp), allocatable :: payload_llrs(:)      !< Payload-channel LLRs passed to the decoder.
    type(soft_demod_result) :: soft_detection     !< Soft-detector diagnostics.
    integer :: acquisition_sample = 0             !< Estimated zero-based frame-start sample.
    real(dp) :: acquisition_score = 0.0_dp        !< Preamble correlation score.
    integer :: payload_bit_count       = 0        !< Recovered payload length in bits.
    integer :: payload_crc_expected    = 0        !< CRC read from the decoded frame header.
    integer :: payload_crc_observed    = 0        !< CRC computed from recovered payload bits.
    integer :: refinement_iters_used   = 0        !< Number of feedback iterations performed.
    logical :: crc_ok                  = .false.  !< True when expected and observed CRCs match.
  end type frame_recovery_result

contains

  !> Compute the integer number of samples per PN chip.
  !!
  !! @param[in] self Waveform configuration object.
  !! @return Samples per chip, rounded to the nearest integer and clamped to one.
  pure integer function samples_per_chip(self)
    class(waveform_config), intent(in) :: self
    samples_per_chip = max(1, nint(self%sample_rate / self%chip_rate))
  end function samples_per_chip

  !> Compute the number of samples occupied by one coded bit.
  !!
  !! @param[in] self Waveform configuration object.
  !! @return Samples per coded bit after spreading.
  pure integer function samples_per_bit(self)
    class(waveform_config), intent(in) :: self
    samples_per_bit = self%samples_per_chip() * self%spreading_factor
  end function samples_per_bit

  !> Estimate the occupied signal bandwidth.
  !!
  !! @param[in] self Waveform configuration object.
  !! @return Approximate chip-rate-limited bandwidth in hertz.
  pure real(dp) function bandwidth_hz(self)
    class(waveform_config), intent(in) :: self
    bandwidth_hz = self%chip_rate
  end function bandwidth_hz
end module cdss_types
