!> @file cdss_fft.f90
!! @brief FFTW3 bindings and STFT helper routines.
!!
!! @details
!! This module isolates FFTW plan management from receiver algorithms. Forward
!! and inverse complex transforms use cached per-thread plans, while fft_stft
!! provides the power spectrogram used by chip-level interference estimation.
module cdss_fft
  use cdss_kinds, only: dp
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
  private
  public :: fft_forward, fft_inverse, fft_stft

  !> Cached FFTW plan and aligned work buffers for one transform length.
  type :: fft_plan_cache
    integer       :: n     = 0                  !< Transform length associated with the plan.
    type(c_ptr)   :: plan  = c_null_ptr         !< Opaque FFTW plan handle.
    complex(c_double_complex), pointer :: in(:)  => null() !< Fortran view of input buffer.
    complex(c_double_complex), pointer :: out(:) => null() !< Fortran view of output buffer.
    type(c_ptr)   :: in_raw  = c_null_ptr       !< FFTW-owned input allocation.
    type(c_ptr)   :: out_raw = c_null_ptr       !< FFTW-owned output allocation.
  end type fft_plan_cache

  type(fft_plan_cache), save :: g_fwd, g_inv
  !$omp threadprivate(g_fwd, g_inv)

contains

  !> Ensure that a cached FFTW plan exists for the requested transform.
  !!
  !! @param[inout] cache Thread-local plan cache to update.
  !! @param[in]    n     Transform length.
  !! @param[in]    sign  FFTW_FORWARD or FFTW_BACKWARD.
  subroutine ensure_plan(cache, n, sign)
    type(fft_plan_cache), intent(inout) :: cache
    integer, intent(in) :: n, sign
    if (cache%n == n .and. c_associated(cache%plan)) return
    call release_plan(cache)
    cache%n = n
    cache%in_raw  = fftw_alloc_complex(int(n, c_size_t))
    cache%out_raw = fftw_alloc_complex(int(n, c_size_t))
    call c_f_pointer(cache%in_raw,  cache%in,  [n])
    call c_f_pointer(cache%out_raw, cache%out, [n])
    cache%plan = fftw_plan_dft_1d(int(n, c_int), cache%in, cache%out, &
                                  sign, FFTW_ESTIMATE)
  end subroutine ensure_plan

  !> Release an FFTW plan and its aligned buffers.
  !!
  !! @param[inout] cache Cache object to reset.
  subroutine release_plan(cache)
    type(fft_plan_cache), intent(inout) :: cache
    if (c_associated(cache%plan))    call fftw_destroy_plan(cache%plan)
    if (c_associated(cache%in_raw))  call fftw_free(cache%in_raw)
    if (c_associated(cache%out_raw)) call fftw_free(cache%out_raw)
    cache%plan    = c_null_ptr
    cache%in_raw  = c_null_ptr
    cache%out_raw = c_null_ptr
    cache%n       = 0
    nullify(cache%in)
    nullify(cache%out)
  end subroutine release_plan

  !> Compute an unnormalized complex forward FFT.
  !!
  !! @param[in]  input  Complex time-domain input.
  !! @param[out] output Allocated frequency-domain output.
  subroutine fft_forward(input, output)
    complex(dp), intent(in)  :: input(:)
    complex(dp), allocatable, intent(out) :: output(:)
    integer :: n
    n = size(input)
    call ensure_plan(g_fwd, n, FFTW_FORWARD)
    g_fwd%in(1:n) = input
    call fftw_execute_dft(g_fwd%plan, g_fwd%in, g_fwd%out)
    allocate(output(n))
    output = g_fwd%out(1:n)
  end subroutine fft_forward

  !> Compute a normalized complex inverse FFT.
  !!
  !! @param[in]  input  Complex frequency-domain input.
  !! @param[out] output Allocated time-domain output scaled by 1/N.
  subroutine fft_inverse(input, output)
    complex(dp), intent(in)  :: input(:)
    complex(dp), allocatable, intent(out) :: output(:)
    integer :: n
    n = size(input)
    call ensure_plan(g_inv, n, FFTW_BACKWARD)
    g_inv%in(1:n) = input
    call fftw_execute_dft(g_inv%plan, g_inv%in, g_inv%out)
    allocate(output(n))
    output = g_inv%out(1:n) / real(n, dp)
  end subroutine fft_inverse

  !> Compute a Hann-windowed short-time power spectrogram.
  !!
  !! @param[in]  signal      Complex input signal.
  !! @param[in]  nfft        Window length and FFT length.
  !! @param[in]  hop         Hop size in samples.
  !! @param[out] spectrogram Allocated matrix of power values, indexed as
  !!                         (frame, frequency bin).
  !! @param[out] n_frames    Number of STFT frames returned.
  subroutine fft_stft(signal, nfft, hop, spectrogram, n_frames)
    complex(dp), intent(in) :: signal(:)
    integer, intent(in)     :: nfft, hop
    real(dp), allocatable, intent(out) :: spectrogram(:, :)
    integer, intent(out)    :: n_frames
    complex(dp), allocatable :: frame_in(:), frame_out(:)
    real(dp), allocatable :: window(:)
    integer :: i, k, idx

    n_frames = max(0, (size(signal) - nfft) / max(hop, 1) + 1)
    if (n_frames <= 0) then
      allocate(spectrogram(0, 0))
      return
    end if

    allocate(spectrogram(n_frames, nfft))
    allocate(window(nfft))
    do k = 1, nfft
      window(k) = 0.5_dp - 0.5_dp * &
                  cos(2.0_dp * acos(-1.0_dp) * real(k - 1, dp) / real(nfft - 1, dp))
    end do

    do i = 1, n_frames
      allocate(frame_in(nfft))
      idx = (i - 1) * hop
      do k = 1, nfft
        frame_in(k) = signal(idx + k) * window(k)
      end do
      call fft_forward(frame_in, frame_out)
      spectrogram(i, :) = abs(frame_out)**2 / real(nfft, dp)
      deallocate(frame_in, frame_out)
    end do
    deallocate(window)
  end subroutine fft_stft
end module cdss_fft
