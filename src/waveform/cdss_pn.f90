!> @file cdss_pn.f90
!! @brief Pseudo-noise sequence generation for DSSS spreading.
!!
!! @details
!! The modem uses maximum-length LFSR sequences and Gold-code families to obtain
!! deterministic chip signs with bounded cross-correlation. Tap positions are
!! documented as one-based register positions, where position 1 is the most
!! significant (oldest) bit in the internal state representation.
module cdss_pn
  use cdss_kinds, only: dp, i32
  implicit none
  private
  public :: pn_msequence, pn_gold, pn_to_npolar

  ! Preferred-pair m-sequence feedback taps for supported LFSR sizes.
  integer, parameter :: TAPS_5_A(2)  = [5, 3]        !< First n=5 primitive polynomial taps.
  integer, parameter :: TAPS_5_B(4)  = [5, 4, 3, 2] !< Second n=5 primitive polynomial taps.
  integer, parameter :: TAPS_6_A(2)  = [6, 5]        !< First n=6 primitive polynomial taps.
  integer, parameter :: TAPS_6_B(4)  = [6, 5, 3, 2] !< Second n=6 primitive polynomial taps.
  integer, parameter :: TAPS_7_A(2)  = [7, 6]        !< First n=7 primitive polynomial taps.
  integer, parameter :: TAPS_7_B(4)  = [7, 6, 5, 4] !< Second n=7 primitive polynomial taps.
  integer, parameter :: TAPS_8_A(4)  = [8, 6, 5, 4] !< First n=8 primitive polynomial taps.
  integer, parameter :: TAPS_8_B(4)  = [8, 6, 5, 3] !< Second n=8 primitive polynomial taps.
  integer, parameter :: TAPS_9_A(2)  = [9, 5]        !< First n=9 primitive polynomial taps.
  integer, parameter :: TAPS_9_B(4)  = [9, 6, 4, 3] !< Second n=9 primitive polynomial taps.
  integer, parameter :: TAPS_10_A(2) = [10, 7]       !< First n=10 primitive polynomial taps.
  integer, parameter :: TAPS_10_B(4) = [10, 8, 3, 2] !< Second n=10 primitive polynomial taps.

contains

  !> Generate a binary maximum-length sequence with a Fibonacci LFSR.
  !!
  !! @param[in]  stages Number of LFSR stages. Must be at least two.
  !! @param[in]  taps   Feedback tap positions, one-based from the most
  !!                    significant bit.
  !! @param[in]  seed   Initial state. A zero seed is coerced to a nonzero state.
  !! @param[in]  length Number of output chips.
  !! @param[out] out    Allocated sequence containing values 0 or 1.
  subroutine pn_msequence(stages, taps, seed, length, out)
    integer, intent(in)  :: stages
    integer, intent(in)  :: taps(:)
    integer, intent(in)  :: seed
    integer, intent(in)  :: length
    integer, allocatable, intent(out) :: out(:)
    integer :: state, i, fb, j
    if (stages < 2) error stop "pn_msequence: need at least 2 stages"
    if (length <= 0) error stop "pn_msequence: length must be positive"
    state = iand(seed, ishft(1, stages) - 1)
    if (state == 0) state = 1                    ! avoid all-zero lock-up
    ! Force the MSB bit on so the feedback computation in the first iteration
    ! cannot collapse to zero; this guarantees entry into the full LFSR cycle.
    state = ior(state, ishft(1, stages - 1))
    allocate(out(length))
    do i = 1, length
      ! Output is the LSB of the state.
      out(i) = iand(state, 1)
      ! XOR the tapped bits to compute new feedback.
      fb = 0
      do j = 1, size(taps)
        fb = ieor(fb, iand(ishft(state, -(stages - taps(j))), 1))
      end do
      ! Shift right; insert feedback into the MSB of the n-stage register.
      state = ior(ishft(state, -1), ishft(fb, stages - 1))
    end do
  end subroutine pn_msequence

  !> Generate a Gold-code sequence from two preferred-pair m-sequences.
  !!
  !! @param[in]  stages      LFSR stage count. Supported values are 5 through 10.
  !! @param[in]  length      Number of output chips to generate.
  !! @param[in]  seed_offset Cyclic offset applied to the second m-sequence.
  !! @param[out] out         Allocated Gold-code chips containing values 0 or 1.
  subroutine pn_gold(stages, length, seed_offset, out)
    integer, intent(in)  :: stages, length, seed_offset
    integer, allocatable, intent(out) :: out(:)
    integer, allocatable :: a(:), b(:), b_shifted(:)
    integer :: i, period, off

    period = ishft(1, stages) - 1
    off = modulo(seed_offset, period)

    select case (stages)
    case (5)
      call pn_msequence(stages, TAPS_5_A, 1, length, a)
      call pn_msequence(stages, TAPS_5_B, 1, length + period, b)
    case (6)
      call pn_msequence(stages, TAPS_6_A, 1, length, a)
      call pn_msequence(stages, TAPS_6_B, 1, length + period, b)
    case (7)
      call pn_msequence(stages, TAPS_7_A, 1, length, a)
      call pn_msequence(stages, TAPS_7_B, 1, length + period, b)
    case (8)
      call pn_msequence(stages, TAPS_8_A, 1, length, a)
      call pn_msequence(stages, TAPS_8_B, 1, length + period, b)
    case (9)
      call pn_msequence(stages, TAPS_9_A, 1, length, a)
      call pn_msequence(stages, TAPS_9_B, 1, length + period, b)
    case (10)
      call pn_msequence(stages, TAPS_10_A, 1, length, a)
      call pn_msequence(stages, TAPS_10_B, 1, length + period, b)
    case default
      error stop "pn_gold: unsupported stage count (5..10 supported)"
    end select

    allocate(b_shifted(length))
    do i = 1, length
      b_shifted(i) = b(off + i)
    end do
    allocate(out(length))
    do i = 1, length
      out(i) = ieor(a(i), b_shifted(i))
    end do
    deallocate(a, b, b_shifted)
  end subroutine pn_gold

  !> Map binary PN chips to the BPSK sign convention used by the waveform.
  !!
  !! @param[in]  bits  Input PN chips represented as 0 or 1.
  !! @param[out] signs Output signs where 1 maps to +1 and 0 maps to -1.
  pure subroutine pn_to_npolar(bits, signs)
    integer, intent(in) :: bits(:)
    real(dp), intent(out) :: signs(:)
    integer :: i
    do i = 1, size(bits)
      if (iand(bits(i), 1) == 1) then
        signs(i) = 1.0_dp
      else
        signs(i) = -1.0_dp
      end if
    end do
  end subroutine pn_to_npolar
end module cdss_pn
