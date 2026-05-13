!> @file cdss_rng.f90
!! @brief Reproducible random-number utilities for simulation experiments.
!!
!! @details
!! The implementation wraps the Fortran intrinsic random_number interface and
!! stores a seed vector per rng_state. This gives deterministic Monte Carlo
!! trials when the same compiler/runtime and seed are used.
module cdss_rng
  use cdss_kinds, only: dp, i32
  use cdss_constants, only: two_pi
  implicit none
  private
  public :: rng_state, rng_seed, rng_uniform, rng_normal, rng_bits, rng_int_range

  !> Seedable pseudo-random generator state.
  type :: rng_state
    integer, allocatable :: seed_vec(:)       !< Runtime-specific seed vector.
    logical :: initialized = .false.          !< True after random_seed has been set.
  end type rng_state

contains

  !> Initialize the intrinsic generator with a deterministic integer seed.
  !!
  !! @param[inout] self Generator state to initialize.
  !! @param[in]    seed User-visible scalar seed.
  subroutine rng_seed(self, seed)
    type(rng_state), intent(inout) :: self
    integer, intent(in) :: seed
    integer :: n, i
    call random_seed(size=n)
    if (allocated(self%seed_vec)) deallocate(self%seed_vec)
    allocate(self%seed_vec(n))
    do i = 1, n
      self%seed_vec(i) = seed + 9973 * (i - 1) + 17
    end do
    call random_seed(put=self%seed_vec)
    self%initialized = .true.
  end subroutine rng_seed

  !> Draw a uniform random number in the interval [0, 1).
  !!
  !! @param[inout] self Generator state. It is initialized automatically if needed.
  !! @return Uniform variate in double precision.
  function rng_uniform(self) result(value)
    type(rng_state), intent(inout) :: self
    real(dp) :: value
    if (.not. self%initialized) call rng_seed(self, 0)
    call random_number(value)
  end function rng_uniform

  !> Draw a standard normal variate using the Box-Muller transform.
  !!
  !! @param[inout] self Generator state.
  !! @return Normally distributed sample with mean zero and variance one.
  function rng_normal(self) result(value)
    type(rng_state), intent(inout) :: self
    real(dp) :: value, u1, u2
    u1 = max(rng_uniform(self), 1.0e-300_dp)
    u2 = rng_uniform(self)
    value = sqrt(-2.0_dp * log(u1)) * cos(two_pi * u2)
  end function rng_normal

  !> Draw independent Bernoulli(0.5) bits.
  !!
  !! @param[inout] self Generator state.
  !! @param[in]    n    Number of bits to generate.
  !! @param[out]   bits Allocated output bit vector.
  subroutine rng_bits(self, n, bits)
    type(rng_state), intent(inout) :: self
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: bits(:)
    integer :: i
    allocate(bits(n))
    do i = 1, n
      bits(i) = merge(1, 0, rng_uniform(self) >= 0.5_dp)
    end do
  end subroutine rng_bits

  !> Draw an integer uniformly from an inclusive interval.
  !!
  !! @param[inout] self Generator state.
  !! @param[in]    lo   Lower inclusive bound.
  !! @param[in]    hi   Upper inclusive bound.
  !! @return Integer sample in [lo, hi].
  integer function rng_int_range(self, lo, hi)
    type(rng_state), intent(inout) :: self
    integer, intent(in) :: lo, hi
    rng_int_range = lo + int(rng_uniform(self) * real(hi - lo + 1, dp))
    rng_int_range = min(rng_int_range, hi)
  end function rng_int_range
end module cdss_rng
