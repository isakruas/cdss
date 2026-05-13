!> @file cdss_polar.f90
!! @brief Polar encoder and successive-cancellation list decoder.
!!
!! @details
!! The implementation supports short binary Polar codes with power-of-two
!! block length. The information set is selected by a classical Bhattacharyya
!! construction for a binary erasure channel design point, the encoder applies
!! the standard in-place Polar transform, and the decoder uses an LLR-domain
!! successive-cancellation list (SCL) search.
!!
!! LLR sign convention: positive LLR means bit 0 is more likely; negative LLR
!! means bit 1 is more likely.
module cdss_polar
  use cdss_kinds, only: dp, i32
  use cdss_constants, only: eps
  implicit none
  private
  public :: polar_code
  public :: polar_encode_stream, polar_decode_stream
  public :: polar_information_set

  !> Polar code descriptor shared by encoder and decoder.
  type :: polar_code
    integer :: n = 128          !< Codeword length; must be a power of two.
    integer :: k = 64           !< Number of information bits per block.
    integer :: list_size = 8    !< Maximum number of active SCL paths.
  end type polar_code

contains

  !> Construct the Polar information set using Bhattacharyya reliability values.
  !!
  !! @param[in]  n         Codeword length. Must be a power of two.
  !! @param[in]  k         Number of information-bit positions to select.
  !! @param[out] info_mask Allocated mask: 1 for information positions, 0 for
  !!                       frozen positions.
  subroutine polar_information_set(n, k, info_mask)
    integer, intent(in) :: n, k
    integer, allocatable, intent(out) :: info_mask(:)   ! 1 if info, 0 if frozen
    real(dp), allocatable :: z(:), znext(:)
    integer, allocatable :: order(:)
    integer :: i, levels, lvl, half, p
    real(dp) :: eps0

    levels = nint(log(real(n, dp)) / log(2.0_dp))
    if (ishft(1, levels) /= n) error stop "polar_information_set: N must be power of 2"

    allocate(z(n), znext(n))
    eps0 = 0.32_dp
    z = 0.0_dp
    z(1) = eps0

    do lvl = 1, levels
      half = ishft(1, lvl - 1)
      do i = 1, half
        znext(2 * i - 1) = 2.0_dp * z(i) - z(i)**2
        znext(2 * i)     = z(i)**2
      end do
      z(1:2 * half) = znext(1:2 * half)
    end do

    ! Sort indices by Z parameter ascending; the best channels are first.
    allocate(order(n))
    do i = 1, n
      order(i) = i
    end do
    call sort_indices_by_value(order, z)

    allocate(info_mask(n))
    info_mask = 0
    do i = 1, k
      info_mask(order(i)) = 1
    end do
    deallocate(z, znext, order)
  end subroutine polar_information_set

  !> Encode one Polar block.
  !!
  !! @param[in]  code  Polar code descriptor.
  !! @param[in]  info  Information bits. The first code%k values are consumed.
  !! @param[out] coded Preallocated output codeword with length code%n.
  subroutine polar_encode_block(code, info, coded)
    type(polar_code), intent(in) :: code
    integer, intent(in)  :: info(:)
    integer, intent(out) :: coded(:)
    integer, allocatable :: u(:), info_mask(:)
    integer :: i, info_count, half, step, j

    call polar_information_set(code%n, code%k, info_mask)
    allocate(u(code%n))
    u = 0
    info_count = 0
    do i = 1, code%n
      if (info_mask(i) == 1) then
        info_count = info_count + 1
        u(i) = iand(info(info_count), 1)
      end if
    end do

    ! In-place polar transform on u.
    half = 1
    do while (half < code%n)
      step = half * 2
      do i = 1, code%n, step
        do j = 0, half - 1
          u(i + j) = ieor(u(i + j), u(i + j + half))
        end do
      end do
      half = half * 2
    end do

    coded = u
    deallocate(u, info_mask)
  end subroutine polar_encode_block

  !> Encode an arbitrary-length bit stream using repeated Polar blocks.
  !!
  !! The final block is zero-padded when the input length is not a multiple of
  !! code%k.
  !!
  !! @param[in]  code  Polar code descriptor.
  !! @param[in]  info  Input information bit stream.
  !! @param[out] coded Allocated concatenated codewords.
  subroutine polar_encode_stream(code, info, coded)
    type(polar_code), intent(in) :: code
    integer, intent(in) :: info(:)
    integer, allocatable, intent(out) :: coded(:)
    integer :: blocks, b
    integer, allocatable :: info_block(:), coded_block(:)
    blocks = (size(info) + code%k - 1) / code%k
    allocate(coded(blocks * code%n))
    allocate(info_block(code%k), coded_block(code%n))
    coded = 0
    do b = 1, blocks
      info_block = 0
      info_block(1:min(code%k, size(info) - (b - 1) * code%k)) = &
           info((b - 1) * code%k + 1:min(b * code%k, size(info)))
      call polar_encode_block(code, info_block, coded_block)
      coded((b - 1) * code%n + 1:b * code%n) = coded_block
    end do
    deallocate(info_block, coded_block)
  end subroutine polar_encode_stream

  !> Decode one Polar block with an LLR-domain SCL decoder.
  !!
  !! @param[in]  code     Polar code descriptor.
  !! @param[in]  llrs     Channel LLRs with positive values favoring bit 0.
  !! @param[out] info_out Preallocated decoded information bits of length code%k.
  subroutine polar_decode_block(code, llrs, info_out)
    type(polar_code), intent(in)  :: code
    real(dp), intent(in)  :: llrs(:)
    integer,  intent(out) :: info_out(:)
    integer, allocatable :: info_mask(:)
    integer :: levels, n, l_max
    real(dp), allocatable :: llr_tree(:, :, :)         ! (depth, idx, path)
    integer,  allocatable :: u_tree  (:, :, :)         ! decoded bits per path
    real(dp), allocatable :: path_metric(:)
    logical,  allocatable :: path_active(:)
    integer,  allocatable :: parent_path(:)
    real(dp), allocatable :: scratch_metric(:)
    integer :: phi, l, l2, info_count, i, depth, idx, h, p, best
    real(dp) :: ml_0, ml_1, metric_split

    n = code%n
    levels = nint(log(real(n, dp)) / log(2.0_dp))
    l_max = code%list_size

    call polar_information_set(code%n, code%k, info_mask)

    allocate(llr_tree(0:levels, n, l_max))
    allocate(u_tree(0:levels, n, l_max))
    allocate(path_metric(l_max), path_active(l_max))
    path_metric = 0.0_dp
    path_active = .false.
    path_active(1) = .true.
    do l = 1, l_max
      llr_tree(0, 1:n, l) = llrs(1:n)
      u_tree(:, :, l) = 0
    end do

    ! Walk leaves left-to-right.
    do phi = 0, n - 1
      ! Compute LLRs down to leaf for every active path.
      do l = 1, l_max
        if (.not. path_active(l)) cycle
        call recurse_llr(phi, levels, llr_tree(:, :, l), u_tree(:, :, l))
      end do

      if (info_mask(phi + 1) == 0) then
        ! Frozen: decode 0, update metric.
        do l = 1, l_max
          if (.not. path_active(l)) cycle
          path_metric(l) = path_metric(l) + log1p_exp(-llr_tree(levels, phi + 1, l))
          u_tree(levels, phi + 1, l) = 0
        end do
      else
        ! Info bit: split each active path into two candidates, keep best L_max.
        call split_paths(l_max, levels, phi, info_mask, &
             llr_tree, u_tree, path_metric, path_active)
      end if

      ! Propagate decided bit up the partial-sum tree.
      do l = 1, l_max
        if (.not. path_active(l)) cycle
        call recurse_partial_sum(phi, levels, u_tree(:, :, l))
      end do
    end do

    ! Pick path with minimum metric; export info-bit positions.
    best = 1
    do l = 2, l_max
      if (.not. path_active(l)) cycle
      if (path_metric(l) < path_metric(best) .or. .not. path_active(best)) best = l
    end do

    info_count = 0
    do i = 1, n
      if (info_mask(i) == 1) then
        info_count = info_count + 1
        info_out(info_count) = u_tree(levels, i, best)
      end if
    end do

    deallocate(llr_tree, u_tree, path_metric, path_active, info_mask)
  end subroutine polar_decode_block

  !> Decode a concatenated stream of complete Polar codewords.
  !!
  !! @param[in]  code      Polar code descriptor.
  !! @param[in]  llrs      Concatenated channel LLRs. Extra samples beyond a
  !!                       whole number of codewords are ignored by integer division.
  !! @param[out] info_bits Allocated decoded information stream.
  subroutine polar_decode_stream(code, llrs, info_bits)
    type(polar_code), intent(in) :: code
    real(dp), intent(in) :: llrs(:)
    integer, allocatable, intent(out) :: info_bits(:)
    integer :: blocks, b
    real(dp), allocatable :: chunk(:)
    integer,  allocatable :: info_chunk(:)
    blocks = size(llrs) / code%n
    allocate(info_bits(blocks * code%k))
    allocate(chunk(code%n), info_chunk(code%k))
    do b = 1, blocks
      chunk = llrs((b - 1) * code%n + 1:b * code%n)
      call polar_decode_block(code, chunk, info_chunk)
      info_bits((b - 1) * code%k + 1:b * code%k) = info_chunk
    end do
    deallocate(chunk, info_chunk)
  end subroutine polar_decode_stream

  !> Recursively compute LLR messages down to the current leaf.
  !!
  !! @param[in]    phi    Zero-based bit index currently being decoded.
  !! @param[in]    levels log2(codeword length).
  !! @param[inout] llr    LLR tree for one active SCL path.
  !! @param[inout] u      Partial-sum tree for one active SCL path.
  recursive subroutine recurse_llr(phi, levels, llr, u)
    integer, intent(in) :: phi, levels
    real(dp), intent(inout) :: llr(0:, :)
    integer,  intent(inout) :: u(0:, :)
    integer :: depth, span, half, idx, j, bit
    real(dp) :: la, lb, sign_a, sign_b

    ! Walk down: for each level from 0 to levels-1, depending on whether phi
    ! sits in the upper or lower half of the corresponding partition, compute
    ! either the f-function (XOR LLR combine) or the g-function (sum LLR
    ! combine using the already-decoded partial sum).
    do depth = 0, levels - 1
      span = ishft(1, levels - depth)
      half = span / 2
      idx  = (phi / span) * span + 1
      if (mod(phi / half, 2) == 0) then
        ! f-step on (idx ... idx+span-1) into the first half at depth+1.
        do j = 0, half - 1
          la = llr(depth, idx + j)
          lb = llr(depth, idx + j + half)
          sign_a = sign(1.0_dp, la)
          sign_b = sign(1.0_dp, lb)
          llr(depth + 1, idx + j) = sign_a * sign_b * min(abs(la), abs(lb))
        end do
      else
        ! g-step using the partial sum stored in u at depth+1.
        do j = 0, half - 1
          la = llr(depth, idx + j)
          lb = llr(depth, idx + j + half)
          bit = u(depth + 1, idx + j)
          if (bit == 0) then
            llr(depth + 1, idx + j + half) = lb + la
          else
            llr(depth + 1, idx + j + half) = lb - la
          end if
        end do
      end if
    end do
  end subroutine recurse_llr

  !> Propagate a newly decided leaf bit upward through the partial-sum tree.
  !!
  !! @param[in]    phi    Zero-based bit index just decided.
  !! @param[in]    levels log2(codeword length).
  !! @param[inout] u      Partial-sum tree for one active SCL path.
  recursive subroutine recurse_partial_sum(phi, levels, u)
    integer, intent(in) :: phi, levels
    integer, intent(inout) :: u(0:, :)
    integer :: depth, span, half, idx, j

    ! After deciding leaf at index phi+1 (level=levels), propagate XOR up.
    do depth = levels - 1, 0, -1
      span = ishft(1, levels - depth)
      half = span / 2
      idx  = (phi / span) * span + 1
      if (mod(phi / half, 2) == 1) then
        ! The right child has been decided, so the parent can be finalized.
        do j = 0, half - 1
          u(depth, idx + j)        = ieor(u(depth + 1, idx + j), &
                                          u(depth + 1, idx + j + half))
          u(depth, idx + j + half) = u(depth + 1, idx + j + half)
        end do
      else
        return                                     ! still need right child
      end if
    end do
  end subroutine recurse_partial_sum

  !> Split active SCL paths on an information bit and retain the best candidates.
  !!
  !! @param[in]    l_max       Maximum number of list paths.
  !! @param[in]    levels      log2(codeword length).
  !! @param[in]    phi         Zero-based information-bit index.
  !! @param[in]    info_mask   Polar information-position mask.
  !! @param[inout] llr_tree    LLR state for all paths.
  !! @param[inout] u_tree      Partial-sum state for all paths.
  !! @param[inout] path_metric Accumulated path metrics.
  !! @param[inout] path_active Active-path flags.
  subroutine split_paths(l_max, levels, phi, info_mask, llr_tree, u_tree, &
                          path_metric, path_active)
    integer, intent(in) :: l_max, levels, phi
    integer, intent(in) :: info_mask(:)
    real(dp), intent(inout) :: llr_tree(0:, :, :)
    integer,  intent(inout) :: u_tree(0:, :, :)
    real(dp), intent(inout) :: path_metric(:)
    logical,  intent(inout) :: path_active(:)

    real(dp), allocatable :: cand_metric(:)
    integer,  allocatable :: cand_path(:), cand_bit(:), order(:)
    logical,  allocatable :: parent_was(:)
    real(dp), allocatable :: new_llr(:, :, :)
    integer,  allocatable :: new_u(:, :, :)
    real(dp), allocatable :: new_metric(:)
    logical,  allocatable :: new_active(:)

    integer :: n_active, k, l, candidate_count, n_keep, i, p, parent, bit, n, levels_local
    real(dp) :: m0, m1, llr_v

    n_active = count(path_active)
    candidate_count = 2 * n_active
    allocate(cand_metric(candidate_count), cand_path(candidate_count), cand_bit(candidate_count))

    k = 0
    do l = 1, l_max
      if (.not. path_active(l)) cycle
      llr_v = llr_tree(levels, phi + 1, l)
      m0 = path_metric(l) + log1p_exp(-llr_v)
      m1 = path_metric(l) + log1p_exp( llr_v)
      k = k + 1; cand_metric(k) = m0; cand_path(k) = l; cand_bit(k) = 0
      k = k + 1; cand_metric(k) = m1; cand_path(k) = l; cand_bit(k) = 1
    end do

    allocate(order(candidate_count))
    do i = 1, candidate_count
      order(i) = i
    end do
    call sort_indices_by_value(order, cand_metric)

    n_keep = min(l_max, candidate_count)

    ! Materialize the surviving candidates into a NEW slot bank, then copy back.
    n = size(llr_tree, 2)
    levels_local = levels
    allocate(new_llr(0:levels_local, n, l_max))
    allocate(new_u(0:levels_local, n, l_max))
    allocate(new_metric(l_max), new_active(l_max))
    new_metric = 0.0_dp
    new_active = .false.
    new_llr    = 0.0_dp
    new_u      = 0

    do i = 1, n_keep
      parent = cand_path(order(i))
      bit    = cand_bit(order(i))
      new_llr(:, :, i)        = llr_tree(:, :, parent)
      new_u(:, :, i)          = u_tree(:, :, parent)
      new_u(levels_local, phi + 1, i) = bit
      new_metric(i)           = cand_metric(order(i))
      new_active(i)           = .true.
    end do

    llr_tree    = new_llr
    u_tree      = new_u
    path_metric = new_metric
    path_active = new_active

    deallocate(cand_metric, cand_path, cand_bit, order)
    deallocate(new_llr, new_u, new_metric, new_active)
  end subroutine split_paths

  !> Evaluate log(1 + exp(x)) without overflow.
  !!
  !! @param[in] x Real-valued argument.
  !! @return Stable softplus value.
  pure real(dp) function log1p_exp(x) result(y)
    real(dp), intent(in) :: x
    if (x > 30.0_dp) then
      y = x
    else if (x < -30.0_dp) then
      y = 0.0_dp
    else
      y = log(1.0_dp + exp(x))
    end if
  end function log1p_exp

  !> Sort an index permutation according to associated real values.
  !!
  !! @param[inout] order  Index permutation to sort in ascending value order.
  !! @param[in]    values Values addressed by the indices in order.
  subroutine sort_indices_by_value(order, values)
    integer,  intent(inout) :: order(:)
    real(dp), intent(in) :: values(:)
    integer :: i, j, key
    do i = 2, size(order)
      key = order(i)
      j = i - 1
      do while (j >= 1)
        if (values(order(j)) <= values(key)) exit
        order(j + 1) = order(j)
        j = j - 1
      end do
      order(j + 1) = key
    end do
  end subroutine sort_indices_by_value
end module cdss_polar
