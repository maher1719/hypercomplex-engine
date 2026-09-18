from functools import lru_cache
import math


class HolographicMultiplier:
    """
    Holographic O(n) Cayley–Dickson multiplier with tuple-based interface.

    Accepts and returns tuples:
      standard / split  →  (sign, index)
      dual              →  (sign, index, eps_flag)

    Complexity: O(n) where n = ⌈log₂(max(i, j) + 1)⌉
    """

    def __init__(self):
        self._clear_cache()

    def _clear_cache(self):
        self._get_cell_cached.cache_clear()

    # ============================================================
    # VALIDATION  (matches FastMultiplier style)
    # ============================================================

    @staticmethod
    def validate_tuple(data):
        if not isinstance(data, tuple):
            raise TypeError(f"Expected a tuple, but got {type(data).__name__}")
        if len(data) < 2 or len(data) > 3:
            raise ValueError(f"Tuple must have 2 or 3 elements, but has {len(data)}")
        if not isinstance(data[0], int) or isinstance(data[0], bool):
            raise TypeError(f"sign must be int, got {type(data[0]).__name__}")
        if not isinstance(data[1], int) or isinstance(data[1], bool):
            raise TypeError(f"index must be int, got {type(data[1]).__name__}")
        if data[1] < 0:
            raise ValueError(f"index must be >= 0, got {data[1]}")
        if data[0] not in (-1, 1):
            raise ValueError(f"sign must be -1 or 1, got {data[0]}")
        return True

    # ============================================================
    # INTERNAL: STANDARD CD RECURSIVE DESCENT  (unchanged logic)
    # ============================================================

    @lru_cache(maxsize=None)
    def _get_cell_cached(self, n: int, i: int, j: int) -> tuple:
        if n == 0:
            return (1, 0)

        half = 1 << (n - 1)

        if i < half and j < half:
            return self._get_cell_cached(n - 1, i, j)

        if i < half and j >= half:
            quadrant = 'b'; i_loc, j_loc = i, j % half
        elif i >= half and j < half:
            quadrant = 'c'; i_loc, j_loc = i % half, j
        else:
            quadrant = 'd'; i_loc, j_loc = i % half, j % half

        ancestor_sign, ancestor_val = self._get_cell_cached(n - 1, i_loc, j_loc)
        final_val = (ancestor_val + half) if quadrant in ('b', 'c') else ancestor_val

        if self._is_structural(n, i, j, quadrant):
            final_sign = self._assign_sign_structural(n, i, j, quadrant)
        else:
            final_sign = -ancestor_sign

        return (final_sign, final_val)

    def _is_structural(self, n, i, j, quadrant):
        half = 1 << (n - 1)
        i_loc, j_loc = i % half, j % half

        if i == 0 or j == 0 or i == j:
            return True
        if quadrant == 'b':
            return j_loc == 0 or i_loc == j_loc
        if quadrant == 'c':
            return i_loc == 0 or i_loc == j_loc
        if quadrant == 'd':
            return i_loc == 0 or j_loc == 0 or i_loc == j_loc
        return False

    def _assign_sign_structural(self, n, i, j, quadrant):
        half = 1 << (n - 1)
        i_loc, j_loc = i % half, j % half

        if i == 0 or j == 0:
            return 1
        if i == j:
            return -1

        if quadrant == 'b':
            if j_loc == 0: return 1
            if i_loc == j_loc: return -1
        elif quadrant == 'c':
            if i_loc == 0:
                return 1 if j_loc == 0 else -1
            if i_loc == j_loc: return 1
        elif quadrant == 'd':
            if i_loc == 0: return 1
            if j_loc == 0: return -1
            if i_loc == j_loc: return -1

        raise ValueError("Invalid structural position.")

    @staticmethod
    def _dimension_level(max_index: int) -> int:
        if max_index == 0:
            return 0
        return math.floor(math.log2(max_index)) + 1

    def _multiply_sd_unsigned(self, i: int, j: int) -> tuple:
        """Standard CD sign + index for non-negative i, j."""
        if i == 0 and j == 0:
            return (1, 0)
        n = self._dimension_level(max(i, j))
        return self._get_cell_cached(n, i, j)

    # ============================================================
    # PUBLIC API — tuple in, tuple out
    # ============================================================

    def standard(self, t1: tuple, t2: tuple) -> tuple:
        """
        Standard Cayley–Dickson O(n) multiplication.
        (sign₁, idx₁) × (sign₂, idx₂) → (final_sign, idx₁ ⊕ idx₂)
        """
        self.validate_tuple(t1)
        self.validate_tuple(t2)

        s1, i = t1[0], t1[1]
        s2, j = t2[0], t2[1]
        outside_sign = s1 * s2

        sign, idx = self._multiply_sd_unsigned(i, j)
        return (outside_sign * sign, idx)

    def split(self, t1: tuple, t2: tuple, dim: int) -> tuple:
        """
        Split Cayley–Dickson O(n) multiplication for A_dim.
        Standard parent A_{dim-1} + one split doubling at top.
        """
        self.validate_tuple(t1)
        self.validate_tuple(t2)

        s1, i = t1[0], t1[1]
        s2, j = t2[0], t2[1]
        outside_sign = s1 * s2

        if dim < 0:
            raise ValueError("dim must be >= 0")

        size = 1 << dim
        if i >= size or j >= size:
            raise ValueError(f"indices must be in [0, {size - 1}] for dim={dim}")

        if dim == 0:
            return (outside_sign, 0)

        if i == 0 or j == 0:
            return (outside_sign, i ^ j)

        half = 1 << (dim - 1)

        # Global diagonal
        if i == j:
            base_sign = 1 if i >= half else -1
            return (outside_sign * base_sign, 0)

        # Top-level Block D
        if i >= half and j >= half:
            i_loc = i - half
            j_loc = j - half

            if i_loc == j_loc:
                return (outside_sign, 0)           # split diagonal: +1
            if i_loc == 0:
                return (-outside_sign, i ^ j)      # split first row: -1
            if j_loc == 0:
                return (outside_sign, i ^ j)       # split first col: +1

            # Block D interior = standard Block A interior
            sign, _ = self._multiply_sd_unsigned(i_loc, j_loc)
            return (outside_sign * sign, i ^ j)

        # Blocks A, B, C: identical to standard
        sign, idx = self._multiply_sd_unsigned(i, j)
        return (outside_sign * sign, idx)

    def dual(self, t1: tuple, t2: tuple, dim: int, split: bool = False) -> tuple:
        """
        Dual Cayley–Dickson O(n) multiplication for A_dim[ε]/(ε²).
        Total dimension 2^(dim+1); bit `dim` selects the ε-component.

        Returns: (final_sign, value_index, eps_flag)
        Convention: sign == 0 means the product is ZERO (ε-nilpotency).
        """
        self.validate_tuple(t1)
        self.validate_tuple(t2)

        if dim < 0:
            raise ValueError("dim must be >= 0")

        s1, i = t1[0], t1[1]
        s2, j = t2[0], t2[1]
        outside_sign = s1 * s2

        size = 1 << (dim + 1)
        if i >= size or j >= size:
            raise ValueError(f"indices must be in [0, {size - 1}] for dim={dim}")

        half = 1 << dim
        i_eps = i >= half
        j_eps = j >= half
        i_loc = i & (half - 1)
        j_loc = j & (half - 1)

        # Block d: ε-nilpotency  (εe_i)(εe_j) = 0
        if i_eps and j_eps:
            return (0, 0, 1)

        # Blocks a, b, c: delegate base multiplication
        if split:
            base_sign, idx = self.split((1, i_loc), (1, j_loc), dim)
        else:
            base_sign, idx = self.standard((1, i_loc), (1, j_loc))

        final_sign = outside_sign * base_sign
        eps_flag = 1 if (i_eps or j_eps) else 0

        return (final_sign, idx, eps_flag)