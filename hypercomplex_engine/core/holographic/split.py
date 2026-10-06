from ..validation import Validation
from .standard import StandardHolographic


class SplitHolographic:
    """
    Split Cayley-Dickson holographic O(n) multiplier.

    Architecture:
        Standard parent A_{dim-1} + one split doubling at the top.

    Blocks A, B, C:
        delegate to StandardHolographic.

    Block D:
        handled by split rules.
    """

    def __init__(self):
        self._standard = StandardHolographic()

    def multiply(self, t1: tuple, t2: tuple, dim: int) -> tuple:
        """
        Multiply two basis elements in split A_dim.

        t1 = (sign1, index1)
        t2 = (sign2, index2)

        Returns:
            (final_sign, index1 XOR index2)
        """
        Validation.basis_tuple(t1, allow_zero=False, allow_eps=False)
        Validation.basis_tuple(t2, allow_zero=False, allow_eps=False)

        s1, i = int(t1[0]), int(t1[1])
        s2, j = int(t2[0]), int(t2[1])

        dim = Validation.dimension(dim)
        Validation.index_in_range(i, dim)
        Validation.index_in_range(j, dim)

        sign, idx = self.multiply_indices(i, j, dim)

        return (s1 * s2 * sign, idx)

    def multiply_indices(self, i: int, j: int, dim: int) -> tuple:
        """
        Core split multiplication for non-negative indices.

        Returns:
            (sign, i XOR j)
        """
        i = int(i)
        j = int(j)
        dim = Validation.dimension(dim)

        if i < 0 or j < 0:
            raise ValueError("indices must be >= 0")

        Validation.index_in_range(i, dim)
        Validation.index_in_range(j, dim)

        if dim == 0:
            return (1, 0)

        half = 1 << (dim - 1)

        # --------------------------------------------------------------
        # Top-level Block D:
        # both indices are in the upper half.
        # --------------------------------------------------------------
        if i >= half and j >= half:
            i_loc = i - half
            j_loc = j - half

            # Split Block D diagonal:
            # (e_i ℓ)^2 = +e_0 for all i.
            if i_loc == j_loc:
                return (1, 0)

            # Split Block D first row:
            # i_loc == 0, j_loc > 0 -> sign -1.
            if i_loc == 0:
                return (-1, i ^ j)

            # Split Block D first column:
            # j_loc == 0, i_loc > 0 -> sign +1.
            if j_loc == 0:
                return (1, i ^ j)

            # Split Block D interior:
            # sign = sigma_a(i_loc, j_loc).
            sign, _ = self._standard.multiply_indices(i_loc, j_loc)
            return (sign, i ^ j)

        # --------------------------------------------------------------
        # Blocks A, B, C:
        # identical to standard.
        # --------------------------------------------------------------
        return self._standard.multiply_indices(i, j)