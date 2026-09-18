from ..validation import Validation
from .standard import StandardHolographic
from .split import SplitHolographic


class DualHolographic:
    """
    Dual Cayley-Dickson holographic multiplier.

    Input index convention:
        The total dual dimension is 2^(dim + 1).
        Bit `dim` marks the epsilon component.

        lower half: 0 ... 2^dim - 1
        upper half: epsilon times lower half

    Output:
        (sign, local_index, eps_flag)

    Nilpotency:
        (epsilon e_i) * (epsilon e_j) = 0
        returned as (0, 0, 1)
    """

    def __init__(self, split: bool = False):
        self._split_mode = bool(split)
        self._standard = StandardHolographic()
        self._split = SplitHolographic() if self._split_mode else None

    def multiply(self, t1: tuple, t2: tuple, dim: int) -> tuple:
        """
        Multiply two dual basis elements.

        t1, t2 may be:
            (sign, global_index)
            (sign, global_index, eps_flag)

        The epsilon component is determined by bit `dim` of global_index.

        Returns:
            (final_sign, local_index, eps_flag)
        """
        Validation.basis_tuple(t1, allow_zero=True, allow_eps=True)
        Validation.basis_tuple(t2, allow_zero=True, allow_eps=True)

        dim = Validation.dimension(dim)

        s1, i = int(t1[0]), int(t1[1])
        s2, j = int(t2[0]), int(t2[1])

        # Zero propagates.
        if s1 == 0 or s2 == 0:
            return (0, 0, 0)

        total_size = 1 << (dim + 1)

        if i < 0 or i >= total_size:
            raise ValueError(f"index must be in [0, {total_size - 1}] for dual dim={dim}, got {i}")

        if j < 0 or j >= total_size:
            raise ValueError(f"index must be in [0, {total_size - 1}] for dual dim={dim}, got {j}")

        half = 1 << dim

        i_eps = i >= half
        j_eps = j >= half

        i_loc = i & (half - 1)
        j_loc = j & (half - 1)

        # --------------------------------------------------------------
        # Epsilon nilpotency:
        # (epsilon e_i) * (epsilon e_j) = 0
        # --------------------------------------------------------------
        if i_eps and j_eps:
            return (0, 0, 1)

        # --------------------------------------------------------------
        # Base multiplication.
        # --------------------------------------------------------------
        if self._split_mode:
            base_sign, local_idx = self._split.multiply_indices(i_loc, j_loc, dim)
        else:
            base_sign, local_idx = self._standard.multiply_indices(i_loc, j_loc)

        eps_flag = 1 if (i_eps or j_eps) else 0

        return (s1 * s2 * base_sign, local_idx, eps_flag)