from ..validation import Validation
from .fast_standard import FastStandard
from .fast_split import FastSplit


class FastDual:
    """
    Dual Cayley-Dickson O(1) multiplier.

    Input index convention:
        Total dual dimension is 2^(dim + 1).
        Bit `dim` marks the epsilon component.

        lower half: base elements
        upper half: epsilon * base elements

    Output:
        (sign, local_index, eps_flag)

    The output index is the local base index, so it can be formatted
    directly by CDFormat.

    Zero products are normalized to:
        (0, 0, 0)
    """

    def __init__(self, split: bool = False):
        self._split_mode = bool(split)
        self._standard = FastStandard()
        self._split = FastSplit() if self._split_mode else None

    # ==================================================================
    # Index helpers
    # ==================================================================

    @staticmethod
    def local_to_global(index: int, eps: int, dim: int) -> int:
        """
        Convert local dual representation to global index.

        local:
            index in [0, 2^dim - 1]
            eps = 0 or 1

        global:
            index in [0, 2^(dim+1) - 1]
        """
        dim = Validation.dimension(dim)
        index = int(index)
        eps = int(eps)

        half = 1 << dim

        if index < 0 or index >= half:
            raise ValueError(f"local index must be in [0, {half - 1}], got {index}")

        if eps:
            return index + half

        return index

    @staticmethod
    def global_to_local(global_index: int, dim: int) -> tuple:
        """
        Convert global dual index to:
            (local_index, eps_flag)
        """
        dim = Validation.dimension(dim)
        global_index = int(global_index)

        half = 1 << dim
        total = half << 1

        if global_index < 0 or global_index >= total:
            raise ValueError(f"global index must be in [0, {total - 1}], got {global_index}")

        return (global_index & (half - 1), 1 if global_index >= half else 0)

    # ==================================================================
    # Tuple input helper
    # ==================================================================

    def _extract_global(self, t: tuple, dim: int) -> tuple:
        """
        Accepts:
            (sign, global_index)
            (sign, local_index, eps_flag)

        Returns:
            (sign, global_index)
        """
        Validation.basis_tuple(t, allow_zero=True, allow_eps=True)

        sign = int(t[0])
        idx = int(t[1])

        if sign == 0:
            return (0, 0)

        half = 1 << Validation.dimension(dim)

        if len(t) >= 3:
            eps = int(t[2])

            # If the tuple carries an epsilon flag and the index is still
            # local, promote it to the global upper half.
            if eps == 1 and idx < half:
                idx += half

        return (sign, idx)

    # ==================================================================
    # Public API
    # ==================================================================

    def multiply(self, t1: tuple, t2: tuple, dim: int) -> tuple:
        """
        Multiply two dual basis element tuples.

        Accepted inputs:
            (sign, global_index)
            (sign, local_index, eps_flag)

        Returns:
            (final_sign, local_index, eps_flag)
        """
        dim = Validation.dimension(dim)
        Validation.basis_tuple(t1, allow_zero=True, allow_eps=True)
        Validation.basis_tuple(t2, allow_zero=True, allow_eps=True)

        

        s1, i = int(t1[0]), int(t1[1])
        s2, j = int(t2[0]), int(t2[1])
        i_eps_flag = int(t1[2]) if len(t1) == 3 else 0
        j_eps_flag = int(t2[2]) if len(t2) == 3 else 0

        if s1 == 0 or s2 == 0:
            return (0, 0, 0)

        sign, idx, eps = self.multiply_indices(i, j,dim,i_eps_flag,j_eps_flag)

        return (s1 * s2 * sign, idx, eps)

    def multiply_indices(self, i: int, j: int,dim:int, i_eps_flag=0 ,j_eps_flag=0) -> tuple:
        """
        Core dual O(1) multiplication using global indices.

        Returns:
            (sign, local_index, eps_flag)
        """
        i = int(i)
        j = int(j)
        dim = Validation.dimension(dim)

        half = 1 << dim
        total = half << 1

        if i < 0 or i >= total:
            raise ValueError(f"index must be in [0, {total - 1}] for dual dim={dim}, got {i}")

        if j < 0 or j >= total:
            raise ValueError(f"index must be in [0, {total - 1}] for dual dim={dim}, got {j}")

        i_eps = i >= half 
        j_eps = j >= half 

        i_loc = i & (half - 1)
        j_loc = j & (half - 1)

        # --------------------------------------------------------------
        # Epsilon nilpotency:
        # (epsilon e_i)(epsilon e_j) = 0
        # --------------------------------------------------------------
        if i_eps and j_eps:
            return (0, 0, 0)

        # --------------------------------------------------------------
        # Base multiplication
        # --------------------------------------------------------------
        if self._split_mode:
            base_sign, local_idx = self._split.multiply_indices(i_loc, j_loc, dim)
        else:
            base_sign, local_idx = self._standard.multiply_indices(i_loc, j_loc)

        eps_flag = 1 if (i_eps or j_eps) else 0

        return (base_sign, local_idx, eps_flag)