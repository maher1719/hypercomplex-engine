from ..validation import Validation
from .fast_standard import FastStandard
from .fast_split import FastSplit

class FastDual:
    """
    Dual Cayley-Dickson O(1) multiplier.

    STRICT 3-TUPLE INPUT:
        (sign, local_index, eps_flag)
        
    local_index must be strictly < 2^dim.
    eps_flag must be 0 or 1.
    """

    def __init__(self, split: bool = False):
        self._split_mode = bool(split)
        self._standard = FastStandard()
        self._split = FastSplit() if self._split_mode else None

    def multiply(self, t1: tuple, t2: tuple, dim: int) -> tuple:
        dim = Validation.dimension(dim)
        Validation.basis_tuple(t1, allow_zero=True, allow_eps=True)
        Validation.basis_tuple(t2, allow_zero=True, allow_eps=True)

        # STRICT 3-TUPLE UNPACKING
        s1, i, i_eps = int(t1[0]), int(t1[1]), int(t1[2])
        s2, j, j_eps = int(t2[0]), int(t2[1]), int(t2[2])

        if s1 == 0 or s2 == 0:
            return (0, 0, 0)

        sign, idx, eps = self.multiply_indices(i, j, dim, i_eps, j_eps)
        return (s1 * s2 * sign, idx, eps)

    def multiply_indices(self, i: int, j: int, dim: int, i_eps_flag: int = 0, j_eps_flag: int = 0) -> tuple:
        dim = Validation.dimension(dim)
        half = 1 << dim
        
        # STRICT LOCAL ENFORCEMENT: No more global overflow!
        if i < 0 or i >= half:
            raise ValueError(f"local index must be in [0, {half - 1}] for dual dim={dim}, got {i}")
        if j < 0 or j >= half:
            raise ValueError(f"local index must be in [0, {half - 1}] for dual dim={dim}, got {j}")
            
        i_eps = bool(i_eps_flag)
        j_eps = bool(j_eps_flag)

        # Nilpotency: eps^2 = 0
        if i_eps and j_eps:
            return (0, 0, 0)

        # Base multiplication (parent algebra)
        if self._split_mode:
            base_sign, local_idx = self._split.multiply_indices(i, j, dim)
        else:
            base_sign, local_idx = self._standard.multiply_indices(i, j)

        eps_flag = 1 if (i_eps or j_eps) else 0
        return (base_sign, local_idx, eps_flag)