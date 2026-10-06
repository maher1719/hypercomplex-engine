import numpy as np

from ..validation import Validation
from .common import index_dtype
from .standard import StandardTableBuilder
from .split import SplitTableBuilder


class DualTableBuilder:
    """
    Dual Cayley-Dickson table builder.

    Adds epsilon with epsilon^2 = 0.

    The returned table uses:
        signs[i, j]
        indices[i, j]   local parent index
        eps[i, j]       epsilon flag

    For dual tables, index is the local parent basis index, and eps_flag
    tells whether the product carries epsilon.

    This matches the tuple convention:
        (sign, index, eps_flag)
    """

    def __init__(self):
        self._standard = StandardTableBuilder()
        self._split = SplitTableBuilder()

    def build(self, n: int, split: bool = False):
        """
        Build the dual extension of A_n.

        Parameters
        ----------
        n:
            Parent algebra exponent.
            Parent dimension is 2^n.
            Dual table dimension is 2^(n+1).

        split:
            False -> dual over standard A_n
            True  -> dual over split A_n

        Returns
        -------
        signs, indices, eps
        """
        n = Validation.dimension(n)

        if split:
            base_signs, base_indices = self._split.build(n)
        else:
            base_signs, base_indices = self._standard.build(n)

        N = base_signs.shape[0]
        full = 2 * N
        dtype = index_dtype(full)

        signs = np.zeros((full, full), dtype=np.int8)
        indices = np.zeros((full, full), dtype=dtype)
        eps = np.zeros((full, full), dtype=np.uint8)

        # --------------------------------------------------------------
        # Block a:
        # e_i * e_j = parent product
        # --------------------------------------------------------------
        signs[:N, :N] = base_signs
        indices[:N, :N] = base_indices
        eps[:N, :N] = 0

        # --------------------------------------------------------------
        # Block b:
        # e_i * (epsilon e_j) = epsilon (e_i e_j)
        # --------------------------------------------------------------
        signs[:N, N:] = base_signs
        indices[:N, N:] = base_indices
        eps[:N, N:] = 1

        # --------------------------------------------------------------
        # Block c:
        # (epsilon e_i) * e_j = epsilon (e_i e_j)
        # --------------------------------------------------------------
        signs[N:, :N] = base_signs
        indices[N:, :N] = base_indices
        eps[N:, :N] = 1

        # --------------------------------------------------------------
        # Block d:
        # (epsilon e_i) * (epsilon e_j) = epsilon^2 (e_i e_j) = 0
        # --------------------------------------------------------------
        signs[N:, N:] = 0
        indices[N:, N:] = 0
        eps[N:, N:] = 0

        return signs, indices, eps