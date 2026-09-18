import numpy as np

from ..validation import Validation
from .common import index_dtype, real_table
from .standard import StandardTableBuilder


class SplitTableBuilder:
    """
    Split Cayley-Dickson table builder.

    Architecture:
        Standard parent A_{n-1} + one split doubling at the top.

    OPMT Theorem 2.1:
        Blocks a, b, c are identical to standard.
        Block d is the split Block d:
            diagonal = +1
            first column = +1
            first row j > 0 = -1
            interior = +sigma_a
    """

    def __init__(self):
        self._standard = StandardTableBuilder()

    def build(self, n: int):
        """
        Build the split Cayley-Dickson table A_n.

        Returns
        -------
        signs:
            signs[i, j] = ±1
        indices:
            indices[i, j] = k
        """
        n = Validation.dimension(n)

        if n == 0:
            return real_table()

        # Standard parent A_{n-1}
        parent_signs, parent_indices = self._standard.build(n - 1)

        half = parent_signs.shape[0]
        full = 2 * half
        dtype = index_dtype(full)

        signs = np.zeros((full, full), dtype=np.int8)
        indices = np.zeros((full, full), dtype=dtype)

        # --------------------------------------------------------------
        # Block a:
        # Inherits the standard parent.
        # --------------------------------------------------------------
        signs[:half, :half] = parent_signs
        indices[:half, :half] = parent_indices

        # --------------------------------------------------------------
        # Block b:
        # Same as standard.
        # --------------------------------------------------------------
        signs[:half, half:] = parent_signs.T
        np.add(
            parent_indices.T,
            dtype(half),
            out=indices[:half, half:],
            casting="unsafe",
        )

        # --------------------------------------------------------------
        # Block c:
        # Same as standard.
        # --------------------------------------------------------------
        block_c = parent_signs.copy()
        block_c[:, 1:] = -block_c[:, 1:]

        signs[half:, :half] = block_c
        np.add(
            parent_indices,
            dtype(half),
            out=indices[half:, :half],
            casting="unsafe",
        )

        # --------------------------------------------------------------
        # Block d:
        # Split Block d.
        #
        # (0, e_i)(0, e_j) = (+e_j* e_i, 0)
        #
        # No leading minus sign, but conjugation still flips columns j > 0.
        # --------------------------------------------------------------
        block_d = parent_signs.T.copy()
        block_d[:, 1:] = -block_d[:, 1:]

        signs[half:, half:] = block_d
        indices[half:, half:] = parent_indices.T

        return signs, indices