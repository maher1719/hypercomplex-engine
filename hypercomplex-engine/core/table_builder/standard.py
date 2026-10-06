import numpy as np

from ..validation import Validation
from .common import index_dtype, real_table


class StandardTableBuilder:
    """
    Standard Cayley-Dickson table builder.

    OPMT Theorem 1.2:
        Block a: inherits parent
        Block b: transpose, index shifted
        Block c: conjugation on column j
        Block d: standard sign rules
    """

    def build(self, n: int):
        """
        Build the standard Cayley-Dickson table A_n.

        Returns
        -------
        signs:
            signs[i, j] = ±1
        indices:
            indices[i, j] = k, where e_i e_j = signs[i, j] * e_k
        """
        n = Validation.dimension(n)

        signs, indices = real_table()

        for _ in range(n):
            half = signs.shape[0]
            full = 2 * half
            dtype = index_dtype(full)

            new_signs = np.zeros((full, full), dtype=np.int8)
            new_indices = np.zeros((full, full), dtype=dtype)

            # ----------------------------------------------------------
            # Block a:
            # (e_i, 0)(e_j, 0) = (e_i e_j, 0)
            # ----------------------------------------------------------
            new_signs[:half, :half] = signs
            new_indices[:half, :half] = indices

            # ----------------------------------------------------------
            # Block b:
            # (e_i, 0)(0, e_j) = (0, e_j e_i)
            # ----------------------------------------------------------
            new_signs[:half, half:] = signs.T
            np.add(
                indices.T,
                dtype(half),
                out=new_indices[:half, half:],
                casting="unsafe",
            )

            # ----------------------------------------------------------
            # Block c:
            # (0, e_i)(e_j, 0) = (0, e_i e_j*)
            # Conjugation flips columns j > 0.
            # ----------------------------------------------------------
            block_c = signs.copy()
            block_c[:, 1:] = -block_c[:, 1:]

            new_signs[half:, :half] = block_c
            np.add(
                indices,
                dtype(half),
                out=new_indices[half:, :half],
                casting="unsafe",
            )

            # ----------------------------------------------------------
            # Block d:
            # (0, e_i)(0, e_j) = (-e_j* e_i, 0)
            # Start from -signs.T, then conjugation flips columns j > 0.
            # ----------------------------------------------------------
            block_d = -signs.T
            block_d[:, 1:] = -block_d[:, 1:]

            new_signs[half:, half:] = block_d
            new_indices[half:, half:] = indices.T

            signs = new_signs
            indices = new_indices

        return signs, indices