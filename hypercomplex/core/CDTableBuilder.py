import numpy as np




class CDTableBuilder:

    # ==================================================================
    # 1. TABLE CONSTRUCTION
    #    (UNCHANGED - faithful to OPMT)
    # ==================================================================

    @staticmethod
    def standard(n):
        signs = np.array([[1]], dtype=np.int8)
        indices = np.array([[0]], dtype=np.int64)
        
        for _ in range(1, n + 1):
            half = signs.shape[0]
            full = 2 * half
            new_signs = np.zeros((full, full), dtype=np.int8)
            new_indices = np.zeros((full, full), dtype=np.int64)
            
            # Upper-Left Block (li, lj)
            new_signs[:half, :half] = signs
            new_indices[:half, :half] = indices
            
            # Upper-Right Block (lj, li) + half
            new_signs[:half, half:] = signs.T
            new_indices[:half, half:] = indices.T + half
            
            # Lower-Left Block (li, lj) with sign flipping for lj > 0
            s_ll = signs.copy()
            s_ll[:, 1:] = -s_ll[:, 1:]  # Flip signs where col index (j) > 0
            new_signs[half:, :half] = s_ll
            new_indices[half:, :half] = indices + half
            
            # Lower-Right Block (lj, li) with systematic sign flipping
            s_lr = -signs.T
            s_lr[:, 1:] = -s_lr[:, 1:]  # FIX: Flip signs where col index (j) > 0
            new_signs[half:, half:] = s_lr
            new_indices[half:, half:] = indices.T
            
            signs, indices = new_signs, new_indices
            
        return signs, indices
    @staticmethod
    def split(n):
        if n == 0:
            return (
                np.array([[1]], dtype=np.int8),
                np.array([[0]], dtype=np.int64),
            )

        parent_signs, parent_indices = CDTableBuilder.standard(n - 1)

        half = parent_signs.shape[0]
        full = 2 * half

        signs = np.zeros((full, full), dtype=np.int8)
        indices = np.zeros((full, full), dtype=np.int64)

        for i in range(full):
            for j in range(full):
                row_high = i >= half
                col_high = j >= half

                li = i % half
                lj = j % half

                if not row_high and not col_high:
                    signs[i, j] = parent_signs[li, lj]
                    indices[i, j] = parent_indices[li, lj]

                elif not row_high and col_high:
                    signs[i, j] = parent_signs[lj, li]
                    indices[i, j] = parent_indices[lj, li] + half

                elif row_high and not col_high:
                    s = parent_signs[li, lj]

                    if lj > 0:
                        s = -s

                    signs[i, j] = s
                    indices[i, j] = parent_indices[li, lj] + half

                else:
                    s = parent_signs[lj, li]

                    if lj > 0:
                        s = -s

                    signs[i, j] = s
                    indices[i, j] = parent_indices[lj, li]

        idx = np.arange(full, dtype=np.int64)

        signs[:, 0] = 1
        indices[:, 0] = idx

        signs[0, :] = 1
        indices[0, :] = idx

        return signs, indices

    @staticmethod
    def dual(n, split=False):
        if split:
            base_signs, base_indices = CDTableBuilder.split(n)
        else:
            base_signs, base_indices = CDTableBuilder.standard(n)

        N = base_signs.shape[0]
        full = 2 * N

        signs = np.zeros((full, full), dtype=np.int8)
        indices = np.zeros((full, full), dtype=np.int64)
        eps = np.zeros((full, full), dtype=np.int8)

        for i in range(full):
            for j in range(full):
                i_eps = i >= N
                j_eps = j >= N

                li = i % N
                lj = j % N

                if i_eps and j_eps:
                    signs[i, j] = 0
                    indices[i, j] = 0
                    eps[i, j] = 0

                elif i_eps or j_eps:
                    signs[i, j] = base_signs[li, lj]
                    indices[i, j] = base_indices[li, lj]
                    eps[i, j] = 1

                else:
                    signs[i, j] = base_signs[i, j]
                    indices[i, j] = base_indices[i, j]
                    eps[i, j] = 0

        return signs, indices, eps

    # ==================================================================
    # 2. SMALL FAITHFUL ACCESSORS / VALIDATION
    # ==================================================================

    @staticmethod
    def basis_product(signs, indices, i, j):
        """
        Returns (sign, index) for e_i * e_j.
        """
        return int(signs[i, j]), int(indices[i, j])

    @staticmethod
    def validate(signs, indices):
        """
        OPMT-based sanity checks (read-only, no construction changes):

        1. e0 identity row and column
        2. index law k = i XOR j for nonzero products
        3. anti-commutativity for i, j > 0, i != j
        """
        dim = signs.shape[0]
        ar = np.arange(dim, dtype=np.int64)

        if not np.all(signs[0, :] == 1):
            return False
        if not np.all(signs[:, 0] == 1):
            return False
        if not np.all(indices[0, :] == ar):
            return False
        if not np.all(indices[:, 0] == ar):
            return False

        i_grid = ar[:, None]
        j_grid = ar[None, :]

        nz = signs != 0
        if not np.array_equal(indices[nz], (i_grid ^ j_grid)[nz]):
            return False

        off = (i_grid > 0) & (j_grid > 0) & (i_grid != j_grid)
        if not np.all(signs[off] == -signs.T[off]):
            return False

        return True

    