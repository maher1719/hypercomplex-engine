import numpy as np



class CDTableBuilder:
    """
    Builds full multiplication tables for Cayley–Dickson algebras using the
    OPMT (Ordered-Pair Multiplication Table) block decomposition.

    The doubling formula (a, b)(c, d) = (ac - d*b, da + bc*) decomposes into
    four blocks per Theorem 1.2 / Theorem 2.1 of the OPMT sign law:

        Block a: (e_i, 0)(e_j, 0) = (e_i e_j, 0)
        Block b: (e_i, 0)(0, e_j) = (0, e_j e_i)
        Block c: (0, e_i)(e_j, 0) = (0, e_i e_j*)
        Block d: (0, e_i)(0, e_j) = (-e_j* e_i, 0)   [standard]
                 (0, e_i)(0, e_j) = (+e_j* e_i, 0)   [split]

    Returns arrays: signs[i,j] = ±1, indices[i,j] = basis index k,
    so that e_i * e_j = signs[i,j] * e_{indices[i,j]}.
    """

    @staticmethod
    def standard(n: int):
        """
        Standard Cayley–Dickson algebra A_n of dimension 2^n.
        OPMT Theorem 1.2 (standard sign law).

        n=0: reals, n=1: complex, n=2: quaternions, n=3: octonions.
        """
        signs = np.array([[1]], dtype=np.int8)
        indices = np.array([[0]], dtype=np.int64)

        for _ in range(1, n + 1):
            half = signs.shape[0]
            full = 2 * half
            new_signs = np.zeros((full, full), dtype=np.int8)
            new_indices = np.zeros((full, full), dtype=np.int64)

            # ------------------------------------------------------------------
            # BLOCK A: (e_i, 0)(e_j, 0) = (e_i e_j, 0)
            # OPMT Theorem 1.2, Eq. (1): Block a inherits the parent algebra A_n.
            # All structured and interior signs come directly from the parent.
            # ------------------------------------------------------------------
            new_signs[:half, :half] = signs
            new_indices[:half, :half] = indices

            # ------------------------------------------------------------------
            # BLOCK B: (e_i, 0)(0, e_j) = (0, e_j e_i)
            # OPMT Theorem 1.2, Eq. (2): product order is REVERSED (e_j e_i),
            # so we use the transpose. Index shifts into the upper half (+half).
            # Interior sign = -σ_a, which arises automatically from
            # anti-commutativity (e_j e_i = -e_i e_j for i≠j, i,j>0).
            # Structured diagonal/row/col signs also emerge correctly from
            # the transpose without explicit overrides.
            # ------------------------------------------------------------------
            new_signs[:half, half:] = signs.T
            new_indices[:half, half:] = indices.T + half

            # ------------------------------------------------------------------
            # BLOCK C: (0, e_i)(e_j, 0) = (0, e_i e_j*)
            # OPMT Theorem 1.2, Eq. (3): conjugation applies to e_j.
            # Conjugation rule: e_0* = e_0, e_k* = -e_k for k>0.
            # Therefore we flip the sign for all COLUMNS j>0 (the conjugated
            # factor is the column index j, not the row index i).
            # Interior sign = -σ_a for j>0. Diagonal (j=i>0) becomes +1
            # because -σ_a(i,i) = -(-1) = +1, matching Theorem 1.2's rule
            # that Block c diagonal is always positive.
            # ------------------------------------------------------------------
            s_ll = signs.copy()
            s_ll[:, 1:] = -s_ll[:, 1:]  # conjugation: flip columns j>0
            new_signs[half:, :half] = s_ll
            new_indices[half:, :half] = indices + half

            # ------------------------------------------------------------------
            # BLOCK D: (0, e_i)(0, e_j) = (-e_j* e_i, 0)
            # OPMT Theorem 1.2, Eq. (4): conjugation on e_j, then negate.
            # Start with -signs.T (the leading minus from the doubling formula
            # combined with the reversed product order e_j e_i).
            # Then apply conjugation: for columns j>0, e_j* = -e_j introduces
            # a second sign flip, so we negate those columns again.
            # Conjugation applies to COLUMNS (j index), not rows.
            # Result: diagonal i=j=0 -> -1, diagonal i=j>0 -> -1,
            #         first row (i=0, j>0) -> +1, first col (j=0, i>0) -> -1,
            # all matching Theorem 1.2's Block d structured signs.
            # ------------------------------------------------------------------
            s_lr = -signs.T
            s_lr[:, 1:] = -s_lr[:, 1:]  # conjugation: flip columns j>0
            new_signs[half:, half:] = s_lr
            new_indices[half:, half:] = indices.T

            signs, indices = new_signs, new_indices

        return signs, indices

    @staticmethod
    def split(n: int):
        """
        Split Cayley–Dickson algebra of dimension 2^n.
        OPMT Theorem 2.1 (split sign law).

        Uses the split doubling formula: (a,b)(c,d) = (ac + d*b, da + bc*).
        The ONLY difference from standard is Block d:
          - Block d interior sign = +σ_a (not -σ_a).
          - Block d structured signs are inverted vs standard.
        Blocks a, b, c are identical to the standard construction.
        """
        signs = np.array([[1]], dtype=np.int8)
        indices = np.array([[0]], dtype=np.int64)

        for _ in range(1, n + 1):
            half = signs.shape[0]
            full = 2 * half
            new_signs = np.zeros((full, full), dtype=np.int8)
            new_indices = np.zeros((full, full), dtype=np.int64)

            # ------------------------------------------------------------------
            # BLOCK A: inherits parent (identical to standard).
            # OPMT Theorem 2.1, Eq. (5).
            # ------------------------------------------------------------------
            new_signs[:half, :half] = signs
            new_indices[:half, :half] = indices

            # ------------------------------------------------------------------
            # BLOCK B: (e_i, 0)(0, e_j) = (0, e_j e_i)
            # OPMT Theorem 2.1, Eq. (6): identical to standard Block b.
            # ------------------------------------------------------------------
            new_signs[:half, half:] = signs.T
            new_indices[:half, half:] = indices.T + half

            # ------------------------------------------------------------------
            # BLOCK C: (0, e_i)(e_j, 0) = (0, e_i e_j*)
            # OPMT Theorem 2.1, Eq. (7): identical to standard Block c.
            # Conjugation flips columns j>0.
            # ------------------------------------------------------------------
            s_ll = signs.copy()
            s_ll[:, 1:] = -s_ll[:, 1:]  # conjugation: flip columns j>0
            new_signs[half:, :half] = s_ll
            new_indices[half:, :half] = indices + half

            # ------------------------------------------------------------------
            # BLOCK D: (0, e_i)(0, e_j) = (+e_j* e_i, 0)
            # OPMT Theorem 2.1, Eq. (8): the KEY DIFFERENCE from standard.
            # The split doubling has +d*b instead of -d*b, so there is NO
            # leading negation. We start with +signs.T (not -signs.T).
            # Conjugation still flips columns j>0.
            # Result: diagonal always +1, first col +1, first row (j>0) -1,
            # interior = +σ_a. This is the complete sign inversion of Block d
            # relative to standard, per Corollary 3.10.1.
            # ------------------------------------------------------------------
            s_lr = signs.T.copy()
            s_lr[:, 1:] = -s_lr[:, 1:]  # conjugation: flip columns j>0
            new_signs[half:, half:] = s_lr
            new_indices[half:, half:] = indices.T

            signs, indices = new_signs, new_indices

        return signs, indices

    @staticmethod
    def dual(n: int, split: bool = False):
        """
        Dual extension of A_n: adjoins ε with ε² = 0.
        Dimension doubles from 2^n to 2^(n+1).

        This is NOT part of the OPMT block theorem (Theorems 1.2 / 2.1);
        it is a separate central-nilpotent extension. The basis becomes
        {e_0,...,e_{N-1}, εe_0,...,εe_{N-1}} where N = 2^n.

        Multiplication rules (ε is central, ε² = 0):
            e_i · e_j       = (parent product)          -> no ε
            e_i · (εe_j)    = ε(e_i e_j)                -> ε factor
            (εe_i) · e_j    = ε(e_i e_j)                -> ε factor
            (εe_i) · (εe_j) = ε²(e_i e_j) = 0           -> zero

        Returns: signs, indices (local parent indices), eps (0/1 flag).
        The index is the LOCAL parent index (not shifted), because the ε
        factor is tracked separately in the eps array.
        """
        if split:
            base_signs, base_indices = CDTableBuilder.split(n)
        else:
            base_signs, base_indices = CDTableBuilder.standard(n)

        N = base_signs.shape[0]
        full = 2 * N

        signs = np.zeros((full, full), dtype=np.int8)
        indices = np.zeros((full, full), dtype=np.int64)
        eps = np.zeros((full, full), dtype=np.int8)

        # ------------------------------------------------------------------
        # BLOCK A (lower-left): e_i · e_j = parent product, no ε.
        # Inherits the parent algebra directly.
        # ------------------------------------------------------------------
        signs[:N, :N] = base_signs
        indices[:N, :N] = base_indices
        # eps[:N, :N] = 0  (already zero)

        # ------------------------------------------------------------------
        # BLOCK B (upper-right): e_i · (εe_j) = ε(e_i e_j).
        # ε is central, so the product is ε times the parent product.
        # Index is the LOCAL parent index; eps flag marks the ε factor.
        # ------------------------------------------------------------------
        signs[:N, N:] = base_signs
        indices[:N, N:] = base_indices
        eps[:N, N:] = 1

        # ------------------------------------------------------------------
        # BLOCK C (lower-left of upper half): (εe_i) · e_j = ε(e_i e_j).
        # Same as Block b because ε is central (commutes with all elements).
        # ------------------------------------------------------------------
        signs[N:, :N] = base_signs
        indices[N:, :N] = base_indices
        eps[N:, :N] = 1

        # ------------------------------------------------------------------
        # BLOCK D (upper-right): (εe_i) · (εe_j) = ε²(e_i e_j) = 0.
        # Nilpotency kills the entire block. signs, indices, eps all zero.
        # ------------------------------------------------------------------
        # signs[N:, N:] = 0  (already zero)
        # indices[N:, N:] = 0  (already zero)
        # eps[N:, N:] = 0  (already zero)

        return signs, indices, eps