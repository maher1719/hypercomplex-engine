from ..validation import Validation


class StandardHolographic:
    """
    Standard Cayley-Dickson holographic O(n) multiplier.

    Input / output convention:
        (sign, index)
    """

    def multiply(self, t1: tuple, t2: tuple) -> tuple:
        """
        Multiply two basis elements.

        t1 = (sign1, index1)
        t2 = (sign2, index2)

        Returns:
            (final_sign, index1 XOR index2)
        """
        Validation.basis_tuple(t1, allow_zero=False, allow_eps=False)
        Validation.basis_tuple(t2, allow_zero=False, allow_eps=False)

        s1, i = int(t1[0]), int(t1[1])
        s2, j = int(t2[0]), int(t2[1])

        sign, idx = self.multiply_indices(i, j)

        return (s1 * s2 * sign, idx)

    def multiply_indices(self, i: int, j: int) -> tuple:
        """
        Core O(n) descent for non-negative integer indices.

        Returns:
            (sign, i XOR j)
        """
        i = int(i)
        j = int(j)

        if i < 0 or j < 0:
            raise ValueError("indices must be >= 0")

        if i == 0 and j == 0:
            return (1, 0)

        sign = 1
        ic = i
        jc = j

        n = max(i, j).bit_length()

        for level in range(n, 0, -1):
            half = 1 << (level - 1)

            if ic < half and jc < half:
                q = "a"
            elif ic < half and jc >= half:
                q = "b"
                jc -= half
            elif ic >= half and jc < half:
                q = "c"
                ic -= half
            else:
                q = "d"
                ic -= half
                jc -= half

            # Structural position inside the current quadrant.
            if ic == 0 or jc == 0 or ic == jc:
                structural_sign = self._structural_sign(q, ic, jc, split=False)
                return (sign * structural_sign, i ^ j)

            # Non-structural sign flip.
            if self._flips(q, split=False):
                sign = -sign

        # Should normally terminate via structural detection.
        return (sign, i ^ j)

    @staticmethod
    def _flips(q: str, split: bool) -> bool:
        """
        Non-structural sign flip rule.

        Standard:
            flips in b, c, d

        Split:
            flips in b, c only
        """
        if split:
            return q in ("b", "c")
        return q in ("b", "c", "d")

    @staticmethod
    def _structural_sign(q: str, i_loc: int, j_loc: int, split: bool) -> int:
        """
        Structural sign tables from OPMT Theorem 1.2 and Theorem 2.1.
        """

        # --------------------------------------------------------------
        # Diagonal
        # --------------------------------------------------------------
        if i_loc == j_loc:
            if i_loc == 0:
                if q in ("a", "b", "c"):
                    return 1
                # Block d, (0,0)
                return 1 if split else -1

            # i_loc == j_loc > 0
            if q in ("a", "b"):
                return -1
            if q == "c":
                return 1
            # Block d diagonal
            return 1 if split else -1

        # --------------------------------------------------------------
        # First row: i_loc == 0, j_loc > 0
        # --------------------------------------------------------------
        if i_loc == 0:
            if q in ("a", "b"):
                return 1
            if q == "c":
                return -1
            # Block d first row
            return -1 if split else 1

        # --------------------------------------------------------------
        # First column: j_loc == 0, i_loc > 0
        # --------------------------------------------------------------
        if j_loc == 0:
            if q in ("a", "b", "c"):
                return 1
            # Block d first column
            return 1 if split else -1

        raise ValueError("Not a structural position.")