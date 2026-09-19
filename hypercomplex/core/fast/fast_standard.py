from ..validation import Validation
from .bit_utils import nu2, popcount


class FastStandard:
    """
    Standard Cayley-Dickson O(1) multiplier.

    Input / output convention:
        (sign, index)

    The sign is computed in O(1) Word-RAM time using:
        t  = nu2(i XOR j)
        ti = nu2(i)
        tj = nu2(j)

    The basis index is always:
        i XOR j
    """

    def multiply(self, t1: tuple, t2: tuple) -> tuple:
        """
        Multiply two standard basis element tuples.

        t1 = (sign1, index1)
        t2 = (sign2, index2)

        Returns:
            (final_sign, index1 XOR index2)
        """
        Validation.basis_tuple(t1, allow_zero=True, allow_eps=False)
        Validation.basis_tuple(t2, allow_zero=True, allow_eps=False)

        s1, i = int(t1[0]), int(t1[1])
        s2, j = int(t2[0]), int(t2[1])

        if s1 == 0 or s2 == 0:
            return (0, 0)

        sign, idx = self.multiply_indices(i, j)

        return (s1 * s2 * sign, idx)

    def multiply_indices(self, i: int, j: int) -> tuple:
        """
        Core O(1) multiplication for non-negative integer indices.

        Returns:
            (sign, i XOR j)
        """
        i = int(i)
        j = int(j)

        if i < 0 or j < 0:
            raise ValueError("indices must be >= 0")

        # Identity
        if i == 0 or j == 0:
            return (1, i ^ j)

        # Main diagonal
        if i == j:
            return (-1, 0)

        # Structural break levels
        t = nu2(i ^ j)
        ti = nu2(i)
        tj = nu2(j)

        # --------------------------------------------------------------
        # Diagonal break:
        #   t >= ti and t >= tj
        # --------------------------------------------------------------
        if t >= ti and t >= tj:
            k = t + 1
            s = 2 * ((i >> t) & 1) - 1

        # --------------------------------------------------------------
        # i-axis break:
        #   ti > t and ti >= tj
        # --------------------------------------------------------------
        elif ti >= tj:
            k = ti + 1

            # Quadrant bit:
            #   0 -> quadrant c
            #   1 -> quadrant d
            b = (j >> ti) & 1

            # Local j below the break level
            j_loc = j & ((1 << ti) - 1)

            # delta(j_loc):
            #   +1 if j_loc == 0
            #   -1 otherwise
            delta = 1 if j_loc == 0 else -1

            if b == 0:
                # quadrant c
                s = delta
            else:
                # quadrant d
                s = -delta

        # --------------------------------------------------------------
        # j-axis break:
        #   tj > t and tj > ti
        # --------------------------------------------------------------
        else:
            k = tj + 1

            # Quadrant bit:
            #   0 -> quadrant b, sign +1
            #   1 -> quadrant d, sign -1
            s = 1 - 2 * ((i >> tj) & 1)

        # --------------------------------------------------------------
        # Higher-bit parity correction:
        # standard flips in quadrants b, c, d
        # --------------------------------------------------------------
        if popcount((i | j) >> k) & 1:
            s = -s

        return (s, i ^ j)