from ..validation import Validation
from .bit_utils import nu2, popcount   # kept: public helpers, other modules use them

# Unbound builtins bound once at import: one module-global lookup instead of a
# per-call attribute lookup on the int object.
_bl = int.bit_length
_bc = int.bit_count


class FastStandard:
    """
    Standard Cayley-Dickson O(1) multiplier.

    Input / output convention: (sign, index), index = i XOR j.

    Optimized evaluation (verified exhaustively against the reference form):
      * the break type is decided by (ti, tj) alone -- no comparison against t;
      * t = nu2(i XOR j) is computed only on the tie ti == tj (~1/3 of pairs);
      * the i-axis branch is collapsed to 2*b - 1 (Corollary cor:iaxis-collapse);
      * nu2 / popcount are inlined (their zero-guards are provably dead here).
    """

    def multiply(self, t1: tuple, t2: tuple) -> tuple:
        """Multiply two signed basis tuples (validated)."""
        Validation.basis_tuple(t1, allow_zero=False, allow_eps=False)
        Validation.basis_tuple(t2, allow_zero=False, allow_eps=False)
        s1, i = int(t1[0]), int(t1[1])
        s2, j = int(t2[0]), int(t2[1])
        sign, idx = self.multiply_indices(i, j)
        return (s1 * s2 * sign, idx)

    def multiply_indices(self, i: int, j: int) -> tuple:
        """
        Core O(1) multiplication for non-negative integer indices.
        Returns (sign, i XOR j).
        """
        if type(i) is not int:
            i = int(i)
        if type(j) is not int:
            j = int(j)
        # Explicit guard: do NOT fold into (i ^ j) < 0, which misses i == j < 0.
        if i < 0 or j < 0:
            raise ValueError("indices must be >= 0")

        idx = i ^ j
        if not (i and j):          # identity: e_0 * e_j = e_j
            return (1, idx)
        if i == j:                 # e_i^2 = -e_0 for i > 0
            return (-1, 0)

        ti = _bl(i & -i) - 1
        tj = _bl(j & -j) - 1

        if ti < tj:
            # t = ti < tj  ->  j-axis break
            k = tj + 1
            s = 1 - 2 * ((i >> tj) & 1)
        elif ti > tj:
            # t = tj < ti  ->  i-axis break; delta == -1 unconditionally
            k = ti + 1
            s = 2 * ((j >> ti) & 1) - 1
        else:
            # ti == tj  ->  t > ti  ->  diagonal break (only case needing t)
            t = _bl(idx & -idx) - 1
            k = t + 1
            s = 2 * ((i >> t) & 1) - 1

        if _bc((i | j) >> k) & 1:
            s = -s
        return (s, idx)