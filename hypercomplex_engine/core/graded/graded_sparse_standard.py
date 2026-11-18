"""
GradedSparseStandard — arbitrary-width O(p) evaluator on the graded form.

An index m is a sorted tuple of 0-INDEXED set-bit positions, m = sum(2**b).
The Wilmot generator level is k = b + 1 (display only, graded_to_notation).

Optimizations vs. the first release (verified against FastStandard):
  * SINGLE fused two-pointer merge replaces three separate O(p) scans
    (_nu2_xor_gr + _pop_or_shift_gr + _xor_gr). Valid because every element
    before the first difference is common to both lists and lies below t < k,
    so it contributes nothing to the union count.
  * the diagonal-break bit comes free from the merge (no bisect)
  * i-axis branch uses 2*b - 1 (no integer exponentiation)
  * hybrid int_to_graded: peel a few low bits with m & -m (wins for very
    sparse indices), then a C-level byte scan for the remainder.

Measured: 1.8x-2.0x on multiply_graded; up to 9x on wide int_to_graded.
The legacy helpers are retained for reference and backward compatibility.
"""
import bisect
import numpy as np

# Peel budget for int_to_graded, scaled by width: peeling one bit costs
# O(B/w) big-int work, the byte scan costs O(B/8) C work plus O(p) Python.
_PEEL_NUM = 500_000
_PEEL_CAP = 24


class GradedSparseStandard:
    """Arbitrary-width sparse (graded-form) standard Cayley-Dickson multiplier."""

    # ---------- representation ----------
    @staticmethod
    def int_to_graded(m):
        """Python int -> sorted tuple of 0-indexed set-bit positions."""
        m = int(m)
        if m == 0:
            return ()
        pos = []
        budget = max(1, min(_PEEL_CAP, _PEEL_NUM // max(1, m.bit_length())))
        for _ in range(budget):
            if not m:
                return tuple(pos)
            lo = m & -m
            pos.append(lo.bit_length() - 1)
            m ^= lo
        if not m:
            return tuple(pos)
        # dense remainder: C-level scan, then per-nonzero-byte extraction
        data = m.to_bytes((m.bit_length() + 7) // 8, 'little')
        arr = np.frombuffer(data, dtype=np.uint8)
        for bi in np.flatnonzero(arr):
            v = int(arr[bi])
            base = int(bi) * 8
            while v:
                lo = v & -v
                pos.append(base + lo.bit_length() - 1)
                v ^= lo
        return tuple(pos)

    @staticmethod
    def graded_to_int(pos):
        """Sorted tuple of bit positions -> Python int."""
        m = 0
        for b in pos:
            m |= 1 << b
        return m

    @staticmethod
    def graded_to_notation(pos):
        """Graded form -> Wilmot string 'o{k1},{k2},...' (1-indexed levels)."""
        return "1" if not pos else "o" + ",".join(str(b + 1) for b in pos)

    # ---------- primitives (kept for reference / external use) ----------
    @staticmethod
    def _nu2_gr(pos):
        return pos[0]

    @staticmethod
    def _bit_gr(pos, t):
        q = bisect.bisect_left(pos, t)
        return 1 if (q < len(pos) and pos[q] == t) else 0

    @staticmethod
    def _nu2_xor_gr(pos_i, pos_j):
        a = b = 0
        while a < len(pos_i) and b < len(pos_j):
            if pos_i[a] < pos_j[b]:
                return pos_i[a]
            if pos_j[b] < pos_i[a]:
                return pos_j[b]
            a += 1; b += 1
        if a < len(pos_i):
            return pos_i[a]
        if b < len(pos_j):
            return pos_j[b]
        return None

    @staticmethod
    def _xor_gr(pos_i, pos_j):
        out = []; a = b = 0
        while a < len(pos_i) and b < len(pos_j):
            if pos_i[a] < pos_j[b]:
                out.append(pos_i[a]); a += 1
            elif pos_j[b] < pos_i[a]:
                out.append(pos_j[b]); b += 1
            else:
                a += 1; b += 1
        out.extend(pos_i[a:]); out.extend(pos_j[b:])
        return tuple(out)

    @staticmethod
    def _pop_or_shift_gr(pos_i, pos_j, k):
        cnt = 0; a = b = 0
        while a < len(pos_i) and b < len(pos_j):
            if pos_i[a] < pos_j[b]:
                cnt += (pos_i[a] >= k); a += 1
            elif pos_j[b] < pos_i[a]:
                cnt += (pos_j[b] >= k); b += 1
            else:
                cnt += (pos_i[a] >= k); a += 1; b += 1
        while a < len(pos_i):
            cnt += (pos_i[a] >= k); a += 1
        while b < len(pos_j):
            cnt += (pos_j[b] >= k); b += 1
        return cnt

    # ---------- fused evaluator ----------
    def multiply_graded(self, pos_i, pos_j):
        """
        Core O(p) evaluator on graded forms, single fused merge.

        Returns
        -------
        (sign, result_pos) : int in {-1,+1}, sorted tuple (graded form of i ^ j)
        """
        ni, nj = len(pos_i), len(pos_j)
        if ni == 0:
            return (1, pos_j)                      # e_0 * e_j = e_j
        if nj == 0:
            return (1, pos_i)
        if pos_i == pos_j:
            return (-1, ())                        # e_i^2 = -e_0  (i > 0)

        ti, tj = pos_i[0], pos_j[0]
        a = b = 0
        # Common prefix: all these positions are shared and < t < k, so they
        # contribute nothing to the union count. Skipping them yields t at once
        # (zero iterations when ti != tj, since then t = min(ti, tj)).
        while a < ni and b < nj and pos_i[a] == pos_j[b]:
            a += 1; b += 1
        if a < ni and b < nj:
            if pos_i[a] < pos_j[b]:
                t, t_from_i = pos_i[a], 1
            else:
                t, t_from_i = pos_j[b], 0
        elif a < ni:
            t, t_from_i = pos_i[a], 1
        else:
            t, t_from_i = pos_j[b], 0

        k = (t if t > ti else ti)
        k = (k if k > tj else tj)
        k += 1

        out = []
        cnt = 0
        while a < ni and b < nj:
            x, y = pos_i[a], pos_j[b]
            if x < y:
                out.append(x); cnt += (x >= k); a += 1
            elif y < x:
                out.append(y); cnt += (y >= k); b += 1
            else:
                cnt += (x >= k); a += 1; b += 1     # union: count once
        while a < ni:
            x = pos_i[a]; out.append(x); cnt += (x >= k); a += 1
        while b < nj:
            y = pos_j[b]; out.append(y); cnt += (y >= k); b += 1

        if t >= ti and t >= tj:
            s = 2 * t_from_i - 1                     # bit known from the merge
        elif ti >= tj:
            s = 2 * self._bit_gr(pos_j, ti) - 1      # delta == -1 (collapsed)
        else:
            s = 1 - 2 * self._bit_gr(pos_i, tj)

        if cnt & 1:
            s = -s
        return (s, tuple(out))

    # ---------- convenience ----------
    def multiply_indices(self, i, j):
        """
        Python ints in -> (sign, i ^ j) out. Mirrors FastStandard.

        NOTE: this wrapper pays O(p * B/w) for the int -> graded conversion.
        The O(p) claim applies to multiply_graded on graded-form inputs.
        """
        sign, pos = self.multiply_graded(self.int_to_graded(i), self.int_to_graded(j))
        return (sign, self.graded_to_int(pos))