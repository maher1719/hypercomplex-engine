import bisect

class GradedSparseStandard:
    """
    Arbitrary-width O(1)-per-primitive multiplier on the GRADED (sparse
    set-of-bit-positions) representation. Memory-optimal: O(p) not O(B).

    Internal convention: an index m is a sorted tuple of 0-INDEXED bit
    positions, m = sum(2**b for b in pos). The Wilmot generator level is
    k = b + 1 (used only for display via graded_to_notation).
    """

    # ---------- representation: int <-> graded ----------
    @staticmethod
    def int_to_graded(m):
        """Python int -> sorted tuple of 0-indexed set-bit positions. O(p)."""
        m = int(m)
        pos = []
        while m:
            low = m & -m                    # isolate lowest set bit
            pos.append(low.bit_length() - 1)
            m ^= low                        # clear it
        return tuple(pos)                   # already ascending

    @staticmethod
    def graded_to_int(pos):
        """Sorted tuple of bit positions -> Python int. O(p)."""
        m = 0
        for b in pos:
            m |= (1 << b)
        return m

    @staticmethod
    def graded_to_notation(pos):
        """Graded form -> Wilmot string 'o{k1},{k2},...' (1-indexed)."""
        return "1" if not pos else "o" + ",".join(str(b + 1) for b in pos)

    # ---------- faithful graded primitives ----------
    @staticmethod
    def _nu2_gr(pos):
        """nu2 = lowest set bit = first position. Requires pos non-empty."""
        return pos[0]

    @staticmethod
    def _bit_gr(pos, t):
        """(m >> t) & 1 via membership test on sorted positions."""
        idx = bisect.bisect_left(pos, t)
        return 1 if (idx < len(pos) and pos[idx] == t) else 0

    @staticmethod
    def _nu2_xor_gr(pos_i, pos_j):
        """Lowest bit of i^j = lowest position in the symmetric difference."""
        a = b = 0
        while a < len(pos_i) and b < len(pos_j):
            if pos_i[a] < pos_j[b]:   return pos_i[a]
            elif pos_j[b] < pos_i[a]: return pos_j[b]
            else:                     a += 1; b += 1
        if a < len(pos_i): return pos_i[a]
        if b < len(pos_j): return pos_j[b]
        return None                       # i == j (handled upstream)

    @staticmethod
    def _xor_gr(pos_i, pos_j):
        """Full symmetric difference = graded form of i^j."""
        out = []; a = b = 0
        while a < len(pos_i) and b < len(pos_j):
            if pos_i[a] < pos_j[b]:   out.append(pos_i[a]); a += 1
            elif pos_j[b] < pos_i[a]: out.append(pos_j[b]); b += 1
            else:                     a += 1; b += 1
        out.extend(pos_i[a:]); out.extend(pos_j[b:])
        return tuple(out)

    @staticmethod
    def _pop_or_shift_gr(pos_i, pos_j, k):
        """popcount((i|j) >> k) = #positions in union(pos_i,pos_j) that are >= k."""
        count = 0; a = b = 0
        while a < len(pos_i) and b < len(pos_j):
            if pos_i[a] < pos_j[b]:
                if pos_i[a] >= k: count += 1
                a += 1
            elif pos_j[b] < pos_i[a]:
                if pos_j[b] >= k: count += 1
                b += 1
            else:
                if pos_i[a] >= k: count += 1      # union: count once
                a += 1; b += 1
        while a < len(pos_i):
            if pos_i[a] >= k: count += 1
            a += 1
        while b < len(pos_j):
            if pos_j[b] >= k: count += 1
            b += 1
        return count

    # ---------- the O(1) evaluator on graded forms ----------
    def multiply_graded(self, pos_i, pos_j):
        """Core evaluator on graded forms -> (sign, result_pos).
        Uses the i-axis collapse (Cor. cor:iaxis-collapse): delta = -1."""
        if not pos_i or not pos_j:                 # identity
            return (1, self._xor_gr(pos_i, pos_j))
        if pos_i == pos_j:                         # main diagonal
            return (-1, ())

        t  = self._nu2_xor_gr(pos_i, pos_j)
        ti = self._nu2_gr(pos_i)
        tj = self._nu2_gr(pos_j)
        k  = max(t, ti, tj) + 1

        if t >= ti and t >= tj:                    # diagonal break
            s = 2 * self._bit_gr(pos_i, t) - 1
        elif ti >= tj:                             # i-axis break; delta=-1 always
            b = self._bit_gr(pos_j, ti)
            s = -((-1) ** b)                       # = 2b - 1
        else:                                      # j-axis break
            s = 1 - 2 * self._bit_gr(pos_i, tj)

        if self._pop_or_shift_gr(pos_i, pos_j, k) & 1:
            s = -s
        return (s, self._xor_gr(pos_i, pos_j))

    def multiply_indices(self, i, j):
        """Convenience: ints in -> (sign, i^j) out. Mirrors FastStandard."""
        sign, result_pos = self.multiply_graded(
            self.int_to_graded(i), self.int_to_graded(j))
        return (sign, self.graded_to_int(result_pos))