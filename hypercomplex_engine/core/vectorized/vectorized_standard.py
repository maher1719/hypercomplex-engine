"""
VectorizedStandard — uint64 SIMD lane-wise O(1) Cayley-Dickson evaluator.

Optimizations vs. the first release (all verified against FastStandard):
  * single `idx = i ^ j`, reused for the trailing-zero count (was 3 XOR passes)
  * nu2 guards zero internally, so the *_safe temporaries are gone
  * i-axis branch collapsed to 2*b - 1 (Corollary cor:iaxis-collapse):
    removes the variable shift `1 << ti`, the `j_loc` AND, and two np.where
  * explicit wide-shift guard for k == 64 (see _GUARD_WIDE_SHIFT)

Measured: 71.2 -> 39.4 ns/pair (guard on, 1.81x) / 34.7 ns/pair (guard off, 2.05x).
"""
import numpy as np

_U64 = np.uint64

# NumPy returns 0 for `x >> width` (verified on 2.2.6 for uint64/int64/uint32/
# int32/uint8, both array and scalar paths). That behaviour is not formally
# documented, and `>> 64` on a 64-bit word is UB in C/C++. The guard makes the
# k == 64 lane explicit and portable; it costs ~4.7 ns/pair. Set False to skip.
_GUARD_WIDE_SHIFT = True


def _get_np_popcount():
    """Hardware population count, with a SWAR fallback for NumPy < 2.0."""
    try:
        return np.bitwise_count
    except AttributeError:
        def _swar(x):
            x = np.asarray(x, dtype=np.uint64)
            x = x - ((x >> _U64(1)) & _U64(0x5555555555555555))
            x = (x & _U64(0x3333333333333333)) + ((x >> _U64(2)) & _U64(0x3333333333333333))
            x = (x + (x >> _U64(4))) & _U64(0x0F0F0F0F0F0F0F0F)
            return ((x * _U64(0x0101010101010101)) >> _U64(56)).astype(np.uint8)
        return _swar


_popcount = _get_np_popcount()


def _nu2_vec(x):
    """
    Lane-wise nu2 via popcount((x & -x) - 1).
    x == 0 maps to 0; those lanes are overridden by mask_zero downstream.
    """
    xs = np.where(x == 0, _U64(1), x)
    return _popcount((xs & -xs) - _U64(1)).astype(_U64)


class VectorizedStandard:
    """Standard Cayley-Dickson O(1) multiplier on uint64 SIMD lanes (indices <= 64 bits)."""

    @staticmethod
    def _get_np_popcount():
        return _popcount

    def multiply_indices_vec(self, i_arr, j_arr):
        """
        Lane-wise O(1) evaluation of Definition def:o1.

        Returns
        -------
        (signs, idx) : (N,) int8 in {-1,+1}, (N,) uint64 with idx = i ^ j
        """
        i = np.asarray(i_arr, dtype=np.uint64)
        j = np.asarray(j_arr, dtype=np.uint64)

        idx = i ^ j                                    # single XOR, reused below
        mask_zero = (i == 0) | (j == 0)
        mask_diag = (i == j) & ~mask_zero

        t  = _nu2_vec(idx)
        ti = _nu2_vec(i)
        tj = _nu2_vec(j)
        k = np.maximum(np.maximum(t, ti), tj) + _U64(1)

        cond1 = (t >= ti) & (t >= tj)                  # diagonal break
        cond2 = (~cond1) & (ti >= tj)                  # i-axis break

        s1 = 2 * ((i >> t)  & _U64(1)).astype(np.int8) - 1
        # cond2 forces ti > tj strictly (ti == tj would make cond1 fire), so
        # j_loc > 0 and delta == -1 unconditionally: s2 = 2*b - 1.
        s2 = 2 * ((j >> ti) & _U64(1)).astype(np.int8) - 1
        s3 = 1 - 2 * ((i >> tj) & _U64(1)).astype(np.int8)

        s = np.where(cond1, s1, np.where(cond2, s2, s3))

        ij_or = i | j
        if _GUARD_WIDE_SHIFT:
            # k == 64 means no bits above the word: parity must be 0.
            parity = np.where(
                k < _U64(64),
                _popcount(ij_or >> np.minimum(k, _U64(63))) & np.uint8(1),
                np.uint8(0),
            )
        else:
            parity = _popcount(ij_or >> k) & np.uint8(1)

        s = np.where(parity == 1, -s, s)
        s = np.where(mask_zero, 1, s)
        s = np.where(mask_diag, -1, s)
        return s.astype(np.int8), idx

    def multiply_arrays(self, arr1, arr2):
        """
        Sparse twisted convolution of two real coefficient arrays
        (Eq. eq:sparse-conv): support extraction -> active grid -> lane-wise
        signs -> coefficient weighting -> scatter-add.

        NOTE: only the sign kernel above was optimized; this combiner is
        unchanged. The real-part factorization x0*y0 + x0*yt + y0*xt is a
        separate future optimization.
        """
        arr1 = np.asarray(arr1, dtype=np.float64)
        arr2 = np.asarray(arr2, dtype=np.float64)
        size1, size2 = len(arr1), len(arr2)
        max_size = max(size1, size2)
        if max_size <= 1:
            if size1 and size2:
                return np.array([arr1[0] * arr2[0]], dtype=np.float64)
            return np.zeros(1, dtype=np.float64)

        dim = 1 << (max_size - 1).bit_length()
        nz_A = np.nonzero(arr1)[0]
        nz_B = np.nonzero(arr2)[0]
        if len(nz_A) == 0 or len(nz_B) == 0:
            return np.zeros(dim, dtype=np.float64)

        I, J = np.meshgrid(nz_A, nz_B, indexing='ij')
        S, V = self.multiply_indices_vec(I, J)
        M = np.outer(arr1[nz_A], arr2[nz_B])
        weighted = S.astype(np.float64) * M

        result = np.zeros(dim, dtype=np.float64)
        np.add.at(result, V.ravel(), weighted.ravel())   # unbuffered scatter-add
        return result