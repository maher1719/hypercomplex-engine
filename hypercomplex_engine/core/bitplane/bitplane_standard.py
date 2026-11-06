import numpy as np


class BitPlaneStandard:
    """
    Arbitrary-width SIMD (bit-sliced) standard Cayley-Dickson O(1) multiplier.

    Breaks the 64-bit SIMD-lane ceiling of VectorizedStandard by representing
    each index as a column of bits (bit-slicing). Every word operation of the
    O(1) evaluator (Definition def:o1) lifts to a lane-parallel boolean
    operation over the bit-planes.

    Representation:
        A batch of N indices of width B is stored as a bit-plane matrix of
        shape (N, B), dtype uint8, where entry (r, b) holds bit b of index r
        (b = 0 is the LSB).

    Complexity (paper Table tab:complexity):
        Per-pair:   O(B) bit-plane work
        Batch of N: O(N*B / L) wall-clock, O(N*B) space, L = SIMD lanes
        Bit-width:  unbounded (no 64-bit ceiling)

    Input / output convention:
        multiply_indices(i_arr, j_arr) -> (signs_array, indices_list)
        signs_array[k] in {-1, +1}, indices_list[k] = i_arr[k] ^ j_arr[k]
    """

    # ---------- representation: big-int <-> bit-plane ----------
    @staticmethod
    def ints_to_bitplane(vals, B):
        """
        Convert a sequence of Python big-ints to an (N, B) uint8 bit-plane.
        Entry (n, b) = bit b of vals[n]; b = 0 is the LSB.
        """
        N = len(vals)
        out = np.zeros((N, B), dtype=np.uint8)
        for n, v in enumerate(vals):
            v = int(v)
            if v == 0:
                continue
            nbytes = (min(B, v.bit_length()) + 7) // 8
            bits = np.unpackbits(
                np.frombuffer(v.to_bytes(nbytes, 'little'), dtype=np.uint8),
                bitorder='little')
            m = min(B, len(bits))
            out[n, :m] = bits[:m]
        return out

    @staticmethod
    def bitplane_to_ints(bits):
        """
        Convert an (N, B) bit-plane back to a list of Python big-ints.
        """
        N, B = bits.shape
        out = []
        nbytes = (B + 7) // 8
        for n in range(N):
            packed = np.packbits(bits[n], bitorder='little')
            out.append(int.from_bytes(packed.tobytes()[:nbytes], 'little'))
        return out

    @staticmethod
    def _required_width(i_arr, j_arr):
        """Minimal bit-width B that holds every index in both batches."""
        B = 1
        for v in i_arr:
            bl = int(v).bit_length()
            if bl > B:
                B = bl
        for v in j_arr:
            bl = int(v).bit_length()
            if bl > B:
                B = bl
        return B

    # ---------- core O(1) evaluator over bit-planes ----------
    def multiply_indices_bitplane(self, i_bits, j_bits):
        """
        Vectorized arbitrary-width O(1) multiplication on bit-planes.

        Args:
            i_bits: (N, B) uint8 bit-plane for the first indices
            j_bits: (N, B) uint8 bit-plane for the second indices

        Returns:
            (signs, idx_bits):
                signs   : (N,) int8 array, +1 / -1
                idx_bits: (N, B) uint8 bit-plane of the result indices i^j
        """
        i_bits = np.asarray(i_bits, dtype=np.uint8)
        j_bits = np.asarray(j_bits, dtype=np.uint8)
        N, B = i_bits.shape
        ar = np.arange(N)

        idx = i_bits ^ j_bits                       # result index = XOR

        # structural masks
        mask_zero = (~i_bits.any(axis=1)) | (~j_bits.any(axis=1))
        mask_diag = ((i_bits == j_bits).all(axis=1)) & ~mask_zero

        # nu2 via first-set-bit (argmax returns the index of the first 1)
        t  = np.argmax(idx,    axis=1)
        ti = np.argmax(i_bits, axis=1)
        tj = np.argmax(j_bits, axis=1)

        # branch selection
        cond1 = (t >= ti) & (t >= tj)              # diagonal break
        cond2 = (~cond1) & (ti >= tj)              # i-axis break
        cond3 = ~cond1 & ~cond2                    # j-axis break

        # branch signs
        s1 = 2 * i_bits[ar, t].astype(np.int64) - 1
        # i-axis: delta = -1 always here (Corollary cor:iaxis-collapse);
        # kept in general form to mirror the scalar evaluator exactly.
        delta = np.where(tj == ti, 1, -1)
        s2 = np.where(j_bits[ar, ti] == 0, delta, -delta)
        s3 = 1 - 2 * i_bits[ar, tj].astype(np.int64)
        s = np.where(cond1, s1, np.where(cond2, s2, s3))

        # higher-bit parity via reverse cumulative sum
        k = np.where(cond1, t + 1, np.where(cond2, ti + 1, tj + 1))
        or_bits = i_bits | j_bits
        rev_cumsum = np.cumsum(or_bits[:, ::-1], axis=1)[:, ::-1]
        parity_count = np.zeros(N, dtype=np.int64)
        valid = k < B
        parity_count[valid] = rev_cumsum[ar[valid], k[valid]]
        s = np.where((parity_count & 1) == 1, -s, s)

        # structural overrides
        s = np.where(mask_zero, 1, s)
        s = np.where(mask_diag, -1, s)

        return s.astype(np.int8), idx

    # ---------- convenience: Python ints in, ints out ----------
    def multiply_indices(self, i_arr, j_arr):
        """
        Convenience wrapper: sequences of Python big-ints in, integer results out.
        Mirrors FastStandard.multiply_indices but vectorized and arbitrary-width.

        Args:
            i_arr: sequence of first indices (Python ints, any width)
            j_arr: sequence of second indices (Python ints, any width)

        Returns:
            (signs, indices):
                signs   : (N,) int8 array, +1 / -1
                indices : list of N Python ints, indices[k] = i_arr[k] ^ j_arr[k]
        """
        i_arr = list(i_arr)
        j_arr = list(j_arr)
        if len(i_arr) != len(j_arr):
            raise ValueError("i_arr and j_arr must have the same length")
        if len(i_arr) == 0:
            return np.zeros(0, dtype=np.int8), []

        B = self._required_width(i_arr, j_arr)
        i_bits = self.ints_to_bitplane(i_arr, B)
        j_bits = self.ints_to_bitplane(j_arr, B)
        signs, idx_bits = self.multiply_indices_bitplane(i_bits, j_bits)
        indices = self.bitplane_to_ints(idx_bits)
        return signs, indices