"""
ADVERSARIAL TEST SUITE for BitPlaneStandard.
Cross-validates against proven FastStandard at extreme bit-widths,
proving the 64-bit SIMD ceiling is broken.

Run:  pytest tests/test_bitplane_standard.py -v
"""
import numpy as np
import pytest
import time
import random

from hypercomplex_engine.core.fast.fast_standard import FastStandard
from hypercomplex_engine.core.bitplane.bitplane_standard import BitPlaneStandard


@pytest.fixture
def fast():
    return FastStandard()

@pytest.fixture
def bp():
    return BitPlaneStandard()


# =============================================================================
# 1. CORRECTNESS vs proven scalar engine
# =============================================================================
class TestCorrectness:
    
    @pytest.mark.parametrize("bits", [8, 64, 65, 130, 300, 1000, 4096, 10000, 100000])
    def test_random_pairs_vs_scalar(self, fast, bp, bits):
        """Random batches at extreme widths must match FastStandard exactly."""
        random.seed(42 + bits)
        n_pairs = 1000
        i_arr = [random.getrandbits(bits) for _ in range(n_pairs)]
        j_arr = [random.getrandbits(bits) for _ in range(n_pairs)]

        signs, indices = bp.multiply_indices(i_arr, j_arr)

        for k in range(n_pairs):
            gt_sign, gt_idx = fast.multiply_indices(i_arr[k], j_arr[k])
            assert int(signs[k]) == gt_sign, f"sign mismatch B={bits} pair {k}"
            assert indices[k] == gt_idx,     f"index mismatch B={bits} pair {k}"
    # you could add more
    @pytest.mark.parametrize("bits", [4, 8, 10])
    def test_exhaustive_small_widths(self, fast, bp, bits):
        """Exhaustive check for small widths."""
        max_val = 1 << bits
        i_arr = [i for i in range(max_val) for j in range(max_val)]
        j_arr = [j for i in range(max_val) for j in range(max_val)]
        signs, indices = bp.multiply_indices(i_arr, j_arr)
        for k in range(len(i_arr)):
            gt_sign, gt_idx = fast.multiply_indices(i_arr[k], j_arr[k])
            assert int(signs[k]) == gt_sign
            assert indices[k] == gt_idx


# =============================================================================
# 2. EDGE CASES
# =============================================================================
class TestEdgeCases:

    def test_identity_and_diagonal(self, bp):
        """e_0*e_j = e_j ; e_i*e_i = -e_0."""
        signs, indices = bp.multiply_indices([0, 5, 7], [9, 5, 7])
        assert int(signs[0]) == 1 and indices[0] == 9      # 0 * 9
        assert int(signs[1]) == -1 and indices[1] == 0     # 5 * 5
        assert int(signs[2]) == -1 and indices[2] == 0     # 7 * 7

    def test_anticommutativity(self, bp):
        """e_i*e_j = -(e_j*e_i) for distinct positive indices."""
        random.seed(7)
        i_arr = [random.getrandbits(100) | 1 for _ in range(50)]
        j_arr = [random.getrandbits(100) | 1 for _ in range(50)]
        # ensure distinct
        j_arr = [j if j != i else j + 1 for i, j in zip(i_arr, j_arr)]

        s_ij, idx_ij = bp.multiply_indices(i_arr, j_arr)
        s_ji, idx_ji = bp.multiply_indices(j_arr, i_arr)

        for k in range(len(i_arr)):
            assert idx_ij[k] == idx_ji[k], "XOR must be symmetric"
            assert int(s_ij[k]) == -int(s_ji[k]), "anticommutativity violated"

    def test_empty_batch(self, bp):
        signs, indices = bp.multiply_indices([], [])
        assert len(signs) == 0 and len(indices) == 0

    def test_length_mismatch_raises(self, bp):
        with pytest.raises(ValueError):
            bp.multiply_indices([1, 2], [3])


# =============================================================================
# 3. REPRESENTATION roundtrip
# =============================================================================
class TestRepresentation:

    def test_bitplane_roundtrip(self, bp):
        vals = [0, 1, 5, 255, 2**64 - 1, 2**200, 2**1000 + 7]
        B = max(v.bit_length() for v in vals) or 1
        plane = bp.ints_to_bitplane(vals, B)
        recovered = bp.bitplane_to_ints(plane)
        assert recovered == vals


# =============================================================================
# 4. THE 64-BIT CEILING IS BROKEN
# =============================================================================
class TestCeilingBroken:

    def test_beyond_uint64(self, fast, bp):
        """Indices wider than 64 bits: NumPy uint64 cannot hold these,
        but BitPlaneStandard must still match the scalar engine."""
        random.seed(64)
        bits = 200   # > 64: impossible for np.uint64
        i_arr = [random.getrandbits(bits) for _ in range(500)]
        j_arr = [random.getrandbits(bits) for _ in range(500)]
        signs, indices = bp.multiply_indices(i_arr, j_arr)
        for k in range(500):
            gt_s, gt_i = fast.multiply_indices(i_arr[k], j_arr[k])
            assert int(signs[k]) == gt_s
            assert indices[k] == gt_i


if __name__ == "__main__":
    pytest.main([__file__, "-v"])