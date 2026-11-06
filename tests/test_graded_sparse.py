"""
ADVERSARIAL TEST SUITE for GradedSparseStandard.
Cross-validates against proven FastStandard at extreme bit-widths.

Categories:
  - Correctness:   Exact match against scalar ground truth
  - EdgeCases:     Zeros, diagonals, identity, max bounds
  - Representation: int <-> graded conversion fidelity
  - Performance:   Sparse-optimal scaling verification

Run:  pytest tests/test_graded_sparse.py -v
"""
import numpy as np
import pytest
import time
import random

from hypercomplex_engine.core.fast.fast_standard import FastStandard
from hypercomplex_engine.core.graded.graded_sparse_standard import GradedSparseStandard


@pytest.fixture
def fast():
    return FastStandard()

@pytest.fixture
def graded():
    return GradedSparseStandard()


# =============================================================================
# 1. CORRECTNESS: Exact match against proven scalar engine
# =============================================================================
class TestCorrectness:
    """Verify GradedSparse matches FastStandard bit-for-bit."""
    # add 100000 for extra test
    @pytest.mark.parametrize("bits", [64, 65, 130, 300, 1000, 4096, 10000])
    def test_random_pairs_vs_scalar(self, fast, graded, bits):
        """
        Random pairs at extreme widths. The scalar engine handles arbitrary
        precision natively; GradedSparse must match exactly.
        """
        random.seed(42 + bits)
        n_pairs = 2000
        
        for _ in range(n_pairs):
            i = random.getrandbits(bits)
            j = random.getrandbits(bits)
            
            # Ground truth from proven scalar engine
            gt_sign, gt_idx = fast.multiply_indices(i, j)
            
            # GradedSparse result
            g_sign, g_idx = graded.multiply_indices(i, j)
            
            assert g_sign == gt_sign, (
                f"Sign mismatch at B={bits}, i={i}, j={j}: "
                f"graded={g_sign}, scalar={gt_sign}"
            )
            assert g_idx == gt_idx, (
                f"Index mismatch at B={bits}, i={i}, j={j}: "
                f"graded={g_idx}, scalar={gt_idx}"
            )
    @pytest.mark.parametrize("bits", [4, 8, 10])
    def test_exhaustive_small_widths(self, fast, graded, bits):
        """Exhaustive check for small widths (computationally feasible)."""
        max_val = 1 << bits
        for i in range(max_val):
            for j in range(max_val):
                gt_sign, gt_idx = fast.multiply_indices(i, j)
                g_sign, g_idx = graded.multiply_indices(i, j)
                assert g_sign == gt_sign
                assert g_idx == gt_idx


# =============================================================================
# 2. EDGE CASES: Structural positions and boundaries
# =============================================================================
class TestEdgeCases:
    """Test structural rules and extreme values."""

    def test_identity_element(self, graded):
        """e_0 * e_j = e_j and e_i * e_0 = e_i."""
        for j in [0, 1, 5, 100, 2**64 - 1, 2**1000]:
            s1, idx1 = graded.multiply_indices(0, j)
            assert s1 == 1 and idx1 == j
            
            s2, idx2 = graded.multiply_indices(j, 0)
            assert s2 == 1 and idx2 == j

    def test_main_diagonal(self, graded):
        """e_i * e_i = -e_0 for i > 0."""
        for i in [1, 5, 100, 2**64 - 1, 2**1000 - 1]:
            sign, idx = graded.multiply_indices(i, i)
            assert sign == -1 and idx == 0, (
                f"Diagonal fail at i={i}: sign={sign}, idx={idx}"
            )

    def test_anticommutativity(self, graded):
        """e_i * e_j = -(e_j * e_i) for distinct positive indices."""
        random.seed(99)
        for bits in [64, 1000, 10000]:
            for _ in range(100):
                i = random.getrandbits(bits) | 1  # ensure non-zero
                j = random.getrandbits(bits) | 1
                if i == j:
                    j += 1  # ensure distinct
                
                s_ij, idx_ij = graded.multiply_indices(i, j)
                s_ji, idx_ji = graded.multiply_indices(j, i)
                
                assert idx_ij == idx_ji, "XOR must be symmetric"
                assert s_ij == -s_ji, (
                    f"Anticommutativity violated at B={bits}: "
                    f"e_{i}*e_{j}={s_ij}, e_{j}*e_{i}={s_ji}"
                )

    def test_powers_of_two(self, fast, graded):
        """Every power of 2 from 2^0 to 2^63."""
        powers = [1 << k for k in range(64)]
        for i in powers:
            for j in powers:
                gt_s, gt_idx = fast.multiply_indices(i, j)
                g_s, g_idx = graded.multiply_indices(i, j)
                assert g_s == gt_s
                assert g_idx == gt_idx


# =============================================================================
# 3. REPRESENTATION: int <-> graded conversion fidelity
# =============================================================================
class TestRepresentation:
    """Test the graded form conversions."""

    def test_int_to_graded_roundtrip(self, graded):
        """int -> graded -> int must be identity."""
        test_values = [
            0, 1, 2, 3, 5, 15, 255, 2**64 - 1,
            2**100, 2**1000, 2**10000
        ]
        for val in test_values:
            pos = graded.int_to_graded(val)
            recovered = graded.graded_to_int(pos)
            assert recovered == val, (
                f"Roundtrip fail for {val}: "
                f"graded={pos}, recovered={recovered}"
            )

    def test_graded_to_notation(self, graded):
        """Test the Wilmot graded notation output."""
        # e_0 = 1
        assert graded.graded_to_notation(()) == "1"
        
        # e_1 = o1 (bit 0 set)
        assert graded.graded_to_notation((0,)) == "o1"
        
        # e_3 = o1,2 (bits 0,1 set)
        assert graded.graded_to_notation((0, 1)) == "o1,2"
        
        # e_5 = o1,3 (bits 0,2 set)
        assert graded.graded_to_notation((0, 2)) == "o1,3"

        # e_37 = o1,3,6 (bits 0,2 set) 0010 0101
        assert graded.graded_to_notation((0, 2, 5)) == "o1,3,6"
        
        # e_1024 = o11 (bit 10 set)
        assert graded.graded_to_notation((10,)) == "o11"

        # e_3405 = o1,3,4,7,9,11,12 (bits 0,2 set) 0000 1101 0100 1101
        assert graded.graded_to_notation((0,2,3,6,8,10,11)) == "o1,3,4,7,9,11,12"

    def test_sparse_representation(self, graded):
        """Verify that sparse indices use O(p) memory, not O(B)."""
        # A million-bit index with only 3 set bits
        huge = (1 << 1000000) | (1 << 500000) | (1 << 100)
        pos = graded.int_to_graded(huge)
        
        # Should have exactly 3 positions
        assert len(pos) == 3, f"Expected 3 positions, got {len(pos)}"
        assert pos == (100, 500000, 1000000)
        
        # Roundtrip must work
        assert graded.graded_to_int(pos) == huge


# =============================================================================
# 4. PERFORMANCE: Sparse-optimal scaling
# =============================================================================
class TestPerformance:
    """Verify O(p) complexity and sparse efficiency."""

    def test_sparse_vs_dense_performance(self, graded):
        """
        Sparse indices (few set bits) should be faster than dense indices
        (many set bits) at the same bit-width.
        """
        bits = 10000
        
        # Dense: all bits set
        dense = (1 << bits) - 1
        
        # Sparse: only 10 bits set
        sparse = sum(1 << (i * 1000) for i in range(10))
        
        # Time dense
        start = time.perf_counter()
        for _ in range(100):
            graded.multiply_indices(dense, dense)
        time_dense = time.perf_counter() - start
        
        # Time sparse
        start = time.perf_counter()
        for _ in range(100):
            graded.multiply_indices(sparse, sparse)
        time_sparse = time.perf_counter() - start
        
        # Sparse should be significantly faster (at least 10x for 10 vs 10000 bits)
        assert time_sparse < time_dense / 5, (
            f"Sparse not faster: dense={time_dense:.3f}s, sparse={time_sparse:.3f}s"
        )

    def test_million_bit_index_performance(self, graded):
        """
        A million-bit index with few set bits should complete quickly.
        This proves the engine doesn't materialize the full bit-width.
        """
        bits = 1_000_000
        # Only 5 set bits in a million-bit index
        i = (1 << 100) | (1 << 50000) | (1 << 250000) | (1 << 500000) | (1 << 999999)
        j = (1 << 200) | (1 << 100000) | (1 << 300000) | (1 << 600000) | (1 << 800000)
        
        start = time.perf_counter()
        sign, idx = graded.multiply_indices(i, j)
        elapsed = time.perf_counter() - start
        
        # Should complete in well under 100ms (typically < 1ms for 5-bit sparse)
        assert elapsed < 0.1, (
            f"Million-bit sparse took {elapsed:.3f}s - not O(p)!"
        )
        assert sign in (-1, 1)
        assert isinstance(idx, int)


# =============================================================================
# 5. STRESS: Extreme cases and collision handling
# =============================================================================
class TestStress:
    """Stress test with extreme inputs."""

    def test_max_uint64_boundary(self, fast, graded):
        """Indices near 2^64 boundary."""
        max64 = 2**64 - 1
        test_pairs = [
            (max64, 1),
            (max64, max64),
            (max64, max64 - 1),
            (max64 ^ 1, max64),
        ]
        for i, j in test_pairs:
            gt_s, gt_idx = fast.multiply_indices(i, j)
            g_s, g_idx = graded.multiply_indices(i, j)
            assert g_s == gt_s
            assert g_idx == gt_idx

    def test_consecutive_indices(self, fast, graded):
        """Consecutive i, i+1 pairs stress nu2 boundaries."""
        for i in range(1000):
            j = i + 1
            gt_s, gt_idx = fast.multiply_indices(i, j)
            g_s, g_idx = graded.multiply_indices(i, j)
            assert g_s == gt_s
            assert g_idx == gt_idx

    def test_all_bits_set_pattern(self, fast, graded):
        """Numbers of form 2^n - 1 (all bits set)."""
        for bits in [8, 16, 32, 64, 128, 256]:
            val = (1 << bits) - 1
            gt_s, gt_idx = fast.multiply_indices(val, val)
            g_s, g_idx = graded.multiply_indices(val, val)
            assert g_s == gt_s
            assert g_idx == gt_idx


if __name__ == "__main__":
    pytest.main([__file__, "-v"])