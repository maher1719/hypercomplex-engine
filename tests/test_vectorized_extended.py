"""
EXTENDED ADVERSARIAL TEST SUITE for VectorizedStandard.
Designed to BREAK the engine. If these pass, the engine is bulletproof.

Categories:
  - AlgebraicProperties:   Mathematical invariants that MUST hold
  - ExhaustiveConsistency: Brute-force comparison vs proven scalar engine
  - BoundaryConditions:    Extreme values, overflow, edge indices
  - ArrayEdgeCases:        Empty, sparse, asymmetric, degenerate inputs
  - NumericalStability:    Float precision, cancellation, extreme magnitudes
  - Stress:                Memory, collisions, performance regression

Run quick:    pytest tests/test_vectorized_extended.py -v -m "not slow"
Run full:     pytest tests/test_vectorized_extended.py -v
"""
import numpy as np
import pytest
import time

from hypercomplex_engine.core.fast.fast_standard import FastStandard
from hypercomplex_engine.core.vectorized.vectorized_standard import VectorizedStandard


@pytest.fixture
def fast():
    return FastStandard()

@pytest.fixture
def vec():
    return VectorizedStandard()


# =============================================================================
# 1. ALGEBRAIC PROPERTIES — Mathematical invariants that MUST hold
# =============================================================================
class TestAlgebraicProperties:
    """These test the MATH, not just the implementation."""

    def test_anticommutativity(self, fast, vec):
        """e_i * e_j = -(e_j * e_i) for i != j, both non-zero."""
        np.random.seed(1)
        n = 5000
        i = np.random.randint(1, 2**20, size=n, dtype=np.uint64)
        j = np.random.randint(1, 2**20, size=n, dtype=np.uint64)
        # Ensure i != j
        mask = i == j
        j[mask] += 1

        s_ij, idx_ij = vec.multiply_indices_vec(i, j)
        s_ji, idx_ji = vec.multiply_indices_vec(j, i)

        # Indices must match (XOR is commutative)
        np.testing.assert_array_equal(idx_ij, idx_ji,
            err_msg="XOR index must be symmetric")
        # Signs must be opposite
        np.testing.assert_array_equal(s_ij, -s_ji,
            err_msg="Anticommutativity violated: e_i*e_j != -(e_j*e_i)")

    def test_square_is_negative_identity(self, vec):
        """e_i * e_i = -e_0 for all i > 0 (sign=-1, index=0)."""
        i = np.arange(1, 100_000, dtype=np.uint64)
        signs, idxs = vec.multiply_indices_vec(i, i)
        np.testing.assert_array_equal(signs, np.full(len(i), -1, dtype=np.int8),
            err_msg="e_i^2 must have sign -1")
        np.testing.assert_array_equal(idxs, np.zeros(len(i), dtype=np.uint64),
            err_msg="e_i^2 must have index 0")

    def test_identity_element(self, vec):
        """e_0 * e_j = e_j and e_i * e_0 = e_i."""
        j = np.arange(0, 50_000, dtype=np.uint64)
        zeros = np.zeros(len(j), dtype=np.uint64)

        # 0 * j
        s1, idx1 = vec.multiply_indices_vec(zeros, j)
        np.testing.assert_array_equal(s1, np.ones(len(j), dtype=np.int8))
        np.testing.assert_array_equal(idx1, j)

        # j * 0
        s2, idx2 = vec.multiply_indices_vec(j, zeros)
        np.testing.assert_array_equal(s2, np.ones(len(j), dtype=np.int8))
        np.testing.assert_array_equal(idx2, j)

    def test_alternativity_octonions(self, vec):
        """
        Octonions (dim<=8) are alternative:
        (e_i * e_i) * e_j == e_i * (e_i * e_j)
        """
        # Test all pairs in octonion range
        for i in range(1, 8):
            for j in range(1, 8):
                # (e_i * e_i) * e_j
                s1, k = vec.multiply_indices_vec(
                    np.array([i], dtype=np.uint64),
                    np.array([i], dtype=np.uint64))
                s2, m_left = vec.multiply_indices_vec(
                    np.array([k[0]], dtype=np.uint64),
                    np.array([j], dtype=np.uint64))
                left = (int(s1[0]) * int(s2[0]), int(m_left[0]))

                # e_i * (e_i * e_j)
                s3, k2 = vec.multiply_indices_vec(
                    np.array([i], dtype=np.uint64),
                    np.array([j], dtype=np.uint64))
                s4, m_right = vec.multiply_indices_vec(
                    np.array([i], dtype=np.uint64),
                    np.array([k2[0]], dtype=np.uint64))
                right = (int(s3[0]) * int(s4[0]), int(m_right[0]))

                assert left == right, (
                    f"Alternativity violated at i={i}, j={j}: {left} != {right}")

    def test_power_associativity(self, vec):
        """(e_i^2)^2 == e_i^4 — powers of a single element associate."""
        i = np.array([7, 15, 31, 63], dtype=np.uint64)
        for idx in i:
            # e_i^2
            s1, k1 = vec.multiply_indices_vec(
                np.array([idx]), np.array([idx]))
            # (e_i^2)^2
            s2, k2 = vec.multiply_indices_vec(
                np.array([k1[0]]), np.array([k1[0]]))
            # Should be (+1, 0) since (-1)^2 = +1
            assert int(k2[0]) == 0, f"e_{idx}^4 index != 0"
            assert int(s1[0])**2 * int(s2[0]) == 1, f"e_{idx}^4 sign != +1"


# =============================================================================
# 2. EXHAUSTIVE CONSISTENCY — Brute force vs proven scalar engine
# =============================================================================
class TestExhaustiveConsistency:
    """Compare every single output against the proven FastStandard."""

    @pytest.mark.parametrize("dim_level", [1, 2, 3, 4, 5, 6])
    def test_exhaustive_vs_scalar(self, fast, vec, dim_level):
        """Every pair in A_n matches scalar FastStandard exactly."""
        dim = 1 << dim_level
        i_all, j_all = np.meshgrid(
            np.arange(dim, dtype=np.uint64),
            np.arange(dim, dtype=np.uint64),
            indexing='ij')
        i_flat = i_all.ravel()
        j_flat = j_all.ravel()

        # Scalar ground truth
        gt_signs = np.zeros(len(i_flat), dtype=np.int8)
        gt_idxs = np.zeros(len(i_flat), dtype=np.uint64)
        for k in range(len(i_flat)):
            s, idx = fast.multiply_indices(int(i_flat[k]), int(j_flat[k]))
            gt_signs[k] = s
            gt_idxs[k] = idx

        # Vectorized
        v_signs, v_idxs = vec.multiply_indices_vec(i_flat, j_flat)

        np.testing.assert_array_equal(v_signs, gt_signs,
            err_msg=f"Sign mismatch in dim 2^{dim_level}")
        np.testing.assert_array_equal(v_idxs, gt_idxs,
            err_msg=f"Index mismatch in dim 2^{dim_level}")

    def test_random_60bit_vs_scalar(self, fast, vec):
        """200k random 60-bit pairs match scalar exactly."""
        np.random.seed(777)
        n = 200_000
        i = np.random.randint(0, 2**60, size=n, dtype=np.uint64)
        j = np.random.randint(0, 2**60, size=n, dtype=np.uint64)

        gt_signs = np.zeros(n, dtype=np.int8)
        gt_idxs = np.zeros(n, dtype=np.uint64)
        for k in range(n):
            s, idx = fast.multiply_indices(int(i[k]), int(j[k]))
            gt_signs[k] = s
            gt_idxs[k] = idx

        v_signs, v_idxs = vec.multiply_indices_vec(i, j)
        np.testing.assert_array_equal(v_signs, gt_signs)
        np.testing.assert_array_equal(v_idxs, gt_idxs)

    @pytest.mark.slow
    def test_exhaustive_dim10(self, fast, vec):
        """Full exhaustive check for dim 2^10 = 1024 (1M pairs)."""
        dim = 1024
        i_all, j_all = np.meshgrid(
            np.arange(dim, dtype=np.uint64),
            np.arange(dim, dtype=np.uint64),
            indexing='ij')
        i_flat = i_all.ravel()
        j_flat = j_all.ravel()

        v_signs, v_idxs = vec.multiply_indices_vec(i_flat, j_flat)

        # Spot-check 10000 random positions against scalar
        np.random.seed(42)
        sample = np.random.choice(len(i_flat), size=10000, replace=False)
        for k in sample:
            s, idx = fast.multiply_indices(int(i_flat[k]), int(j_flat[k]))
            assert int(v_signs[k]) == s, f"Sign mismatch at {i_flat[k]},{j_flat[k]}"
            assert int(v_idxs[k]) == idx, f"Idx mismatch at {i_flat[k]},{j_flat[k]}"


# =============================================================================
# 3. BOUNDARY CONDITIONS — Extreme values designed to trigger overflow/underflow
# =============================================================================
class TestBoundaryConditions:

    def test_max_uint64(self, fast, vec):
        """Indices near 2^64 - 1 must not overflow."""
        max_val = np.uint64(2**64 - 1)
        i = np.array([max_val, max_val - 1, max_val ^ 1], dtype=np.uint64)
        j = np.array([1, max_val, 2], dtype=np.uint64)

        v_signs, v_idxs = vec.multiply_indices_vec(i, j)
        for k in range(len(i)):
            s, idx = fast.multiply_indices(int(i[k]), int(j[k]))
            assert int(v_signs[k]) == s, f"Sign fail at max_uint64 case {k}"
            assert int(v_idxs[k]) == idx, f"Idx fail at max_uint64 case {k}"

    def test_powers_of_two(self, fast, vec):
        """Every power of 2 from 2^0 to 2^63."""
        powers = np.array([1 << k for k in range(64)], dtype=np.uint64)
        i, j = np.meshgrid(powers, powers, indexing='ij')
        i_flat, j_flat = i.ravel(), j.ravel()

        v_signs, v_idxs = vec.multiply_indices_vec(i_flat, j_flat)
        for k in range(len(i_flat)):
            s, idx = fast.multiply_indices(int(i_flat[k]), int(j_flat[k]))
            assert int(v_signs[k]) == s
            assert int(v_idxs[k]) == idx

    def test_all_bits_set(self, fast, vec):
        """Numbers of form 2^n - 1 (all bits set)."""
        vals = np.array([(1 << k) - 1 for k in range(1, 64)], dtype=np.uint64)
        i, j = np.meshgrid(vals, vals, indexing='ij')
        i_flat, j_flat = i.ravel(), j.ravel()

        v_signs, v_idxs = vec.multiply_indices_vec(i_flat, j_flat)
        for k in range(0, len(i_flat), 37):  # sample to keep fast
            s, idx = fast.multiply_indices(int(i_flat[k]), int(j_flat[k]))
            assert int(v_signs[k]) == s
            assert int(v_idxs[k]) == idx

    def test_consecutive_indices(self, fast, vec):
        """Consecutive i, i+1 pairs — stresses nu2 boundary."""
        i = np.arange(0, 100_000, dtype=np.uint64)
        j = i + 1
        v_signs, v_idxs = vec.multiply_indices_vec(i, j)
        for k in range(0, len(i), 997):  # prime stride sampling
            s, idx = fast.multiply_indices(int(i[k]), int(j[k]))
            assert int(v_signs[k]) == s
            assert int(v_idxs[k]) == idx


# =============================================================================
# 4. ARRAY EDGE CASES — Degenerate, empty, sparse, asymmetric inputs
# =============================================================================
class TestArrayEdgeCases:

    def test_single_element_arrays(self, vec):
        """1-D real multiplication."""
        r = vec.multiply_arrays(np.array([3.0]), np.array([4.0]))
        assert len(r) == 1
        assert np.isclose(r[0], 12.0)

    def test_all_zero_arrays(self, vec):
        """Zero vector times anything is zero."""
        dim = 64
        z = np.zeros(dim)
        x = np.random.randn(dim)
        np.testing.assert_array_equal(vec.multiply_arrays(z, x), np.zeros(dim))
        np.testing.assert_array_equal(vec.multiply_arrays(x, z), np.zeros(dim))
        np.testing.assert_array_equal(vec.multiply_arrays(z, z), np.zeros(dim))

    def test_extreme_sparsity(self, vec):
        """1 non-zero in a 2^16 vector — must not OOM, must be correct."""
        dim = 1 << 16
        a = np.zeros(dim); a[5] = 2.0
        b = np.zeros(dim); b[7] = 3.0
        result = vec.multiply_arrays(a, b)
        # e_5 * e_7 -> single contribution
        fast = FastStandard()
        s, idx = fast.multiply_indices(5, 7)
        expected = np.zeros(dim)
        expected[idx] = s * 6.0
        np.testing.assert_allclose(result, expected)

    def test_asymmetric_large_small(self, vec):
        """Tiny array times huge array."""
        small = np.array([1.0, 2.0])           # dim 2
        big = np.random.randn(1024)             # dim 1024
        result = vec.multiply_arrays(small, big)
        assert len(result) == 1024  # padded to next power of 2 >= 1024

    def test_non_power_of_two_sizes(self, vec):
        """Sizes 3, 5, 7 must pad correctly."""
        for size in [3, 5, 7, 9, 15, 17]:
            a = np.random.randn(size)
            b = np.random.randn(size)
            result = vec.multiply_arrays(a, b)
            expected_dim = 1 << (size - 1).bit_length()
            assert len(result) == expected_dim, f"Failed padding for size {size}"

    def test_negative_coefficients(self, fast, vec):
        """Negative coefficients must propagate signs correctly."""
        dim = 16
        a = np.random.randn(dim)
        b = np.random.randn(dim)
        result = vec.multiply_arrays(a, b)

        # Brute force with negatives
        expected = np.zeros(dim)
        for i in range(dim):
            for j in range(dim):
                s, idx = fast.multiply_indices(i, j)
                expected[idx] += s * a[i] * b[j]
        np.testing.assert_allclose(result, expected, rtol=1e-12, atol=1e-12)


# =============================================================================
# 5. NUMERICAL STABILITY — Float precision torture
# =============================================================================
class TestNumericalStability:

    def test_large_coefficients(self, vec):
        """Coefficients near float64 max don't produce NaN unexpectedly."""
        dim = 8
        a = np.array([1e150, 0, 0, 0, 0, 0, 0, 0])
        b = np.array([1e150, 0, 0, 0, 0, 0, 0, 0])
        result = vec.multiply_arrays(a, b)
        # e_0 * e_0 = e_0, coefficient 1e300 (within float64 range)
        assert np.isfinite(result[0])
        assert np.isclose(result[0], 1e300, rtol=1e-10)

    def test_small_coefficients(self, vec):
        """Tiny coefficients don't underflow to zero incorrectly."""
        dim = 8
        a = np.array([1e-150, 0, 0, 0, 0, 0, 0, 0])
        b = np.array([1e-150, 0, 0, 0, 0, 0, 0, 0])
        result = vec.multiply_arrays(a, b)
        assert np.isclose(result[0], 1e-300, rtol=1e-10)

    def test_catastrophic_cancellation(self, fast, vec):
        """Values that nearly cancel must still be accurate."""
        dim = 8
        # Craft inputs where multiple terms hit same index with opposite signs
        a = np.zeros(dim); a[1] = 1.0; a[2] = 1.0
        b = np.zeros(dim); b[1] = 1.0; b[2] = -1.0
        result = vec.multiply_arrays(a, b)

        expected = np.zeros(dim)
        for i in range(dim):
            for j in range(dim):
                s, idx = fast.multiply_indices(i, j)
                expected[idx] += s * a[i] * b[j]
        np.testing.assert_allclose(result, expected, rtol=1e-12, atol=1e-15)

    def test_accumulation_many_terms(self, fast, vec):
        """Many terms accumulating to same index stay precise."""
        dim = 32
        np.random.seed(5)
        a = np.random.randn(dim)
        b = np.random.randn(dim)
        result = vec.multiply_arrays(a, b)

        expected = np.zeros(dim)
        for i in range(dim):
            for j in range(dim):
                s, idx = fast.multiply_indices(i, j)
                expected[idx] += s * a[i] * b[j]
        np.testing.assert_allclose(result, expected, rtol=1e-10, atol=1e-12)


# =============================================================================
# 6. STRESS — Collisions, memory, performance
# =============================================================================
class TestStress:

    def test_massive_xor_collisions(self, fast, vec):
        """
        Force MANY pairs to collide on the same XOR index.
        All pairs (i, i^target) map to 'target'.
        """
        dim = 64
        target = 5
        a = np.random.randn(dim)
        b = np.random.randn(dim)
        result = vec.multiply_arrays(a, b)

        # Manually sum all pairs hitting index 'target'
        expected_target = 0.0
        for i in range(dim):
            j = i ^ target
            if j < dim:
                s, idx = fast.multiply_indices(i, j)
                assert idx == target
                expected_target += s * a[i] * b[j]
        assert np.isclose(result[target], expected_target, rtol=1e-10)

    @pytest.mark.slow
    def test_memory_sparse_huge_dim(self, vec):
        """Sparse multiplication at dim 2^20 must not OOM."""
        dim = 1 << 20
        a = np.zeros(dim); a[3] = 1.5; a[1000] = -2.0
        b = np.zeros(dim); b[7] = 0.5; b[2000] = 3.0
        start = time.perf_counter()
        result = vec.multiply_arrays(a, b)
        elapsed = time.perf_counter() - start
        assert len(result) == dim
        # Should be nearly instant due to sparsity
        assert elapsed < 1.0, f"Sparse 2^20 took {elapsed:.3f}s"

    def test_performance_not_python_loop(self, vec):
        """Dense 512x512 must be far faster than a Python loop would be."""
        dim = 512
        a = np.random.randn(dim)
        b = np.random.randn(dim)
        start = time.perf_counter()
        vec.multiply_arrays(a, b)
        elapsed = time.perf_counter() - start
        # A Python double-loop over 512^2=262k pairs takes >0.5s.
        # SIMD should be well under that.
        assert elapsed < 0.4, f"Took {elapsed:.3f}s — possible Python-loop regression"