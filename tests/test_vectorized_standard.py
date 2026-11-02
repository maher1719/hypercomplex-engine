import numpy as np
import pytest
import time

from hypercomplex_engine.core.fast.fast_standard import FastStandard
from hypercomplex_engine.core.vectorized.vectorized_standard import VectorizedStandard

@pytest.fixture
def engines():
    return FastStandard(), VectorizedStandard()

# =============================================================================
# 1. THE MILLION-PAIR GAUNTLET (Indices Accuracy)
# =============================================================================
def test_multiply_indices_vec_randomized(engines):
    """
    Fires 100,000 random 60-bit indices at the SIMD engine and compares 
    every single sign and XOR index against the scalar FastStandard.
    """
    fast, vec = engines
    np.random.seed(42)
    
    i_arr = np.random.randint(0, 2**60, size=100_000, dtype=np.uint64)
    j_arr = np.random.randint(0, 2**60, size=100_000, dtype=np.uint64)

    gt_signs = np.zeros(100_000, dtype=np.int8)
    gt_idxs = np.zeros(100_000, dtype=np.uint64)
    for idx in range(100_000):
        s, i = fast.multiply_indices(int(i_arr[idx]), int(j_arr[idx]))
        gt_signs[idx] = s
        gt_idxs[idx] = i

    v_signs, v_idxs = vec.multiply_indices_vec(i_arr, j_arr)

    np.testing.assert_array_equal(v_signs, gt_signs, err_msg="SIMD sign mismatch!")
    np.testing.assert_array_equal(v_idxs, gt_idxs, err_msg="SIMD XOR index mismatch!")

# =============================================================================
# 2. THE EDGE-CASE MINEFIELD
# =============================================================================
def test_multiply_indices_vec_edge_cases(engines):
    fast, vec = engines
    i_arr = np.array([0, 1, 5, 10, 0, 7, 2**63 - 1], dtype=np.uint64)
    j_arr = np.array([0, 1, 0, 10, 5, 7, 2**63 - 1], dtype=np.uint64)
    
    gt_signs, gt_idxs = [], []
    for i, j in zip(i_arr, j_arr):
        s, idx = fast.multiply_indices(int(i), int(j))
        gt_signs.append(s)
        gt_idxs.append(idx)

    v_signs, v_idxs = vec.multiply_indices_vec(i_arr, j_arr)
    np.testing.assert_array_equal(v_signs, np.array(gt_signs, dtype=np.int8))
    np.testing.assert_array_equal(v_idxs, np.array(gt_idxs, dtype=np.uint64))

# =============================================================================
# 3. POLYNOMIAL EQUIVALENCE (Dense Array Multiplication)
# =============================================================================
def test_multiply_arrays_dense_equivalence(engines):
    fast, vec = engines
    dim = 64
    np.random.seed(99)
    arr1 = np.random.randn(dim)
    arr2 = np.random.randn(dim)

    res_naive = np.zeros(dim, dtype=np.float64)
    for i in range(dim):
        for j in range(dim):
            s, idx = fast.multiply_indices(i, j)
            res_naive[idx] += s * arr1[i] * arr2[j]

    res_vec = vec.multiply_arrays(arr1, arr2)
    np.testing.assert_allclose(res_vec, res_naive, rtol=1e-12, atol=1e-12)

# =============================================================================
# 4. THE "XOR COLLISION" TRAP (CORRECTED FOR NON-COMMUTATIVITY)
# =============================================================================
def test_multiply_arrays_xor_collisions(engines):
    """
    Ensures np.add.at ACCUMULATES correctly for ALL pairs that map to the same XOR index.
    """
    fast, vec = engines
    arr1 = np.zeros(8, dtype=np.float64)
    arr2 = np.zeros(8, dtype=np.float64)
    
    arr1[1] = 2.0; arr1[2] = 3.0; arr1[0] = 1.5; arr1[3] = -1.0
    arr2[2] = 4.0; arr2[1] = 5.0; arr2[3] = 2.0; arr2[0] = -2.0

    res_vec = vec.multiply_arrays(arr1, arr2)
    
    # Calculate the EXACT expected sum for index 3 by checking ALL 8x8 pairs
    expected_idx3 = 0.0
    for i in range(8):
        for j in range(8):
            s, idx = fast.multiply_indices(i, j)
            if idx == 3:
                expected_idx3 += s * arr1[i] * arr2[j]
                
    # The engine correctly calculated -2.0 (because e2*e1 = -e3, contributing -15.0)
    assert np.isclose(res_vec[3], expected_idx3), f"Expected {expected_idx3}, got {res_vec[3]}"

# =============================================================================
# 5. ASYMMETRIC PADDING & SPARSITY (CORRECTED FOR OVERLAPS)
# =============================================================================
def test_multiply_arrays_asymmetric_padding(engines):
    fast, vec = engines
    arr1 = np.array([1.0, 2.0, 3.0]) 
    arr2 = np.array([4.0, 5.0, 6.0, 7.0, 8.0]) 
    
    res_vec = vec.multiply_arrays(arr1, arr2)
    assert len(res_vec) == 8, "Failed to pad asymmetric arrays to power-of-2!"

    # Calculate exact expected sum for index 1 (0^1, 1^0, 2^3, 3^2, etc.)
    expected_idx1 = 0.0
    for i in range(8):
        for j in range(8):
            s, idx = fast.multiply_indices(i, j)
            if idx == 1:
                # Use padded arrays for calculation
                v1 = arr1[i] if i < len(arr1) else 0.0
                v2 = arr2[j] if j < len(arr2) else 0.0
                expected_idx1 += s * v1 * v2
                
    # The engine correctly calculated 34.0 (because e2*e3 = +e1, contributing +21.0)
    assert np.isclose(res_vec[1], expected_idx1), f"Expected {expected_idx1}, got {res_vec[1]}"

# =============================================================================
# 6. THE PERFORMANCE GATE (RELATIVE SPEEDUP)
# =============================================================================
def test_multiply_arrays_performance_gate(engines):
    """
    Compares SIMD against a small scalar baseline to prove it's not a Python loop.
    """
    fast, vec = engines
    dim = 512  # 512x512 = 262,144 pairs
    arr1 = np.random.randn(dim)
    arr2 = np.random.randn(dim)

    # Time SIMD
    start = time.perf_counter()
    vec.multiply_arrays(arr1, arr2)
    time_simd = time.perf_counter() - start

    # Time a small scalar batch (10,000 pairs) to establish baseline speed
    start = time.perf_counter()
    for i in range(100):
        for j in range(100):
            fast.multiply_indices(i, j)
    time_scalar_10k = time.perf_counter() - start
    
    # Extrapolate scalar time to 262,144 pairs
    time_scalar_extrapolated = time_scalar_10k * (262144 / 10000)

    # SIMD MUST be at least 2x faster than the extrapolated scalar loop
    assert time_simd < (time_scalar_extrapolated / 2.0), \
        f"SIMD took {time_simd:.3f}s, expected < {time_scalar_extrapolated/2.0:.3f}s (half of scalar time)"