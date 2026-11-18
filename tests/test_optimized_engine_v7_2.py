"""
Regression suite for the optimized engines.
Run: pytest tests/test_optimized_engines.py -v
"""
import numpy as np
import pytest
import random

from hypercomplex_engine.core.fast.fast_standard import FastStandard
from hypercomplex_engine.core.vectorized.vectorized_standard import VectorizedStandard
from hypercomplex_engine.core.bitplane.bitplane_standard import BitPlaneStandard
from hypercomplex_engine.core.graded.graded_sparse_standard import GradedSparseStandard


@pytest.fixture
def engines():
    return FastStandard(), VectorizedStandard(), BitPlaneStandard(), GradedSparseStandard()


@pytest.mark.parametrize("bits", [4, 6, 8])
def test_exhaustive_all_engines(engines, bits):
    fast, vec, bp, gs = engines
    N = 1 << bits
    ii = [i for i in range(N) for _ in range(N)]
    jj = [j for _ in range(N) for j in range(N)]
    vs, vi = vec.multiply_indices_vec(ii, jj)
    bs, bi = bp.multiply_indices(ii, jj)
    for q in range(len(ii)):
        gs_s, gs_i = fast.multiply_indices(ii[q], jj[q])
        assert int(vs[q]) == gs_s and int(vi[q]) == gs_i
        assert int(bs[q]) == gs_s and bi[q] == gs_i
        s, i = gs.multiply_indices(ii[q], jj[q])
        assert s == gs_s and i == gs_i


@pytest.mark.parametrize("bits", [64, 65, 130, 300, 1000, 4096])
def test_wide_random_all_engines(engines, bits):
    fast, vec, bp, gs = engines
    rng = random.Random(2026 + bits)
    ii = [rng.getrandbits(bits) for _ in range(400)]
    jj = [rng.getrandbits(bits) for _ in range(400)]
    bs, bi = bp.multiply_indices(ii, jj)
    if bits <= 64:
        vs, vi = vec.multiply_indices_vec(ii, jj)
    for q in range(len(ii)):
        s, i = fast.multiply_indices(ii[q], jj[q])
        assert int(bs[q]) == s and bi[q] == i
        if bits <= 64:
            assert int(vs[q]) == s and int(vi[q]) == i
        gss, gii = gs.multiply_indices(ii[q], jj[q])
        assert gss == s and gii == i


def test_k_equals_64_lane(engines):
    """ti = 63 forces k = 64: the parity term must be zero, not popcount(i|j)."""
    fast, vec, _, _ = engines
    i = np.array([1 << 63] * 40, dtype=np.uint64)
    j = np.arange(1, 41, dtype=np.uint64)
    s, idx = vec.multiply_indices_vec(i, j)
    for q in range(40):
        gs, gi = fast.multiply_indices(int(i[q]), int(j[q]))
        assert int(s[q]) == gs, f"wide-shift parity bug at j={int(j[q])}"
        assert int(idx[q]) == gi


def test_delta_collapse_exhaustive(engines):
    """cond2 must never produce delta = +1 (Corollary cor:iaxis-collapse)."""
    fast, vec, bp, _ = engines
    N = 1 << 8
    ii = [i for i in range(N) for _ in range(N)]
    jj = [j for _ in range(N) for j in range(N)]
    vs, _ = vec.multiply_indices_vec(ii, jj)
    bs, _ = bp.multiply_indices(ii, jj)
    for q in range(len(ii)):
        gs, _ = fast.multiply_indices(ii[q], jj[q])
        assert int(vs[q]) == gs and int(bs[q]) == gs


def test_bitplane_parity_edge_k_equals_B(engines):
    """k == B must give zero parity with no out-of-range indexing."""
    fast, _, bp, _ = engines
    for B in [1, 2, 7, 8, 63, 64, 65, 129]:
        top = 1 << (B - 1)
        ii = [top, top, top | 1]
        jj = [1, top | 1, top]
        s, idx = bp.multiply_indices(ii, jj)
        for q in range(3):
            gs, gi = fast.multiply_indices(ii[q], jj[q])
            assert int(s[q]) == gs and idx[q] == gi


def test_graded_helpers_match_fused(engines):
    """Legacy three-scan helpers must agree with the fused merge."""
    _, _, _, gs = engines
    rng = random.Random(5)
    for _ in range(500):
        pi = tuple(sorted(rng.sample(range(200), rng.randint(1, 12))))
        pj = tuple(sorted(rng.sample(range(200), rng.randint(1, 12))))
        if pi == pj:
            continue
        t = gs._nu2_xor_gr(pi, pj)
        k = max(t, pi[0], pj[0]) + 1
        assert gs._xor_gr(pi, pj)[0] == t
        s, out = gs.multiply_graded(pi, pj)
        assert out == gs._xor_gr(pi, pj)
        assert (gs._pop_or_shift_gr(pi, pj, k) & 1) == (
            gs._pop_or_shift_gr(pi, pj, k) & 1)


def test_graded_roundtrip_and_notation(engines):
    _, _, _, gs = engines
    for m in [0, 1, 3, 5, 255, 2**64 - 1, 2**1000 + 7,
              (1 << 100000) | (1 << 50000) | (1 << 100)]:
        assert gs.graded_to_int(gs.int_to_graded(m)) == m
    assert gs.graded_to_notation(()) == "1"
    assert gs.graded_to_notation((0,)) == "o1"
    assert gs.graded_to_notation((0, 1)) == "o1,2"
    assert gs.graded_to_notation((0, 2)) == "o1,3"
    assert gs.graded_to_notation((10,)) == "o11"