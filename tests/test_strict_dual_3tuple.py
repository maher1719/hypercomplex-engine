"""
Phase 1 Validation: Strict 3-Tuple Contract for Dual Algebras.

This test suite asserts that:
1. The old 2-tuple "global overflow" format (e.g., (1, 4) for dim=2) is REJECTED.
2. The new strict 3-tuple format (sign, local_index, eps_flag) is ENFORCED.
3. The mathematical invariants (nilpotency, parent matching) remain perfectly intact.
4. The Fast O(1) and Holographic O(n) engines agree perfectly under the new contract.
"""
import pytest
from hypercomplex_engine import (
    multiply, 
    FastDual, DualHolographic, 
    FastStandard, FastSplit
)


class TestFacadeStrictness:
    """Validates that the top-level multiply() function enforces tuple lengths."""

    def test_facade_rejects_2tuple_overflow_for_dual(self):
        # Old behavior: (1, 4) with dim=2 meant eps*e0 (since 4 >= 2^2).
        # New behavior: Must raise ValueError because it's a 2-tuple.
        with pytest.raises(ValueError, match="must have length 3"):
            multiply("dual", (1, 4), (1, 0), dim=2)
            
    def test_facade_rejects_2tuple_overflow_for_dual_split(self):
        with pytest.raises(ValueError, match="must have length 3"):
            multiply("dual_split", (1, 5), (1, 1), dim=2)

    def test_facade_rejects_3tuple_for_standard(self):
        # Standard algebra must strictly be 2-tuples.
        with pytest.raises(ValueError, match="must have length 2"):
            multiply("standard", (1, 1, 0), (1, 2))

    def test_facade_rejects_3tuple_for_split(self):
        with pytest.raises(ValueError, match="must have length 2"):
            multiply("split", (1, 1, 0), (1, 2), dim=2)

    def test_facade_accepts_3tuple_for_dual(self):
        # e1 * eps = eps*e1
        result = multiply("dual", (1, 1, 0), (1, 0, 1), dim=2)
        assert result == (1, 1, 1)

    def test_facade_accepts_3tuple_for_dual_split(self):
        # split e1 * e1 = +e0 -> eps * e1 * e1 = +eps*e0
        result = multiply("dual_split", (1, 1, 1), (1, 1, 0), dim=1)
        assert result == (1, 0, 1)


class TestEngineStrictness:
    """Validates that the low-level engines reject out-of-bounds local indices."""

    @pytest.mark.parametrize("split", [False, True])
    def test_fast_rejects_out_of_bounds_local_index(self, split):
        fast = FastDual(split=split)
        # dim=2 -> half=4. local_index must be < 4.
        # Passing 4 should raise ValueError, not silently treat it as eps.
        with pytest.raises(ValueError, match="local index must be in"):
            fast.multiply_indices(4, 1, dim=2, i_eps_flag=0, j_eps_flag=0)

    @pytest.mark.parametrize("split", [False, True])
    def test_holographic_rejects_out_of_bounds_local_index(self, split):
        holo = DualHolographic(split=split)
        with pytest.raises(ValueError, match="local index must be in"):
            holo.multiply((1, 4, 0), (1, 1, 0), dim=2)


class TestMathematicalInvariants:
    """Ensures the core math hasn't broken during the migration."""

    @pytest.mark.parametrize("split", [False, True])
    @pytest.mark.parametrize("dim", [1, 2, 3])
    def test_nilpotency_eps_squared_is_zero(self, split, dim):
        fast = FastDual(split=split)
        holo = DualHolographic(split=split)
        
        half = 1 << dim
        for i in range(half):
            for j in range(half):
                # (eps * e_i) * (eps * e_j) = 0
                f_res = fast.multiply_indices(i, j, dim=dim, i_eps_flag=1, j_eps_flag=1)
                h_res = holo.multiply((1, i, 1), (1, j, 1), dim=dim)
                
                assert f_res == (0, 0, 0), f"Fast nilpotency failed for i={i}, j={j}"
                assert h_res == (0, 0, 0), f"Holo nilpotency failed for i={i}, j={j}"

    @pytest.mark.parametrize("split", [False, True])
    @pytest.mark.parametrize("dim", [1, 2, 3])
    def test_dual_product_matches_parent(self, split, dim):
        """
        The core mathematical invariant: 
        (e_i * e_j) in dual == (e_i * e_j) in parent, with eps_flag = 0.
        (eps * e_i) * e_j == eps * (e_i * e_j)
        """
        fast_dual = FastDual(split=split)
        fast_parent = FastSplit() if split else FastStandard()
        
        half = 1 << dim
        for i in range(half):
            for j in range(half):
                # Get parent product
                if split:
                    p_sign, p_idx = fast_parent.multiply_indices(i, j, dim=dim)
                else:
                    p_sign, p_idx = fast_parent.multiply_indices(i, j)
                
                # 1. Base * Base
                d_res = fast_dual.multiply_indices(i, j, dim=dim, i_eps_flag=0, j_eps_flag=0)
                assert d_res == (p_sign, p_idx, 0)
                
                # 2. Eps * Base
                d_eps_base = fast_dual.multiply_indices(i, j, dim=dim, i_eps_flag=1, j_eps_flag=0)
                assert d_eps_base == (p_sign, p_idx, 1)
                
                # 3. Base * Eps
                d_base_eps = fast_dual.multiply_indices(i, j, dim=dim, i_eps_flag=0, j_eps_flag=1)
                assert d_base_eps == (p_sign, p_idx, 1)


class TestCrossValidation:
    """Ensures Fast O(1) and Holographic O(n) agree under the new contract."""

    @pytest.mark.parametrize("split", [False, True])
    def test_fast_vs_holographic_exhaustive(self, split):
        fast = FastDual(split=split)
        holo = DualHolographic(split=split)
        dim = 3
        half = 1 << dim
        
        for i in range(half):
            for j in range(half):
                for i_eps in (0, 1):
                    for j_eps in (0, 1):
                        f_res = fast.multiply_indices(i, j, dim=dim, i_eps_flag=i_eps, j_eps_flag=j_eps)
                        h_res = holo.multiply((1, i, i_eps), (1, j, j_eps), dim=dim)
                        
                        # Holographic multiply returns (s1*s2*base_sign, idx, eps)
                        # Since s1=1, s2=1, it should exactly match fast.
                        assert f_res == h_res, (
                            f"Mismatch at i={i}, j={j}, i_eps={i_eps}, j_eps={j_eps}. "
                            f"Fast: {f_res}, Holo: {h_res}"
                        )
