import pytest

from hypercomplex.core.table_builder import (
    StandardTableBuilder,
    SplitTableBuilder,
    DualTableBuilder,
)

from hypercomplex.core.fast import (
    FastStandard,
    FastSplit,
    FastDual,
)

from hypercomplex import multiply


# ======================================================================
# STANDARD O(1)
# ======================================================================

class TestFastStandard:

    def test_identity(self):
        fast = FastStandard()

        assert fast.multiply((1, 0), (1, 5)) == (1, 5)
        assert fast.multiply((1, 5), (1, 0)) == (1, 5)

    def test_complex_square(self):
        fast = FastStandard()

        assert fast.multiply((1, 1), (1, 1)) == (-1, 0)

    def test_quaternion_products(self):
        fast = FastStandard()

        # i * j = k
        assert fast.multiply((1, 1), (1, 2)) == (1, 3)

        # j * i = -k
        assert fast.multiply((1, 2), (1, 1)) == (-1, 3)

    def test_cross_validation_with_table(self):
        builder = StandardTableBuilder()
        fast = FastStandard()

        for n in range(4):
            signs, indices = builder.build(n)
            dim = 1 << n

            for i in range(dim):
                for j in range(dim):
                    h_sign, h_idx = fast.multiply_indices(i, j)

                    assert int(signs[i, j]) == h_sign, \
                        f"O(1) sign mismatch at n={n}, ({i},{j})"

                    assert int(indices[i, j]) == h_idx, \
                        f"O(1) index mismatch at n={n}, ({i},{j})"


# ======================================================================
# SPLIT O(1)
# ======================================================================

class TestFastSplit:

    def test_split_complex_square(self):
        fast = FastSplit()

        assert fast.multiply((1, 1), (1, 1), dim=1) == (1, 0)

    def test_split_new_elements_square_positive(self):
        fast = FastSplit()

        assert fast.multiply((1, 2), (1, 2), dim=2) == (1, 0)
        assert fast.multiply((1, 3), (1, 3), dim=2) == (1, 0)

    def test_split_old_elements_square_negative(self):
        fast = FastSplit()

        assert fast.multiply((1, 1), (1, 1), dim=2) == (-1, 0)

    def test_cross_validation_with_table(self):
        builder = SplitTableBuilder()
        fast = FastSplit()

        for n in range(4):
            signs, indices = builder.build(n)
            dim = 1 << n

            for i in range(dim):
                for j in range(dim):
                    h_sign, h_idx = fast.multiply_indices(i, j, dim=n)

                    assert int(signs[i, j]) == h_sign, \
                        f"Split O(1) sign mismatch at n={n}, ({i},{j})"

                    assert int(indices[i, j]) == h_idx, \
                        f"Split O(1) index mismatch at n={n}, ({i},{j})"


# ======================================================================
# DUAL O(1)
# ======================================================================

class TestFastDual:

    def test_nilpotency(self):
        fast = FastDual(split=False)

        # eps*e0 has global index 2 when dim=1
        result = fast.multiply((1, 2), (1, 2), dim=1)

        assert result == (0, 0, 0)

    def test_base_product_no_eps(self):
        fast = FastDual(split=False)

        # e1 * e1 = -e0
        result = fast.multiply((1, 1), (1, 1), dim=1)

        assert result == (-1, 0, 0)

    def test_one_eps_factor(self):
        fast = FastDual(split=False)

        # e0 * (eps*e0) = eps*e0
        result = fast.multiply((1, 0), (1, 2), dim=1)

        assert result == (1, 0, 1)

    def test_local_tuple_input(self):
        fast = FastDual(split=False)

        # local tuple: (sign, local_index, eps_flag)
        eps_e0 = (1, 0, 1)

        result = fast.multiply((1, 0), eps_e0, dim=1)

        assert result == (1, 0, 1)

    def test_cross_validation_with_table(self):
        builder = DualTableBuilder()

        for split in (False, True):
            fast = FastDual(split=split)

            dim = 2
            signs, indices, eps = builder.build(dim, split=split)

            total = 1 << (dim + 1)

            for i in range(total):
                for j in range(total):
                    s, idx, e = fast.multiply_indices(i, j, dim)

                    assert int(signs[i, j]) == s, \
                        f"Dual O(1) sign mismatch split={split}, ({i},{j})"

                    assert int(indices[i, j]) == idx, \
                        f"Dual O(1) index mismatch split={split}, ({i},{j})"

                    assert int(eps[i, j]) == e, \
                        f"Dual O(1) eps mismatch split={split}, ({i},{j})"


# ======================================================================
# FACADE O(1)
# ======================================================================

class TestFacadefast:

    def test_standard(self):
        assert multiply("standard", (1, 1), (1, 2), engine="fast") == (1, 3)

    def test_split(self):
        assert multiply("split", (1, 1), (1, 1), dim=1, engine="fast") == (1, 0)

    def test_dual(self):
        result = multiply("dual", (1, 1), (1, 1), dim=1, engine="fast")
        assert result == (-1, 0, 0)

    def test_dual_local_tuple_chaining(self):
        # eps*e0 as local tuple
        eps_e0 = (1, 0, 1)

        result = multiply("dual", (1, 0), eps_e0, dim=1, engine="fast")

        assert result == (1, 0, 1)

        # Now reuse result directly
        result2 = multiply("dual", result, (1, 0), dim=1, engine="fast")

        assert result2 == (1, 0, 1)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])