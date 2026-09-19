# tests/test_mega_mother.py
"""
MEGA MOTHER TEST
================
Exercises the entire hypercomplex package:

    1. BasisNotation        – notation conversion
    2. Validation           – input checking
    3. Table builders       – standard, split, dual
    4. Holographic multipliers – standard, split, dual
    5. Cross-validation     – table vs holographic
    6. CDFormat             – element formatting
    7. CDTablePrinter       – CSV export
    8. Facade               – simple API
    9. End-to-end           – full workflows

Run:
    pytest tests/test_mega_mother.py -v
"""

import os
import tempfile

import numpy as np
import pytest

from hypercomplex import (
    # Facade
    build_table,
    multiply,
    format_element,
    print_table,
    export_csv,
    # Core
    BasisNotation,
    Validation,
    StandardTableBuilder,
    SplitTableBuilder,
    DualTableBuilder,
    StandardHolographic,
    SplitHolographic,
    DualHolographic,
    # Printer
    CDFormat,
    CDTablePrinter,
)


# ======================================================================
# 1. BASIS NOTATION
# ======================================================================

class TestBasisNotation:
    """Integer ↔ graded ↔ LaTeX round-trips."""

    def test_graded_strings(self):
        assert BasisNotation.to_graded_str(0) == "1"
        assert BasisNotation.to_graded_str(1) == "o1"
        assert BasisNotation.to_graded_str(2) == "o2"
        assert BasisNotation.to_graded_str(3) == "o12"
        assert BasisNotation.to_graded_str(5) == "o13"
        assert BasisNotation.to_graded_str(7) == "o123"

    def test_from_graded_strings(self):
        assert BasisNotation.from_graded_str("1") == 0
        assert BasisNotation.from_graded_str("o1") == 1
        assert BasisNotation.from_graded_str("o12") == 3
        assert BasisNotation.from_graded_str("o13") == 5

    def test_roundtrip(self):
        for k in range(16):
            graded = BasisNotation.to_graded_str(k)
            back = BasisNotation.from_graded_str(graded)
            assert back == k, f"Round-trip failed for k={k}"

    def test_latex_graded(self):
        assert BasisNotation.to_latex(0, mode="graded") == "1"
        assert BasisNotation.to_latex(5, mode="graded") == "o_{13}"

    def test_latex_integer(self):
        assert BasisNotation.to_latex(5, mode="integer") == "e_{5}"


# ======================================================================
# 2. VALIDATION
# ======================================================================

class TestValidation:
    """Input validation for basis tuples and dimensions."""

    def test_valid_standard_tuple(self):
        assert Validation.basis_tuple((1, 0))
        assert Validation.basis_tuple((-1, 7))

    def test_valid_dual_tuple(self):
        assert Validation.basis_tuple((1, 3, 1), allow_eps=True)
        assert Validation.basis_tuple((0, 0, 0), allow_zero=True, allow_eps=True)

    def test_reject_bad_sign(self):
        with pytest.raises(ValueError):
            Validation.basis_tuple((2, 0))

    def test_reject_negative_index(self):
        with pytest.raises(ValueError):
            Validation.basis_tuple((1, -1))

    def test_reject_zero_sign_by_default(self):
        with pytest.raises(ValueError):
            Validation.basis_tuple((0, 0))

    def test_reject_eps_when_not_allowed(self):
        with pytest.raises(ValueError):
            Validation.basis_tuple((1, 0, 0), allow_eps=False)

    def test_dimension(self):
        assert Validation.dimension(0) == 0
        assert Validation.dimension(5) == 5
        with pytest.raises(ValueError):
            Validation.dimension(-1)
        with pytest.raises(TypeError):
            Validation.dimension("3")


# ======================================================================
# 3. TABLE BUILDERS
# ======================================================================

class TestStandardTableBuilder:
    """Standard Cayley-Dickson tables."""

    def test_real(self):
        signs, indices = StandardTableBuilder().build(0)
        assert signs.shape == (1, 1)
        assert signs[0, 0] == 1

    def test_complex(self):
        signs, indices = StandardTableBuilder().build(1)
        assert signs.shape == (2, 2)
        assert signs[1, 1] == -1   # i² = -1
        assert indices[1, 1] == 0

    def test_quaternion_anticommutativity(self):
        signs, indices = StandardTableBuilder().build(2)
        # i*j = k
        assert signs[1, 2] == 1
        assert indices[1, 2] == 3
        # j*i = -k
        assert signs[2, 1] == -1
        assert indices[2, 1] == 3

    def test_octonion_diagonal(self):
        signs, indices = StandardTableBuilder().build(3)
        for i in range(1, 8):
            assert signs[i, i] == -1
            assert indices[i, i] == 0

    def test_identity_row_column(self):
        signs, indices = StandardTableBuilder().build(3)
        dim = 8
        for i in range(dim):
            assert signs[0, i] == 1
            assert signs[i, 0] == 1
            assert indices[0, i] == i
            assert indices[i, 0] == i


class TestSplitTableBuilder:
    """Split Cayley-Dickson tables."""

    def test_split_complex_square(self):
        signs, _ = SplitTableBuilder().build(1)
        assert signs[1, 1] == 1   # j² = +1

    def test_split_quaternion_new_elements(self):
        signs, _ = SplitTableBuilder().build(2)
        # Old elements: e1² = -1
        assert signs[1, 1] == -1
        # New elements: e2² = e3² = +1
        assert signs[2, 2] == 1
        assert signs[3, 3] == 1

    def test_split_block_d_first_row(self):
        """Block d first row (i_loc=0, j_loc>0) should be -1 in split."""
        signs, indices = SplitTableBuilder().build(2)
        half = 2
        # e2 * e3: i_loc=0, j_loc=1 → first row of block d
        assert signs[half, half + 1] == -1


class TestDualTableBuilder:
    """Dual extension tables."""

    def test_dual_real_nilpotency(self):
        signs, indices, eps = DualTableBuilder().build(0, split=False)
        assert signs.shape == (2, 2)
        assert signs[1, 1] == 0   # ε² = 0

    def test_dual_complex_shape(self):
        signs, indices, eps = DualTableBuilder().build(1, split=False)
        assert signs.shape == (4, 4)
        assert eps.shape == (4, 4)

    def test_dual_eps_flags(self):
        signs, indices, eps = DualTableBuilder().build(1, split=False)
        # Lower-left block: no epsilon
        assert eps[0, 0] == 0
        assert eps[1, 1] == 0
        # Upper-right block: has epsilon
        assert eps[0, 2] == 1
        # Upper-upper block: zero product
        assert signs[2, 2] == 0


# ======================================================================
# 4. HOLOGRAPHIC MULTIPLIERS
# ======================================================================

class TestStandardHolographic:
    """O(n) standard multiplier."""

    def test_identity(self):
        holo = StandardHolographic()
        assert holo.multiply((1, 0), (1, 5)) == (1, 5)
        assert holo.multiply((1, 5), (1, 0)) == (1, 5)

    def test_complex_square(self):
        holo = StandardHolographic()
        assert holo.multiply((1, 1), (1, 1)) == (-1, 0)

    def test_quaternion_product(self):
        holo = StandardHolographic()
        assert holo.multiply((1, 1), (1, 2)) == (1, 3)   # i*j = k
        assert holo.multiply((1, 2), (1, 1)) == (-1, 3)  # j*i = -k

    def test_sign_composition(self):
        holo = StandardHolographic()
        # (-i) * j = -(i*j) = -k
        assert holo.multiply((-1, 1), (1, 2)) == (-1, 3)


class TestSplitHolographic:
    """O(n) split multiplier."""

    def test_split_complex_square(self):
        holo = SplitHolographic()
        assert holo.multiply((1, 1), (1, 1), dim=1) == (1, 0)

    def test_split_new_elements_square_positive(self):
        holo = SplitHolographic()
        assert holo.multiply((1, 2), (1, 2), dim=2) == (1, 0)
        assert holo.multiply((1, 3), (1, 3), dim=2) == (1, 0)

    def test_split_old_elements_square_negative(self):
        holo = SplitHolographic()
        assert holo.multiply((1, 1), (1, 1), dim=2) == (-1, 0)


class TestDualHolographic:
    """O(n) dual multiplier."""

    def test_nilpotency(self):
        holo = DualHolographic(split=False)
        result = holo.multiply((1, 1), (1, 1), dim=0)
        assert result[0] == 0   # zero element

    def test_base_product_no_eps(self):
        holo = DualHolographic(split=False)
        result = holo.multiply((1, 1), (1, 1), dim=1)
        assert result[0] == -1  # sign
        assert result[2] == 0   # no epsilon

    def test_one_eps_factor(self):
        holo = DualHolographic(split=False)
        # e0 * (eps*e0) = eps*e0
        result = holo.multiply((1, 0), (1, 2), dim=1)
        assert result[2] == 1   # has epsilon

    def test_zero_propagation(self):
        holo = DualHolographic(split=False)
        result = holo.multiply((0, 0), (1, 1), dim=1)
        assert result[0] == 0


# ======================================================================
# 5. CROSS-VALIDATION: TABLE vs HOLOGRAPHIC
# ======================================================================

class TestCrossValidation:
    """
    The most important test class.
    Proves the O(n) holographic multiplier produces identical results
    to the O(4^n) table builder.
    """

    def test_standard_cross_validation(self):
        builder = StandardTableBuilder()
        holo = StandardHolographic()

        for n in range(4):
            signs, indices = builder.build(n)
            dim = 1 << n

            for i in range(dim):
                for j in range(dim):
                    h_sign, h_idx = holo.multiply_indices(i, j)
                    assert int(signs[i, j]) == h_sign, \
                        f"Sign mismatch at n={n}, ({i},{j})"
                    assert int(indices[i, j]) == h_idx, \
                        f"Index mismatch at n={n}, ({i},{j})"

    def test_split_cross_validation(self):
        builder = SplitTableBuilder()
        holo = SplitHolographic()

        for n in range(1, 4):
            signs, indices = builder.build(n)
            dim = 1 << n

            for i in range(dim):
                for j in range(dim):
                    h_sign, h_idx = holo.multiply_indices(i, j, dim=n)
                    assert int(signs[i, j]) == h_sign, \
                        f"Split sign mismatch at n={n}, ({i},{j})"
                    assert int(indices[i, j]) == h_idx, \
                        f"Split index mismatch at n={n}, ({i},{j})"


# ======================================================================
# 6. CD FORMAT
# ======================================================================

class TestCDFormat:
    """Element formatting across all modes."""

    def test_integer_mode(self):
        assert CDFormat.format_element((1, 3), mode="integer") == "+e3"
        assert CDFormat.format_element((-1, 3), mode="integer") == "-e3"

    def test_graded_mode(self):
        assert CDFormat.format_element((1, 5), mode="graded") == "+o13"
        assert CDFormat.format_element((-1, 7), mode="graded") == "-o123"

    def test_latex_graded_mode(self):
        assert CDFormat.format_element((1, 5), mode="latex_graded") == "+o_{13}"

    def test_latex_integer_mode(self):
        assert CDFormat.format_element((1, 5), mode="latex_integer") == "+e_{5}"

    def test_dual_formatting(self):
        assert CDFormat.format_element((1, 2, 1), mode="integer") == "+e2*eps"
        assert CDFormat.format_element((1, 0, 1), mode="integer") == "+eps"
        assert CDFormat.format_element((1, 2, 1), mode="latex") == "+e_{2}\\epsilon"

    def test_zero_element(self):
        assert CDFormat.format_element((0, 0), mode="integer") == "0"
        assert CDFormat.format_element((0, 0, 1), mode="integer") == "0"

    def test_basis_labels_standard(self):
        labels = CDFormat.basis_labels(4, dual=False, mode="integer")
        assert labels == ["e0", "e1", "e2", "e3"]

    def test_basis_labels_dual(self):
        labels = CDFormat.basis_labels(4, dual=True, mode="integer")
        assert labels == ["e0", "e1", "eps", "e1*eps"]


# ======================================================================
# 7. CD TABLE PRINTER – CSV EXPORT
# ======================================================================

class TestCDTablePrinter:
    """CSV export in matrix and long modes."""

    def test_export_matrix_standard(self):
        signs, indices = StandardTableBuilder().build(2)

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            path = f.name

        try:
            CDTablePrinter.export_csv(
                path, signs, indices,
                mode="integer", csv_mode="matrix",
            )
            assert os.path.exists(path)

            with open(path) as f:
                content = f.read()
            assert "e0" in content
            assert "+e0" in content
        finally:
            os.unlink(path)

    def test_export_long_standard(self):
        signs, indices = StandardTableBuilder().build(1)

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            path = f.name

        try:
            CDTablePrinter.export_csv(
                path, signs, indices,
                mode="integer", csv_mode="long",
            )
            with open(path) as f:
                lines = f.readlines()

            assert lines[0].strip() == "i,j,sign,index"
            assert len(lines) == 5  # header + 4 entries
        finally:
            os.unlink(path)

    def test_export_dual_includes_eps_column(self):
        signs, indices, eps = DualTableBuilder().build(1, split=False)

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            path = f.name

        try:
            CDTablePrinter.export_csv(
                path, signs, indices, eps=eps,
                mode="integer", csv_mode="long",
            )
            with open(path) as f:
                header = f.readline().strip()

            assert header == "i,j,sign,index,eps"
        finally:
            os.unlink(path)


# ======================================================================
# 8. FACADE API
# ======================================================================

class TestFacade:
    """The simple user-facing API."""

    def test_build_table_standard(self):
        signs, indices = build_table("standard", 2)
        assert signs.shape == (4, 4)

    def test_build_table_split(self):
        signs, indices = build_table("split", 2)
        assert signs.shape == (4, 4)

    def test_build_table_dual(self):
        signs, indices, eps = build_table("dual", 2)
        assert signs.shape == (8, 8)

    def test_build_table_dual_split(self):
        signs, indices, eps = build_table("dual_split", 2)
        assert signs.shape == (8, 8)

    def test_multiply_standard(self):
        assert multiply("standard", (1, 1), (1, 2)) == (1, 3)

    def test_multiply_split(self):
        assert multiply("split", (1, 1), (1, 1), dim=1) == (1, 0)

    def test_multiply_dual(self):
        result = multiply("dual", (1, 1), (1, 1), dim=1)
        assert result[0] == -1
        assert result[2] == 0

    def test_format_element(self):
        assert format_element((1, 5), mode="graded") == "+o13"

    def test_export_csv_via_facade(self):
        table = build_table("standard", 2)

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            path = f.name

        try:
            result = export_csv(path, table, mode="integer")
            assert os.path.exists(result)
        finally:
            os.unlink(path)

    def test_invalid_kind_raises(self):
        with pytest.raises(ValueError):
            build_table("bogus", 2)


# ======================================================================
# 9. END-TO-END WORKFLOWS
# ======================================================================

class TestEndToEnd:
    """Full workflows combining every layer."""

    def test_octonion_workflow(self):
        """Build → multiply → format → export."""
        table = build_table("standard", 3)

        result = multiply("standard", (1, 1), (1, 2))
        assert result == (1, 3)

        assert format_element(result, mode="integer") == "+e3"
        assert format_element(result, mode="graded") == "+o3"

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            path = f.name
        try:
            export_csv(path, table, mode="graded", csv_mode="matrix")
            assert os.path.exists(path)
        finally:
            os.unlink(path)

    def test_split_octonion_workflow(self):
        """New elements square to +1."""
        result = multiply("split", (1, 4), (1, 4), dim=3)
        assert result[0] == 1

    def test_dual_nilpotency_workflow(self):
        """ε² = 0 in dual algebra."""
        result = multiply("dual", (1, 4), (1, 4), dim=2)
        assert result[0] == 0

    def test_notation_roundtrip_in_context(self):
        """Multiply, format to graded, parse back."""
        result = multiply("standard", (1, 1), (1, 2))
        graded = format_element(result, mode="graded")
        # Strip the sign prefix to get the basis label
        basis_label = graded[1:]  # remove '+' or '-'
        recovered_index = BasisNotation.from_graded_str(basis_label)
        assert recovered_index == result[1]

    def test_all_modes_produce_output(self):
        """Every formatting mode produces a non-empty string."""
        element = (1, 5)
        for mode in CDFormat.VALID_MODES:
            text = CDFormat.format_element(element, mode=mode)
            assert len(text) > 0
            assert text[0] in ("+", "-", "0")


# ======================================================================
# RUNNER
# ======================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])