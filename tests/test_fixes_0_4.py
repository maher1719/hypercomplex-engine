"""Regression tests for the 0.4.0 fixes."""
import itertools

import pytest

from hypercomplex import (
    DualTableBuilder,
    FastDual,
    DualHolographic,
    build_table,
    estimate_table_bytes,
    format_element,
    multiply,
)


# --- 1. dual nilpotency: one canonical zero across every engine -----------

@pytest.mark.parametrize("split", [False, True])
@pytest.mark.parametrize("n", [1, 2, 3])
def test_dual_engines_agree_everywhere(n, split):
    fast = FastDual(split=split)
    holo = DualHolographic(split=split)
    signs, indices, eps = DualTableBuilder().build(n, split=split)
    size = 1 << (n + 1)

    for i, j in itertools.product(range(size), repeat=2):
        f = tuple(fast.multiply((1, i), (1, j), dim=n))
        h = tuple(holo.multiply((1, i), (1, j), dim=n))
        assert f == h, (i, j, f, h)

        if signs[i, j] == 0:
            assert f == (0, 0, 0)
        else:
            assert f == (int(signs[i, j]), int(indices[i, j]), int(eps[i, j]))


@pytest.mark.parametrize("kind", ["dual", "dual_split"])
def test_facade_eps_squared_is_canonical_zero(kind):
    eps = (1, 0, 1)
    for engine in ("fast", "holographic"):
        assert multiply(kind, eps, eps, dim=1, engine=engine) == (0, 0, 0)
    assert format_element((0, 0, 0), mode="integer") == "0"


# --- 2. build_table: strict n, memory guard --------------------------------

@pytest.mark.parametrize("bad", [2.5, 2.0, "3", None, True])
def test_build_table_rejects_non_integer_n(bad):
    with pytest.raises(TypeError):
        build_table("standard", bad)


def test_build_table_rejects_negative_n():
    with pytest.raises(ValueError):
        build_table("standard", -1)


def test_build_table_guard_blocks_huge_request_without_allocating():
    with pytest.raises(ValueError, match="max_bytes"):
        build_table("standard", 40)
    with pytest.raises(ValueError, match="max_bytes"):
        build_table("dual_split", 30)


def test_build_table_max_bytes_can_be_raised_or_disabled():
    assert build_table("standard", 6, max_bytes=None)[0].shape == (64, 64)
    with pytest.raises(ValueError):
        build_table("standard", 6, max_bytes=100)
    assert build_table("standard", 6, max_bytes=estimate_table_bytes("standard", 6))[0].shape == (64, 64)


@pytest.mark.parametrize("kind", ["standard", "split", "dual", "dual_split"])
@pytest.mark.parametrize("n", [0, 1, 3, 6, 9])
def test_estimate_matches_real_allocation(kind, n):
    table = build_table(kind, n, max_bytes=None)
    actual = sum(arr.nbytes for arr in table)
    assert estimate_table_bytes(kind, n) == actual


def test_default_budget_allows_documented_sizes():
    assert estimate_table_bytes("standard", 13) <= 1 << 28 < estimate_table_bytes("standard", 14)
    assert estimate_table_bytes("dual", 12) <= 1 << 28 < estimate_table_bytes("dual", 13)


# --- 3. engine aliases ------------------------------------------------------

@pytest.mark.parametrize("alias", ["fast", "o1", "bitwise", "constant"])
def test_fast_engine_aliases(alias):
    assert multiply("standard", (1, 3), (1, 5), engine=alias) == multiply(
        "standard", (1, 3), (1, 5), engine="fast"
    )


@pytest.mark.parametrize("alias", ["holographic", "on", "o(n)", "descent"])
def test_holographic_engine_aliases(alias):
    assert multiply("standard", (1, 3), (1, 5), engine=alias) == multiply(
        "standard", (1, 3), (1, 5), engine="holographic"
    )


# --- 4. latex documented behaviour -----------------------------------------

def test_latex_is_alias_of_latex_integer():
    for el in [(1, 3), (-1, 5), (1, 0)]:
        assert format_element(el, mode="latex") == format_element(el, mode="latex_integer")
    assert format_element((1, 3), mode="latex") == "+e_{3}"
    assert format_element((-1, 5), mode="latex_graded") == "-o_{13}"
