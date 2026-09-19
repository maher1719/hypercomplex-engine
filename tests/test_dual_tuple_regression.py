"""Regression tests for dual-engine tuple normalization and facade zero handling.

History:
- 0.4.0: DualHolographic ignored the eps flag while FastDual promoted it
  (engines diverged on identical 3-tuple input).
- 0.4.0 hotfix attempt: DualHolographic read t[2] unconditionally, crashing
  with IndexError on documented 2-tuple input and breaking the facade's
  engine="holographic" dual path (which normalizes to 2-tuples).
- 0.4.1: both fixed; facade validates before zero-propagation.
"""

import pytest

# Adjust these imports to whatever hypercomplex/__init__.py exports.
from hypercomplex import multiply, FastDual, DualHolographic
from hypercomplex.core.validation import Validation


# ---------- tuple normalization: engines must agree on identical input ----------

def _all_global_pairs(dim):
    half = 1 << dim
    for i in range(half << 1):
        for j in range(half << 1):
            yield (1, i, 1 if i >= half else 0), (1, j, 1 if j >= half else 0)


@pytest.mark.parametrize("dim", [0, 1, 2, 3])
def test_engines_agree_global_3tuples(dim):
    fast, holo = FastDual(), DualHolographic()
    for a, b in _all_global_pairs(dim):
        assert fast.multiply(a, b, dim) == holo.multiply(a, b, dim)


@pytest.mark.parametrize("dim", [1, 2, 3])
def test_engines_agree_on_2tuples(dim):
    """2-tuples are documented input; must not raise (the 0.4.0 regression)."""
    fast, holo = FastDual(), DualHolographic()
    for a, b in _all_global_pairs(dim):
        a2, b2 = a[:2], b[:2]                      # strip the eps flag
        assert fast.multiply(a2, b2, dim) == holo.multiply(a2, b2, dim)


@pytest.mark.parametrize("dim", [1, 2, 3])
def test_local_3tuple_matches_promoted_global(dim):
    """(sign, local_index, eps=1) must equal (sign, local_index + half)."""
    fast, holo = FastDual(), DualHolographic()
    half = 1 << dim
    for loc in range(half):
        for j in range(half << 1):
            for engine in (fast, holo):
                assert engine.multiply((1, loc, 1), (1, j), dim) == \
                       engine.multiply((1, loc + half), (1, j), dim)


# ---------- facade routes, including the path that crashed ----------

@pytest.mark.parametrize("engine", ["fast", "holographic"])
def test_facade_dual_both_engines(engine):
    assert multiply("dual", (1, 3), (1, 5), dim=3, engine=engine) == \
           multiply("dual", (1, 3), (1, 5), dim=3, engine="fast")


def test_facade_dual_3tuple_promotion():
    half = 1 << 3
    for engine in ("fast", "holographic"):
        assert multiply("dual", (1, 2, 1), (1, 1), dim=3, engine=engine) == \
               multiply("dual", (1, 2 + half), (1, 1), dim=3, engine=engine)


# ---------- zero handling: validated, shape-correct, convention documented ----------

def test_facade_zero_shapes():
    assert multiply("dual", (0, 2), (1, 3), dim=2) == (0, 0, 0)
    assert multiply("dual_split", (1, 2), (0, 3), dim=2) == (0, 0, 0)
    assert multiply("standard", (0, 2), (1, 3), dim=2) == (0, 0)   # formal zero
    assert multiply("split", (0, 2), (1, 3), dim=2) == (0, 0)


def test_facade_malformed_raises_validation_error():
    with pytest.raises((TypeError, ValueError)):
        multiply("dual", 5, (1, 1), dim=2)              # not a tuple
    with pytest.raises((TypeError, ValueError)):
        multiply("standard", (1, 2, 3, 4), (1, 1))      # wrong arity

# ---------- nilpotency still intact ----------

@pytest.mark.parametrize("engine", ["fast", "holographic"])
def test_nilpotency(engine):
    assert multiply("dual", (1, 1), (1, 1), dim=0, engine=engine) == (0, 0, 0)