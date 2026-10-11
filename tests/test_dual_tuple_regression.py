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
from hypercomplex_engine import multiply, FastDual, DualHolographic
from hypercomplex_engine.core.validation import Validation







# ---------- facade routes, including the path that crashed ----------

@pytest.mark.parametrize("engine", ["fast", "holographic"])
def test_facade_dual_both_engines(engine):
    assert multiply("dual", (1, 3,0), (1, 5,0), dim=3, engine=engine) == \
           multiply("dual", (1, 3,0), (1, 5,0), dim=3, engine="fast")




# ---------- zero handling: validated, shape-correct, convention documented ----------

def test_facade_zero_shapes():
    assert multiply("dual", (0, 2, 0), (1, 3, 0), dim=2) == (0, 0, 0)
    assert multiply("dual_split", (1, 2, 0), (0, 3, 0), dim=2) == (0, 0, 0)
    with pytest.raises((TypeError, ValueError)):
        assert multiply("standard", (0, 2, 0), (1, 3, 0), dim=2) == (0, 0) 
    with pytest.raises((TypeError, ValueError)):  # formal zero
        assert multiply("split", (0, 2, 0), (1, 3, 0), dim=2) == (0, 0)


def test_facade_malformed_raises_validation_error():
    with pytest.raises((TypeError, ValueError)):
        multiply("dual", 5, (1, 1), dim=2)              # not a tuple
    with pytest.raises((TypeError, ValueError)):
        multiply("standard", (1, 2, 3, 4), (1, 1))      # wrong arity

# ---------- nilpotency still intact ----------

@pytest.mark.parametrize("engine", ["fast", "holographic"])
def test_nilpotency(engine):
    assert multiply("dual", (1, 0,1), (1, 0,1), dim=0, engine=engine) == (0, 0, 0)