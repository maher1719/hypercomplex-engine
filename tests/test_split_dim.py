"""split/dual require an explicit dim; split rejects epsilon 3-tuples."""
import pytest
from hypercomplex import multiply


@pytest.mark.parametrize("engine", ["fast", "holographic"])
@pytest.mark.parametrize("kind", ["split", "dual", "dual_split"])
def test_dim_is_required_and_message_names_the_kind(kind, engine):
    with pytest.raises(ValueError, match=f"dim is required for {kind}"):
        multiply(kind, (1, 1), (1, 1), engine=engine)


@pytest.mark.parametrize("kind", ["split", "dual", "dual_split"])
def test_missing_dim_is_reported_even_for_zero_elements(kind):
    with pytest.raises(ValueError, match="dim is required"):
        multiply(kind, (0, 1), (1, 2))


@pytest.mark.parametrize("engine", ["fast", "holographic"])
def test_split_result_depends_on_dim(engine):
    assert multiply("split", (1, 2), (1, 2), dim=2, engine=engine) == (1, 0)
    for dim in (3, 4, 8):
        assert multiply("split", (1, 2), (1, 2), dim=dim, engine=engine) == (-1, 0)


@pytest.mark.parametrize("eps", [0, 1])
def test_split_rejects_epsilon_three_tuples_with_clear_message(eps):
    with pytest.raises(ValueError, match="3-tuple"):
        multiply("split", (1, 1, eps), (1, 2), dim=2)


def test_standard_needs_no_dim():
    assert multiply("standard", (1, 1), (1, 2)) == (1, 3)
