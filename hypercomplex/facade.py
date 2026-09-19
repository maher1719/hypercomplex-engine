from .core.table_builder import (
    StandardTableBuilder,
    SplitTableBuilder,
    DualTableBuilder,
)

from .core.holographic import (
    StandardHolographic,
    SplitHolographic,
    DualHolographic,
)

from .printer import CDFormat, CDTablePrinter


# ----------------------------------------------------------------------
# Internal singletons
# ----------------------------------------------------------------------

_standard_table = StandardTableBuilder()
_split_table = SplitTableBuilder()
_dual_table = DualTableBuilder()

_standard_holo = StandardHolographic()
_split_holo = SplitHolographic()
_dual_holo = DualHolographic(split=False)
_dual_split_holo = DualHolographic(split=True)


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def _normalize_kind(kind: str) -> str:
    if not isinstance(kind, str):
        raise TypeError("kind must be a string")

    kind = kind.strip().lower()

    aliases = {
        "standard": "standard",
        "std": "standard",
        "ordinary": "standard",
        "o": "standard",

        "split": "split",
        "s": "split",

        "dual": "dual",
        "dual_standard": "dual",
        "d": "dual",

        "dual_split": "dual_split",
        "split_dual": "dual_split",
        "ds": "dual_split",
    }

    if kind not in aliases:
        raise ValueError(
            f"Unknown algebra kind '{kind}'. "
            "Valid kinds: standard, split, dual, dual_split"
        )

    return aliases[kind]


def _unpack_table(table):
    """
    Accepts:
        (signs, indices)
        (signs, indices, eps)
    """
    if len(table) == 2:
        signs, indices = table
        return signs, indices, None

    if len(table) == 3:
        signs, indices, eps = table
        return signs, indices, eps

    raise ValueError(
        "table must be (signs, indices) or (signs, indices, eps)"
    )


# ----------------------------------------------------------------------
# Public simple API
# ----------------------------------------------------------------------

def build_table(kind: str, n: int):
    """
    Build a multiplication table.

    kind:
        "standard"
        "split"
        "dual"
        "dual_split"

    n:
        Parent dimension exponent.

        standard:
            dimension = 2^n

        split:
            dimension = 2^n

        dual / dual_split:
            dimension = 2^(n+1)

    Returns:
        standard/split:
            (signs, indices)

        dual/dual_split:
            (signs, indices, eps)
    """
    kind = _normalize_kind(kind)
    n = int(n)

    if kind == "standard":
        return _standard_table.build(n)

    if kind == "split":
        return _split_table.build(n)

    if kind == "dual":
        return _dual_table.build(n, split=False)

    if kind == "dual_split":
        return _dual_table.build(n, split=True)

    raise ValueError(f"Unknown algebra kind '{kind}'")


def multiply(kind: str, a: tuple, b: tuple, dim: int | None = None) -> tuple:
    """
    Multiply two basis elements.

    kind:
        "standard"
        "split"
        "dual"
        "dual_split"

    a, b:
        Basis tuples.

        standard/split:
            (sign, index)

        dual:
            (sign, index) or (sign, index, eps)

    dim:
        Required for dual and dual_split.
        Optional for split.

        For split, if dim is not given, it is inferred from the indices.

    Returns:
        standard/split:
            (sign, index)

        dual/dual_split:
            (sign, local_index, eps_flag)
    """
    kind = _normalize_kind(kind)

    if kind == "standard":
        return _standard_holo.multiply(a, b)

    if kind == "split":
        if dim is None:
            dim = max(int(a[1]), int(b[1])).bit_length()
        return _split_holo.multiply(a, b, dim)

    if kind == "dual":
        if dim is None:
            raise ValueError("dim is required for dual multiplication")
        return _dual_holo.multiply(a, b, dim)

    if kind == "dual_split":
        if dim is None:
            raise ValueError("dim is required for dual_split multiplication")
        return _dual_split_holo.multiply(a, b, dim)

    raise ValueError(f"Unknown algebra kind '{kind}'")


def format_element(element: tuple, mode: str = "integer") -> str:
    """
    Format a basis element tuple.

    element:
        (sign, index)
        (sign, index, eps)

    mode:
        "integer"
        "graded"
        "latex"
        "latex_integer"
        "latex_graded"
    """
    return CDFormat.format_element(element, mode=mode)


def print_table(
    table,
    title: str | None = None,
    limit: int | None = None,
    mode: str = "integer",
):
    """
    Print a table built by build_table().

    table:
        (signs, indices)
        or
        (signs, indices, eps)
    """
    signs, indices, eps = _unpack_table(table)

    return CDTablePrinter.print_table(
        signs,
        indices,
        eps=eps,
        title=title,
        limit=limit,
        mode=mode,
    )


def export_csv(
    path: str,
    table,
    mode: str = "integer",
    csv_mode: str = "matrix",
) -> str:
    """
    Export a table built by build_table().

    table:
        (signs, indices)
        or
        (signs, indices, eps)

    csv_mode:
        "matrix"
        "long"
    """
    signs, indices, eps = _unpack_table(table)

    return CDTablePrinter.export_csv(
        path,
        signs,
        indices,
        eps=eps,
        mode=mode,
        csv_mode=csv_mode,
    )