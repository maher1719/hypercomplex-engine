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

from .core.fast import (
    FastStandard,
    FastSplit,
    FastDual,
)

from .core.validation import Validation

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

_standard_fast = FastStandard()
_split_fast = FastSplit()
_dual_fast = FastDual(split=False)
_dual_split_fast = FastDual(split=True)


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


def _normalize_engine(engine: str) -> str:
    if not isinstance(engine, str):
        raise TypeError("engine must be a string")

    engine = engine.strip().lower()

    aliases = {
        "fast": "fast",
        "constant": "fast",
        "fast": "fast",
        "bitwise": "fast",

        "holographic": "holographic",
        "on": "holographic",
        "o(n)": "holographic",
        "descent": "holographic",
    }

    if engine not in aliases:
        raise ValueError(
            f"Unknown engine '{engine}'. "
            "Valid engines: fast, holographic"
        )

    return aliases[engine]


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


def _as_dual_global(t: tuple, dim: int) -> tuple:
    """
    Converts a dual input tuple into a global 2-tuple.

    Accepts:
        (sign, global_index)
        (sign, local_index, eps_flag)

    Returns:
        (sign, global_index)
    """
    sign = int(t[0])
    idx = int(t[1])

    if sign == 0:
        return (0, 0)

    half = 1 << Validation.dimension(dim)

    if len(t) >= 3:
        eps = int(t[2])

        if eps == 1 and idx < half:
            idx += half

    return (sign, idx)


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
    """
    kind = _normalize_kind(kind)
    n = int(n)

    if kind == "standard":
        return _standard_table.build(n)

    if kind == "split":
        return _split_table.build(n,split=True)

    if kind == "dual":
        return _dual_table.build(n, split=False)

    if kind == "dual_split":
        return _dual_table.build(n, split=True)

    raise ValueError(f"Unknown algebra kind '{kind}'")


def multiply(
    kind: str,
    a: tuple,
    b: tuple,
    dim: int | None = None,
    engine: str = "fast",
) -> tuple:
    """
    Multiply two basis elements.

    kind:
        "standard"
        "split"
        "dual"
        "dual_split"

    engine:
        "fast"          default, fastest
        "holographic" O(n) descent, useful for verification

    dim:
        Required for dual and dual_split.
        Optional for split.
    """
    kind = _normalize_kind(kind)
    engine = _normalize_engine(engine)

    # Zero propagation
    if int(a[0]) == 0 or int(b[0]) == 0:
        if kind in ("dual", "dual_split"):
            return (0, 0, 0)
        return (0, 0)

    # Dual inputs may be local 3-tuples: (sign, local_index, eps)
    if kind in ("dual", "dual_split"):
        if dim is None:
            raise ValueError("dim is required for dual multiplication")

        a = _as_dual_global(a, dim)
        b = _as_dual_global(b, dim)

    # ------------------------------------------------------------------
    # O(1) engine
    # ------------------------------------------------------------------
    if engine == "fast":
        if kind == "standard":
            return _standard_fast.multiply(a, b)

        if kind == "split":
            return _split_fast.multiply(a, b, dim)

        if kind == "dual":
            return _dual_fast.multiply(a, b, dim)

        if kind == "dual_split":
            return _dual_split_fast.multiply(a, b, dim)

    # ------------------------------------------------------------------
    # Holographic O(n) engine
    # ------------------------------------------------------------------
    if engine == "holographic":
        if kind == "standard":
            return _standard_holo.multiply(a, b)

        if kind == "split":
            if dim is None:
                dim = max(int(a[1]), int(b[1])).bit_length()
            return _split_holo.multiply(a, b, dim)

        if kind == "dual":
            return _dual_holo.multiply(a, b, dim)

        if kind == "dual_split":
            return _dual_split_holo.multiply(a, b, dim)

    raise ValueError(f"Unknown engine '{engine}'")


def format_element(element: tuple, mode: str = "integer") -> str:
    """
    Format a basis element tuple.
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