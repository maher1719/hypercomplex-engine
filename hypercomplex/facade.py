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

import numpy as np

from .core.table_builder.common import index_dtype
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
        "bitwise": "fast",
        "o1":"fast",

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

# Default memory budget for build_table (bytes of final table data).
# 256 MiB allows standard/split up to n=13 and dual/dual_split up to n=12.
DEFAULT_MAX_TABLE_BYTES = 1 << 28


def estimate_table_bytes(kind: str, n: int) -> int:
    """
    Estimate the bytes of the final table arrays that build_table(kind, n)
    would allocate. Peak usage while building is roughly 1.5x this value.

    Dual kinds are one doubling larger than the same n for standard/split,
    and carry an extra uint8 epsilon array.
    """
    kind = _normalize_kind(kind)
    n = Validation.dimension(n)

    is_dual = kind in ("dual", "dual_split")
    full = 1 << (n + 1 if is_dual else n)

    per_entry = 1 + np.dtype(index_dtype(full)).itemsize  # signs + indices
    if is_dual:
        per_entry += 1  # eps array

    return full * full * per_entry


def build_table(kind: str, n: int, max_bytes: int | None = DEFAULT_MAX_TABLE_BYTES):
    """
    Build a multiplication table.

    kind:
        "standard"
        "split"
        "dual"
        "dual_split"

    n:
        Dimension exponent (algebra dimension 2**n). Must be an integer.

    max_bytes:
        Memory budget for the final table arrays. If the estimated size
        exceeds it, ValueError is raised instead of attempting the build.
        Pass a larger int on a strong machine, or None to disable the check.
        Use estimate_table_bytes(kind, n) to check a size in advance, or
        multiply(..., engine="fast") to avoid tables entirely.
    """
    kind = _normalize_kind(kind)
    n = Validation.dimension(n)

    if max_bytes is not None:
        needed = estimate_table_bytes(kind, n)
        if needed > max_bytes:
            raise ValueError(
                f"build_table({kind!r}, {n}) needs about {needed / 2**20:,.1f} MiB "
                f"(limit {max_bytes / 2**20:,.1f} MiB). Pass a larger max_bytes, "
                "max_bytes=None to disable the check, or use "
                'multiply(..., engine="fast") which needs no table.'
            )

    if kind == "standard":
        return _standard_table.build(n)

    if kind == "split":
        return _split_table.build(n)

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
        Required for split, dual and dual_split.
        Not used for standard.
    """
    kind = _normalize_kind(kind)
    engine = _normalize_engine(engine)

    # dim is part of the algebra's identity for split/dual kinds (e.g. e2*e2 is
    # +e0 in split dim=2 but -e0 for every larger split algebra), so it is
    # required, never inferred. Validate before anything else.
    if kind in ("split", "dual", "dual_split") and dim is None:
        raise ValueError(f"dim is required for {kind} multiplication")
        # Validate inputs before any early return, so malformed data always
        # raises ValidationError instead of leaking TypeError/IndexError.
    allow_eps = kind in ("dual", "dual_split")
    if allow_eps:
        Validation.basis_tuple(a, allow_zero=True, allow_eps=allow_eps)
        Validation.basis_tuple(b, allow_zero=True, allow_eps=allow_eps)
    else:
        Validation.basis_tuple(a, allow_zero=False, allow_eps=False)
        Validation.basis_tuple(b, allow_zero=False, allow_eps=False)

    # Zero propagation
    if int(a[0]) == 0 or int(b[0]) == 0:
        if kind in ("dual", "dual_split"):
            return (0, 0, 0)
        return (0, 0)

    # Dual inputs may be local 3-tuples: (sign, local_index, eps)
    if kind in ("dual", "dual_split"):
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