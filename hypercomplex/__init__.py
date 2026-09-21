import logging

# Library logging etiquette: attach a NullHandler and let the APPLICATION decide
# whether and how to display messages (never call logging.basicConfig here).
logging.getLogger(__name__).addHandler(logging.NullHandler())

from .facade import (
    build_table,
    estimate_table_bytes,
    multiply,
    format_element,
    print_table,
    export_csv,
)

from .core import (
    BasisElement,
    BasisNotation,
    Validation,

    StandardTableBuilder,
    SplitTableBuilder,
    DualTableBuilder,

    StandardHolographic,
    SplitHolographic,
    DualHolographic,

    FastStandard,
    FastSplit,
    FastDual,
)

from .printer import (
    CDFormat,
    CDTablePrinter,
)

__all__ = [
    # Simple API
    "build_table",
    "estimate_table_bytes",
    "multiply",
    "format_element",
    "print_table",
    "export_csv",

    # Core
    "BasisElement",
    "BasisNotation",
    "Validation",

    # Table builders
    "StandardTableBuilder",
    "SplitTableBuilder",
    "DualTableBuilder",

    # Holographic O(n)
    "StandardHolographic",
    "SplitHolographic",
    "DualHolographic",

    # O(1)
    "FastStandard",
    "FastSplit",
    "FastDual",

    # Printer / formatter
    "CDFormat",
    "CDTablePrinter",
]