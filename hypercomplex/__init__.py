"""
Hypercomplex: table generator 
for ordinary, split, dual, and Cayley-Dickson algebras.
"""



from .core import (
    BasisElement,
    Basis_notation,
    Validation,
    StandardTableBuilder,
    SplitTableBuilder,
    DualTableBuilder,
    StandardHolographic,
    SplitHolographic,
    DualHolographic,
)


from .printer import CDTableFormat, CDTablePrinter


__version__ = "0.2.0"

__all__ = [
    "BasisElement",
    "Basis_notation",
    "Validation",
    "StandardTableBuilder",
    "SplitTableBuilder",
    "DualTableBuilder",
    "StandardHolographic",
    "SplitHolographic",
    "DualHolographic",
    "CDTableFormat",
    "CDTablePrinter",
]