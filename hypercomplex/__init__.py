"""
Hypercomplex: table generator 
for ordinary, split, dual, and Cayley-Dickson algebras.
"""
from .core.CDTableBuilder import CDTableBuilder
from .core.Basis_notation import BasisNotation
from .printer.CDTablePrinter import CDTablePrinter
from .core.holographic_multiplier import HolographicMultiplier

__version__ = "0.1.0"
__all__ = ["CDTableBuilder", "BasisNotation", "CDTablePrinter","HolographicMultiplier"]