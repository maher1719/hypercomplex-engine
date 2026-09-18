import numpy as np

# =============================================================================
# 1. BASIS NOTATION ENGINE (The Translation Layer)
# =============================================================================
class BasisNotation:
    """
    Handles conversions between integer indices, graded (bitmask) notation, 
    and LaTeX strings for Cayley-Dickson basis elements.
    """
    
    @staticmethod
    def int_to_generators(k: int) -> tuple:
        """Converts integer index k to a tuple of 1-based generator indices."""
        if k == 0: return ()
        gens = []
        i = 1
        temp_k = k
        while temp_k > 0:
            if temp_k & 1: gens.append(i)
            temp_k >>= 1
            i += 1
        return tuple(gens)

    @staticmethod
    def generators_to_int(gens) -> int:
        """Converts generators (tuple, list, or string like '13') to integer index."""
        if isinstance(gens, str):
            # Extract digits from strings like "e13" or "e_{145}"
            gens = [int(char) for char in gens if char.isdigit()]
        k = 0
        for g in gens:
            k |= (1 << (g - 1))
        return k

    @staticmethod
    def to_graded_str(k: int) -> str:
        """Integer to graded string (e.g., 5 -> 'e13', 0 -> '1')."""
        if k == 0: return "1"
        gens = BasisNotation.int_to_generators(k)
        return "o" + "".join(str(g) for g in gens)

    @staticmethod
    def to_latex(k: int) -> str:
        """Integer to LaTeX graded string (e.g., 5 -> 'e_{13}', 0 -> '1')."""
        if k == 0: return "1"
        gens = BasisNotation.int_to_generators(k)
        return "e_{" + "".join(str(g) for g in gens) + "}"

    @staticmethod
    def from_graded_str(s: str) -> int:
        """Parses a graded string back to an integer index."""
        s = s.strip()
        if s in ("1", "e0", "e", "e_\\emptyset"): return 0
        s = s.replace("e_", "").replace("e", "").replace("{", "").replace("}", "")
        gens = [int(char) for char in s if char.isdigit()]
        return BasisNotation.generators_to_int(gens)


