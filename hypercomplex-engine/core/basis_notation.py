# basis_notation.py

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
            # Extract digits from strings like "o13" or "o_{145}"
            # NOTE: This assumes single-digit generators (valid up to 512-dim / 9 generators).
            gens = [int(char) for char in gens if char.isdigit()]
        k = 0
        for g in gens:
            k |= (1 << (g - 1))
        return k

    @staticmethod
    def to_graded_str(k: int) -> str:
        """Integer to graded string (e.g., 5 -> 'o13', 0 -> '1')."""
        if k == 0: return "1"
        gens = BasisNotation.int_to_generators(k)
        return "o" + "".join(str(g) for g in gens)

    @staticmethod
    def to_latex(k: int, mode: str = "integer") -> str:
        """Integer to LaTeX string.
        mode="integer": e_{5}, e_{13}   (index as number)
        mode="graded":  o_{13}, o_{123} (generators)
        """
        if k == 0:
            return "1"

        if mode == "integer":
            return f"e_{{{k}}}"
        elif mode == "graded":
            gens = BasisNotation.int_to_generators(k)
            return "o_{" + "".join(str(g) for g in gens) + "}"
        else:
            raise ValueError("mode must be 'integer' or 'graded'")

    @staticmethod
    def from_graded_str(s: str) -> int:
        """Parses a graded string back to an integer index.
        Accepts 'o13', 'o_{13}', '1', and also legacy 'e13'.
        """
        s = s.strip()
        
        # Zero / scalar cases
        if s in ("1", "o", "o0", "e", "e0", "o_\\emptyset"):
            return 0
        
        # Strip prefixes and braces; digit extraction does the real work
        s = (s.replace("o_", "").replace("o", "")
              .replace("e_", "").replace("e", "")
              .replace("{", "").replace("}", ""))
        
        gens = [int(char) for char in s if char.isdigit()]
        return BasisNotation.generators_to_int(gens)

    @staticmethod
    def format_entry(sign: int, k: int, eps: int = 0, mode: str = "integer") -> str:
        """Formats a single table cell or O(1) product."""
        if sign == 0:
            return "0"

        prefix = "+" if sign > 0 else "-"

        if mode == "integer":
            base = f"e{k}"
        elif mode == "graded":
            base = BasisNotation.to_graded_str(k)
        elif mode == "latex":              # alias for graded-latex (backward compat)
            base = BasisNotation.to_latex(k, mode="graded")
        elif mode == "latex_integer":
            base = BasisNotation.to_latex(k, mode="integer")
        elif mode == "latex_graded":
            base = BasisNotation.to_latex(k, mode="graded")
        else:
            raise ValueError("mode must be 'integer', 'graded', 'latex', 'latex_integer', or 'latex_graded'")

        if eps:
            if mode.startswith("latex"):
                eps_str = "\\epsilon"
                if k == 0: return f"{prefix}{eps_str}"
                return f"{prefix}{base}{eps_str}"
            else:
                eps_str = "eps"
                if k == 0: return f"{prefix}{eps_str}"
                return f"{prefix}{base}*{eps_str}"

        return prefix + base

    @staticmethod
    def basis_labels(dim: int, dual: bool = False, mode: str = "integer") -> list:
        """Generates the headers/row labels for the table."""
        base_dim = dim // 2 if dual else dim

        def _base_label(i):
            if mode == "integer": return f"e{i}"
            elif mode == "graded": return BasisNotation.to_graded_str(i)
            elif mode == "latex": return BasisNotation.to_latex(i, mode="graded")
            elif mode == "latex_integer": return BasisNotation.to_latex(i, mode="integer")
            elif mode == "latex_graded": return BasisNotation.to_latex(i, mode="graded")
            else: raise ValueError("Invalid mode")

        labels = [_base_label(i) for i in range(base_dim)]

        if dual:
            is_latex = mode.startswith("latex")
            eps_str = "\\epsilon" if is_latex else "eps"
            labels.append(eps_str)

            for k in range(1, base_dim):
                base = _base_label(k)
                labels.append(f"{base}{eps_str}" if is_latex else f"{base}*{eps_str}")

        return labels