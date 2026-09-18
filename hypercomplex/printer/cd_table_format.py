try:
    from ..core import BasisNotation
except ImportError:
    try:
        from ..core.Basis_notation import BasisNotation
    except ImportError:
        from ..core.basis_notation import BasisNotation


class CDTableFormat:
    """
    Formatting helpers for hypercomplex multiplication tables.

    This class contains no printing and no file IO.
    It only converts signs, indices, and epsilon flags into strings.

    Supported modes:
        "integer"        -> e5
        "graded"         -> o13
        "latex"          -> o_{13}   alias for latex_graded
        "latex_integer"  -> e_{5}
        "latex_graded"   -> o_{13}
    """

    VALID_MODES = {
        "integer",
        "graded",
        "latex",
        "latex_integer",
        "latex_graded",
    }

    @classmethod
    def validate_mode(cls, mode: str) -> None:
        if mode not in cls.VALID_MODES:
            raise ValueError(
                f"Invalid mode '{mode}'. "
                f"Valid modes: {sorted(cls.VALID_MODES)}"
            )

    @staticmethod
    def is_latex(mode: str) -> bool:
        return mode.startswith("latex")

    @classmethod
    def format_basis(cls, index: int, mode: str) -> str:
        """
        Format only the basis element.

        integer:
            e0, e1, e2, ...

        graded:
            1, o1, o2, o12, ...

        latex_integer:
            e_{0}, e_{1}, ...

        latex / latex_graded:
            1, o_{1}, o_{2}, o_{12}, ...
        """
        cls.validate_mode(mode)

        index = int(index)

        if mode == "integer":
            return f"e{index}"

        if mode == "graded":
            return BasisNotation.to_graded_str(index)

        if mode == "latex":
            return BasisNotation.to_latex(index, mode="graded")

        if mode == "latex_integer":
            return BasisNotation.to_latex(index, mode="integer")

        if mode == "latex_graded":
            return BasisNotation.to_latex(index, mode="graded")

        raise ValueError(f"Invalid mode '{mode}'")

    @classmethod
    def format_entry(
        cls,
        sign: int,
        index: int,
        eps: int = 0,
        mode: str = "integer",
    ) -> str:
        """
        Format one table cell.

        Standard / split:
            +e3, -e5, 0

        Dual:
            +e3*eps, -eps, 0

        LaTeX dual:
            +e_{3}\epsilon, -\epsilon, 0
        """
        cls.validate_mode(mode)

        sign = int(sign)
        index = int(index)
        eps = int(eps)

        if sign == 0:
            return "0"

        prefix = "+" if sign > 0 else "-"
        base = cls.format_basis(index, mode)

        if eps == 0:
            return prefix + base

        latex = cls.is_latex(mode)
        eps_label = "\\epsilon" if latex else "eps"

        # epsilon alone: +eps or +\epsilon
        if index == 0:
            return prefix + eps_label

        if latex:
            return f"{prefix}{base}{eps_label}"

        return f"{prefix}{base}*{eps_label}"

    @classmethod
    def basis_labels(cls, dim: int, dual: bool, mode: str) -> list:
        """
        Build row/column labels.

        For standard/split:
            dim labels

        For dual:
            base labels + epsilon labels
        """
        cls.validate_mode(mode)

        base_dim = dim // 2 if dual else dim

        labels = [
            cls.format_basis(i, mode)
            for i in range(base_dim)
        ]

        if not dual:
            return labels

        latex = cls.is_latex(mode)
        eps_label = "\\epsilon" if latex else "eps"

        # epsilon itself: eps*e0
        labels.append(eps_label)

        # epsilon times non-scalar basis elements
        for i in range(1, base_dim):
            base = cls.format_basis(i, mode)

            if latex:
                labels.append(f"{base}{eps_label}")
            else:
                labels.append(f"{base}*{eps_label}")

        return labels