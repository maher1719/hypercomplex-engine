try:
    from ..core import BasisNotation
except ImportError:
    try:
        from ..core.Basis_notation import BasisNotation
    except ImportError:
        from ..core.basis_notation import BasisNotation


class CDFormat:
    """
    Formatting engine for hypercomplex elements.

    This class handles all string formatting for:
        - Single basis elements (from holo / O1 multipliers)
        - Table cells (from table builders)
        - Row/column labels

    It does NOT handle printing or file IO.
    """

    VALID_MODES = {
        "integer",
        "graded",
        "latex",
        "latex_integer",
        "latex_graded",
    }

    # ==================================================================
    # MODE HELPERS
    # ==================================================================

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

    # ==================================================================
    # CORE FORMATTING
    # ==================================================================

    @classmethod
    def format_basis(cls, index: int, mode: str) -> str:
        """
        Format only the basis element.

        integer:        e0, e1, e2, ...
        graded:         1, o1, o2, o12, ...
        latex_integer:  e_{0}, e_{1}, ...
        latex_graded:   1, o_{1}, o_{2}, o_{12}, ...
        """
        cls.validate_mode(mode)
        index = int(index)

        if mode == "integer":
            return f"e{index}"

        if mode == "graded":
            return BasisNotation.to_graded_str(index)

        if mode == "latex":
            return BasisNotation.to_latex(index, mode="integer")

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
        Format a single element from its components.

        This is the low-level formatter used by both
        format_element and the table printer.
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

    # ==================================================================
    # SINGLE ELEMENT FORMATTING (for holo / O1 output)
    # ==================================================================

    @classmethod
    def format_element(cls, element: tuple, mode: str = "integer") -> str:
        """
        Format a basis element tuple returned by holo or O1 multipliers.

        Parameters
        ----------
        element : tuple
            (sign, index)        for standard / split
            (sign, index, eps)   for dual

        mode : str
            integer, graded, latex, latex_integer, latex_graded

        Returns
        -------
        str
            Formatted string.

        Examples
        --------
        >>> CDFormat.format_element((1, 3), mode="integer")
        '+e3'

        >>> CDFormat.format_element((-1, 5), mode="graded")
        '-o13'

        >>> CDFormat.format_element((1, 2, 1), mode="integer")
        '+e2*eps'

        >>> CDFormat.format_element((0, 0, 1), mode="integer")
        '0'
        """
        if not isinstance(element, tuple):
            raise TypeError(f"element must be a tuple, got {type(element).__name__}")

        if len(element) == 2:
            sign, index = element
            return cls.format_entry(sign, index, eps=0, mode=mode)

        if len(element) == 3:
            sign, index, eps = element
            return cls.format_entry(sign, index, eps=eps, mode=mode)

        raise ValueError(
            f"element must be a 2-tuple or 3-tuple, got {len(element)}-tuple"
        )

    # ==================================================================
    # TABLE LABELS
    # ==================================================================

    @classmethod
    def basis_labels(cls, dim: int, dual: bool, mode: str) -> list:
        """
        Build row/column labels for a multiplication table.
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