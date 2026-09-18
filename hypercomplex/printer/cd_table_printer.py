from ..core.basis_notation import BasisNotation



class CDTablePrinter:
    """
    Printer / CSV exporter for hypercomplex multiplication tables.

    This class is presentation-only.
    It does not contain algebra rules.
    It delegates notation to BasisNotation.

    Supported modes:
        "integer"        -> e5
        "graded"         -> o13
        "latex"          -> o_{13}   alias for latex_graded
        "latex_integer"  -> e_{5}
        "latex_graded"   -> o_{13}

    Dual tables:
        Pass eps array as third argument.
        eps[i, j] = 1 means the product carries epsilon.
    """

    _VALID_MODES = {
        "integer",
        "graded",
        "latex",
        "latex_integer",
        "latex_graded",
    }

    # ==================================================================
    # INTERNAL HELPERS
    # ==================================================================

    @classmethod
    def _validate_mode(cls, mode: str) -> None:
        if mode not in cls._VALID_MODES:
            raise ValueError(
                f"Invalid mode '{mode}'. "
                f"Valid modes: {sorted(cls._VALID_MODES)}"
            )

    @staticmethod
    def _is_latex(mode: str) -> bool:
        return mode.startswith("latex")

    @classmethod
    def _format_basis(cls, index: int, mode: str) -> str:
        """
        Formats only the basis element.

        integer:
            e0, e1, e2, ...

        graded:
            1, o1, o2, o12, ...

        latex_integer:
            e_{0}, e_{1}, ...

        latex_graded / latex:
            1, o_{1}, o_{2}, o_{12}, ...
        """
        cls._validate_mode(mode)
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
    def _format_entry(
        cls,
        sign: int,
        index: int,
        eps: int = 0,
        mode: str = "integer",
    ) -> str:
        """
        Formats one table cell.

        Standard / split:
            +e3, -e5, 0

        Dual:
            +e3*eps, -eps, 0

        LaTeX dual:
            +e_{3}\epsilon, -\epsilon, 0
        """
        cls._validate_mode(mode)

        sign = int(sign)
        index = int(index)
        eps = int(eps)

        if sign == 0:
            return "0"

        prefix = "+" if sign > 0 else "-"
        base = cls._format_basis(index, mode)

        if eps == 0:
            return prefix + base

        latex = cls._is_latex(mode)
        eps_label = "\\epsilon" if latex else "eps"

        # epsilon alone: +eps or +\epsilon
        if index == 0:
            return prefix + eps_label

        if latex:
            return f"{prefix}{base}{eps_label}"

        return f"{prefix}{base}*{eps_label}"

    @classmethod
    def _basis_labels(cls, dim: int, dual: bool, mode: str) -> list:
        """
        Builds row/column labels.

        For standard/split:
            dim labels

        For dual:
            base labels + epsilon labels
        """
        cls._validate_mode(mode)

        base_dim = dim // 2 if dual else dim

        labels = [
            cls._format_basis(i, mode)
            for i in range(base_dim)
        ]

        if not dual:
            return labels

        latex = cls._is_latex(mode)
        eps_label = "\\epsilon" if latex else "eps"

        # epsilon itself: eps*e0
        labels.append(eps_label)

        # epsilon times non-scalar basis elements
        for i in range(1, base_dim):
            base = cls._format_basis(i, mode)

            if latex:
                labels.append(f"{base}{eps_label}")
            else:
                labels.append(f"{base}*{eps_label}")

        return labels

    # ==================================================================
    # TERMINAL PRINTING
    # ==================================================================

    @classmethod
    def print_table(
        cls,
        signs,
        indices,
        eps=None,
        title: str | None = None,
        limit: int | None = None,
        mode: str = "integer",
    ):
        """
        Print a multiplication table to the terminal.

        Parameters
        ----------
        signs:
            2D array of signs.

        indices:
            2D array of basis indices.

        eps:
            Optional 2D array for dual tables.
            If provided, the table is treated as dual.

        title:
            Optional title.

        limit:
            Optional maximum number of rows/columns to display.

        mode:
            integer, graded, latex, latex_integer, latex_graded
        """
        cls._validate_mode(mode)

        dim = signs.shape[0]
        lim = dim if limit is None else min(limit, dim)

        dual = eps is not None

        labels = cls._basis_labels(dim, dual=dual, mode=mode)[:lim]

        rows = []

        for i in range(lim):
            row = []

            for j in range(lim):
                e = 0 if eps is None else int(eps[i, j])

                row.append(
                    cls._format_entry(
                        int(signs[i, j]),
                        int(indices[i, j]),
                        e,
                        mode,
                    )
                )

            rows.append(row)

        if title is None:
            title = "Dual Table" if dual else "Table"

        print(title)

        if not rows:
            return rows

        cell_width = max(len(cell) for row in rows for cell in row)
        label_width = max(len(label) for label in labels)

        # Header
        print(
            " " * (label_width + 3)
            + " ".join(f"{label:>{cell_width}}" for label in labels)
        )

        # Rows
        for label, row in zip(labels, rows):
            print(
                f"{label:>{label_width}} | "
                + " ".join(f"{cell:>{cell_width}}" for cell in row)
            )

        print()

        return rows

    # ==================================================================
    # CSV EXPORT
    # ==================================================================

    @classmethod
    def export_csv(
        cls,
        path: str,
        signs,
        indices,
        eps=None,
        mode: str = "integer",
        csv_mode: str = "matrix",
    ) -> str:
        """
        Export table to CSV.

        Parameters
        ----------
        path:
            Output file path.

        signs:
            2D sign array.

        indices:
            2D index array.

        eps:
            Optional epsilon array for dual tables.

        mode:
            integer, graded, latex, latex_integer, latex_graded

        csv_mode:
            "matrix":
                Spreadsheet-style grid.

            "long":
                One row per product:
                    i,j,sign,index[,eps]
        """
        cls._validate_mode(mode)

        if csv_mode not in ("matrix", "long"):
            raise ValueError("csv_mode must be 'matrix' or 'long'")

        dim = signs.shape[0]
        dual = eps is not None

        with open(path, "w", newline="", encoding="utf-8") as f:

            if csv_mode == "matrix":
                labels = cls._basis_labels(dim, dual=dual, mode=mode)

                # Header row
                f.write("," + ",".join(labels) + "\n")

                for i in range(dim):
                    cells = []

                    for j in range(dim):
                        e = 0 if eps is None else int(eps[i, j])

                        cells.append(
                            cls._format_entry(
                                int(signs[i, j]),
                                int(indices[i, j]),
                                e,
                                mode,
                            )
                        )

                    f.write(labels[i] + "," + ",".join(cells) + "\n")

            else:
                header = "i,j,sign,index"

                if dual:
                    header += ",eps"

                f.write(header + "\n")

                for i in range(dim):
                    for j in range(dim):
                        line = (
                            f"{i},"
                            f"{j},"
                            f"{int(signs[i, j])},"
                            f"{int(indices[i, j])}"
                        )

                        if dual:
                            line += f",{int(eps[i, j])}"

                        f.write(line + "\n")

        return path