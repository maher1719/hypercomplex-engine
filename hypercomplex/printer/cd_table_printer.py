from .cd_format import CDFormat


class CDTablePrinter:
    """
    Printer / CSV exporter for hypercomplex multiplication tables.

    This class handles only:
        - terminal printing
        - CSV export

    All formatting is delegated to CDTableFormat.
    """

    @staticmethod
    def print_table(
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
        CDTableFormat.validate_mode(mode)

        dim = signs.shape[0]
        lim = dim if limit is None else min(limit, dim)

        if lim == 0:
            return []

        dual = eps is not None

        labels = CDTableFormat.basis_labels(
            dim,
            dual=dual,
            mode=mode,
        )[:lim]

        rows = []

        for i in range(lim):
            row = []

            for j in range(lim):
                e = 0 if eps is None else int(eps[i, j])

                row.append(
                    CDTableFormat.format_entry(
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

    @staticmethod
    def export_csv(
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
        CDTableFormat.validate_mode(mode)

        if csv_mode not in ("matrix", "long"):
            raise ValueError("csv_mode must be 'matrix' or 'long'")

        dim = signs.shape[0]
        dual = eps is not None

        with open(path, "w", newline="", encoding="utf-8") as f:

            if csv_mode == "matrix":
                labels = CDTableFormat.basis_labels(
                    dim,
                    dual=dual,
                    mode=mode,
                )

                # Header row
                f.write("," + ",".join(labels) + "\n")

                for i in range(dim):
                    cells = []

                    for j in range(dim):
                        e = 0 if eps is None else int(eps[i, j])

                        cells.append(
                            CDTableFormat.format_entry(
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