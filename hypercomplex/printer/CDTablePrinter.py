import numpy as np
from ..core.Basis_notation import BasisNotation




# =============================================================================
# .PRESENTATION ENGINE (The Formatting Layer)
# =============================================================================
class CDTablePrinter:
    """
    Handles terminal printing and CSV export for hypercomplex multiplication tables.
    Completely decoupled from CDTableBuilder. Can also be used to format single 
    O(1) or O(n) multiplier outputs.
    """
    
    @staticmethod
    def format_entry(sign: int, k: int, eps: int = 0, mode: str = "integer") -> str:
        """
        Formats a single table cell or O(1) product.
        mode: "integer" (+e5), "graded" (+e13), or "latex" (+e_{13})
        """
        if sign == 0:
            return "0"
            
        prefix = "+" if sign > 0 else "-"
        
        if mode == "integer":
            base = f"e{k}"
        elif mode == "graded":
            base = BasisNotation.to_graded_str(k)
        elif mode == "latex":
            base = BasisNotation.to_latex(k)
        else:
            raise ValueError("mode must be 'integer', 'graded', or 'latex'")
            
        if eps:
            if mode == "latex":
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
        
        if mode == "integer":
            labels = [f"e{i}" for i in range(base_dim)]
        elif mode == "graded":
            labels = [BasisNotation.to_graded_str(i) for i in range(base_dim)]
        elif mode == "latex":
            labels = [BasisNotation.to_latex(i) for i in range(base_dim)]
        else:
            raise ValueError("Invalid mode")
            
        if dual:
            eps_str = "\\epsilon" if mode == "latex" else "eps"
            labels.append(eps_str)
            
            for k in range(1, base_dim):
                if mode == "integer":
                    labels.append(f"e{k}*{eps_str}")
                elif mode == "graded":
                    base = BasisNotation.to_graded_str(k)
                    labels.append(f"{base}*{eps_str}")
                elif mode == "latex":
                    base = BasisNotation.to_latex(k)
                    labels.append(f"{base}{eps_str}")
                    
        return labels

    @staticmethod
    def print_table(signs, indices, eps=None, title="Table",mode="integer", limit=None ):
        """Prints a formatted table to the terminal."""
        dim = signs.shape[0]
        lim = dim if limit is None else min(limit, dim)
        dual = eps is not None
        
        labels = CDTablePrinter.basis_labels(dim, dual=dual, mode=mode)[:lim]
        
        rows = []
        for i in range(lim):
            row = []
            for j in range(lim):
                e = 0 if eps is None else int(eps[i, j])
                row.append(CDTablePrinter.format_entry(int(signs[i, j]), int(indices[i, j]), e, mode))
            rows.append(row)
            
        print(title)
        if not rows: return
        
        cell_width = max(len(c) for r in rows for c in r)
        cell_width = max(cell_width, max(len(l) for l in labels))
        label_width = max(len(l) for l in labels)
        
        # Header
        print(" " * (label_width + 3) + " ".join(f"{l:>{cell_width}}" for l in labels))
        
        # Rows
        for lab, r in zip(labels, rows):
            print(f"{lab:>{label_width}} | " + " ".join(f"{c:>{cell_width}}" for c in r))
        print()

    @staticmethod
    def export_csv(path, signs, indices, eps=None, mode="integer", csv_mode="matrix"):
        """
        Exports table to CSV.
        csv_mode="matrix": spreadsheet grid with formatted labels.
        csv_mode="long":   raw data rows (i,j,sign,index[,eps]) for analysis.
        """
        dim = signs.shape[0]
        dual = eps is not None
        
        with open(path, "w", newline="") as f:
            if csv_mode == "matrix":
                labels = CDTablePrinter.basis_labels(dim, dual=dual, mode=mode)
                f.write("," + ",".join(labels) + "\n")
                
                for i in range(dim):
                    cells = []
                    for j in range(dim):
                        e = 0 if eps is None else int(eps[i, j])
                        cells.append(CDTablePrinter.format_entry(int(signs[i, j]), int(indices[i, j]), e, mode))
                    f.write(labels[i] + "," + ",".join(cells) + "\n")
                    
            elif csv_mode == "long":
                header = "i,j,sign,index"
                if dual: header += ",eps"
                f.write(header + "\n")
                
                for i in range(dim):
                    for j in range(dim):
                        line = f"{i},{j},{int(signs[i, j])},{int(indices[i, j])}"
                        if dual: line += f",{int(eps[i, j])}"
                        f.write(line + "\n")
            else:
                raise ValueError("csv_mode must be 'matrix' or 'long'")