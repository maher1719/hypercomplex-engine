# cd_table_printer.py
from ..core.Basis_notation import BasisNotation

class CDTablePrinter:
    """
    Handles terminal printing and CSV export. 
    Delegates all string formatting to BasisNotation.
    """
    
    @staticmethod
    def print_table(signs, indices, eps=None, title="Table", limit=None, mode="integer"):
        dim = signs.shape[0]
        lim = dim if limit is None else min(limit, dim)
        dual = eps is not None
        
        # 1. Get labels from BasisNotation
        labels = BasisNotation.basis_labels(dim, dual=dual, mode=mode)[:lim]
        
        # 2. Build rows using BasisNotation
        rows = []
        for i in range(lim):
            row = []
            for j in range(lim):
                e = 0 if eps is None else int(eps[i, j])
                row.append(BasisNotation.format_entry(int(signs[i, j]), int(indices[i, j]), e, mode))
            rows.append(row)
            
        # 3. Print layout
        print(title)
        if not rows: return
        
        cell_width = max(len(c) for r in rows for c in r)
        cell_width = max(cell_width, max(len(l) for l in labels))
        label_width = max(len(l) for l in labels)
        
        print(" " * (label_width + 3) + " ".join(f"{l:>{cell_width}}" for l in labels))
        for lab, r in zip(labels, rows):
            print(f"{lab:>{label_width}} | " + " ".join(f"{c:>{cell_width}}" for c in r))
        print()

    @staticmethod
    def export_csv(path, signs, indices, eps=None, mode="integer", csv_mode="matrix"):
        dim = signs.shape[0]
        dual = eps is not None
        
        with open(path, "w", newline="") as f:
            if csv_mode == "matrix":
                labels = BasisNotation.basis_labels(dim, dual=dual, mode=mode)
                f.write("," + ",".join(labels) + "\n")
                
                for i in range(dim):
                    cells = []
                    for j in range(dim):
                        e = 0 if eps is None else int(eps[i, j])
                        cells.append(BasisNotation.format_entry(int(signs[i, j]), int(indices[i, j]), e, mode))
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
        print(f"File succufuly generated with path: {path}")