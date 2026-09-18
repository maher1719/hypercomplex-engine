class CDTableBuilder:
    
    @staticmethod
    def standard(n):
        """
        Generates a standard, ordinary Cayley-Dickson table of dimension 2^n.
        Used to populate the immutable sub-algebra background for the lower half.
        """
        table = [[(1, 0)]]
        for _ in range(1, n + 1):
            half_dim = len(table)
            full_dim = half_dim * 2
            new_table = [[(1, 0) for _ in range(full_dim)] for _ in range(full_dim)]
            for i in range(full_dim):
                for j in range(full_dim):
                    is_row_high = (i >= half_dim)
                    is_col_high = (j >= half_dim)
                    local_i, local_j = i % half_dim, j % half_dim
                    
                    if not is_row_high and not is_col_high:
                        sign, idx = table[local_i][local_j]
                        new_table[i][j] = (sign, idx)
                    elif not is_row_high and is_col_high:
                        sign, idx = table[local_j][local_i]
                        new_table[i][j] = (sign, idx + half_dim)
                    elif is_row_high and not is_col_high:
                        sign, idx = table[local_i][local_j]
                        if local_j > 0: sign = -sign
                        new_table[i][j] = (sign, idx + half_dim)
                    else:
                        # Standard CD slot uses negative -d*b
                        sign, idx = table[local_j][local_i]
                        if local_j > 0: sign = -sign
                        new_table[i][j] = (-sign, idx)
            table = new_table
        return table
    @staticmethod
    def split(n):
        """
        Iteratively builds the hybrid split table.
        Populates A_(n-1) with ordinary hypercomplex math,
        and applies the split formula for extended higher dimensions.
        """
        if n == 0:
            return [[(1, 0)]]
            
        full_dim = 1 << n
        half_dim = 1 << (n - 1)
        
        # 1. Pre-compute the pure standard hypercomplex parent table of level A_(n-1)
        standard_parent = CDTableBuilder.standard(n - 1)
        
        # Initialize the complete full-dimensional matrix grid
        table = [[(1, 0) for _ in range(full_dim)] for _ in range(full_dim)]
        
        # 2. Populate the entire grid using your hybrid block rules
        for i in range(full_dim):
            for j in range(full_dim):
                is_row_high = (i >= half_dim)
                is_col_high = (j >= half_dim)
                
                local_i = i % half_dim
                local_j = j % half_dim
                
                if not is_row_high and not is_col_high:
                    # Quadrant a: Pure Ordinary Hypercomplex Subalgebra inheritance
                    sign, idx = standard_parent[local_i][local_j]
                    table[i][j] = (sign, idx)
                    
                elif not is_row_high and is_col_high:
                    # Quadrant b: (a, 0)(0, d) = (0, da) -> Uses standard parent states
                    sign, idx = standard_parent[local_j][local_i]
                    table[i][j] = (sign, idx + half_dim)
                    
                elif is_row_high and not is_col_high:
                    # Quadrant c: (0, b)(c, 0) = (0, bc*) -> Standard parent with right conjugation
                    sign, idx = standard_parent[local_i][local_j]
                    if local_j > 0:
                        sign = -sign
                    table[i][j] = (sign, idx + half_dim)
                    
                else:
                    # Quadrant d: 's Split Space -> (0, b)(0, d) = (+d*b, 0)
                    # Directly mirrors the parent signs without standard CD's global inversion
                    sign, idx = standard_parent[local_j][local_i]
                    if local_j > 0:
                        sign = -sign
                    table[i][j] = (sign, idx)

        # 3. Apply Global Positivity Constraints on First Row and Column Bounds
        for i in range(full_dim):
            # Enforce absolute positive matrix edges for e_0 interactions
            table[i][0] = (1, i)
            table[0][i] = (1, i) 
        return table

    @staticmethod
    def dual(n, split=False):
        base = CDTableBuilder.split(n) if split else CDTableBuilder.standard(n)
        N = 1 << n
        full = 2 * N
        table = [[(0, 0, 0) for _ in range(full)] for _ in range(full)]
        for i in range(full):
            for j in range(full):
                i_eps, j_eps = i >= N, j >= N
                i_loc, j_loc = i % N, j % N
                if i_eps and j_eps:
                    table[i][j] = (0, 0, 1)             # eps^2 = 0
                elif i_eps or j_eps:
                    
                    sign, idx = base[i_loc][j_loc]
                    table[i][j] = (sign, idx, 1)         # a*eps or eps*a
                else:
                    
                    sign, idx = base[i][j]
                    table[i][j] = (sign, idx, 0)
        return table
    def print_standard(self,n):
        table=self.standard(n)
        formatted_table = []
        for r in table:
            row_str = []
            for sign, idx in r:
                s = "+" if sign == 1 else "-"
                row_str.append(f"{s}e{idx}")
            formatted_table.append(row_str)
        return formatted_table
