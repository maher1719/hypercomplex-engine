import numpy as np
from ..validation import Validation

class VectorizedStandard:
    """
    SIMD-Accelerated O(1) Vectorized Multiplier for Standard Cayley-Dickson Algebras.
    Replaces nested Python loops with C-Level Twisted Convolution via Scatter-Add.
    """

    @staticmethod
    def _get_np_popcount():
        """Hardware-accelerated popcount for NumPy arrays."""
        try:
            return np.bitwise_count
        except AttributeError:
            def popcount_vec(x):
                x = np.asarray(x, dtype=np.uint64)
                x = x - ((x >> 1) & np.uint64(0x5555555555555555))
                x = (x & np.uint64(0x3333333333333333)) + ((x >> 2) & np.uint64(0x3333333333333333))
                x = (x + (x >> 4)) & np.uint64(0x0F0F0F0F0F0F0F0F)
                return ((x * np.uint64(0x0101010101010101)) >> np.uint64(56)).astype(np.uint8)
            return popcount_vec

    def multiply_indices_vec(self, i_arr: np.ndarray, j_arr: np.ndarray) -> tuple:
        """
        Vectorized Core O(1) multiplication for arrays of indices.
        Returns: (signs_array, indices_array)
        """
        _popcount = self._get_np_popcount()
        i = np.asarray(i_arr, dtype=np.uint64)
        j = np.asarray(j_arr, dtype=np.uint64)
        
        # Target index is always XOR
        idx = i ^ j
        
        # Masks for structural cases
        mask_zero = (i == 0) | (j == 0)
        mask_diag = (i == j) & ~mask_zero
        
        # Safe guards for nu2 (avoid log2(0) or popcount(-1) issues)
        i_safe = np.where(i == 0, 1, i)
        j_safe = np.where(j == 0, 1, j)
        ij_safe = np.where((i ^ j) == 0, 1, i ^ j)
        
        # Vectorized nu2: popcount((x & -x) - 1)
        t  = _popcount((ij_safe & -ij_safe) - np.uint64(1)).astype(np.uint64)
        ti = _popcount((i_safe & -i_safe) - np.uint64(1)).astype(np.uint64)
        tj = _popcount((j_safe & -j_safe) - np.uint64(1)).astype(np.uint64)
        
        # Branch 1: Diagonal break (t >= ti and t >= tj)
        cond1 = (t >= ti) & (t >= tj)
        s1 = 2 * ((i >> t) & np.uint64(1)).astype(np.int8) - 1
        
        # Branch 2: i-axis break (ti > t and ti >= tj)
        cond2 = (~cond1) & (ti >= tj)
        b = (j >> ti) & np.uint64(1)
        j_loc = j & ((np.uint64(1) << ti) - np.uint64(1))
        delta = np.where(j_loc == 0, 1, -1).astype(np.int8)
        s2 = np.where(b == 0, delta, -delta)
        
        # Branch 3: j-axis break (tj > t and tj > ti)
        cond3 = ~cond1 & ~cond2
        s3 = 1 - 2 * ((i >> tj) & np.uint64(1)).astype(np.int8)
        
        # Combine branches
        s = np.where(cond1, s1, np.where(cond2, s2, s3))
        
        # Higher-bit parity correction
        k = np.maximum(np.maximum(t, ti), tj) + np.uint64(1)
        parity = _popcount((i | j) >> k) & np.uint8(1)
        s = np.where(parity == 1, -s, s)
        
        # Apply structural overrides
        s = np.where(mask_zero, 1, s)
        s = np.where(mask_diag, -1, s)
        
        return s.astype(np.int8), idx

    def multiply_arrays(self, arr1: np.ndarray, arr2: np.ndarray) -> np.ndarray:
        """
        Full polynomial multiplication using Twisted Convolution via C-Level Scatter-Add.
        Automatically handles sparsity to prevent RAM explosion on high dimensions.
        """
        size1 = len(arr1)
        size2 = len(arr2)
        max_size = max(size1, size2)
        
        if max_size <= 1:
            return np.array([arr1[0] * arr2[0]], dtype=np.float64) if size1 > 0 and size2 > 0 else np.zeros(1, dtype=np.float64)

        # Calculate required power-of-two dimension
        # Equivalent to 1 << math.ceil(math.log2(max_size))
        dim = 1 << (max_size - 1).bit_length()

        # 1. SPARSE OPTIMIZATION: Isolate non-zero elements
        nz_A = np.nonzero(arr1)[0]
        nz_B = np.nonzero(arr2)[0]
        
        if len(nz_A) == 0 or len(nz_B) == 0:
            return np.zeros(dim, dtype=np.float64)

        # 2. Create the 2D Grid of ONLY the active indices
        I, J = np.meshgrid(nz_A, nz_B, indexing='ij')
        
        # 3. Calculate the entire Sign and Index Matrix in one C-level SIMD sweep
        S, V = self.multiply_indices_vec(I, J)
        
        # 4. Calculate the Outer Product of ONLY the active coefficients
        M = np.outer(arr1[nz_A], arr2[nz_B])
        
        # 5. Apply Signs to the Coefficients
        weighted_M = S.astype(np.float64) * M
        
        # 6. C-Level Scatter-Add into the Result Array
        result = np.zeros(dim, dtype=np.float64)
        np.add.at(result, V.ravel(), weighted_M.ravel())
        
        return result