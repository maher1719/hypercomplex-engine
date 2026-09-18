# tests/test_holographic_vs_table.py
import pytest
import numpy as np
from hypercomplex import CDTableBuilder
from hypercomplex.core.holographic_multiplier import HolographicMultiplier

def test_standard_cross_validation():
    """Proves the O(n) multiplier perfectly matches the O(4^n) table builder."""
    n = 5  # Sedenions (16x16 = 256 entries)
    signs_table, indices_table = CDTableBuilder.standard(n)
    
    holo = HolographicMultiplier()
    
    dim = 1 << n
    for i in range(dim):
        for j in range(dim):
            h_sign, h_idx = holo.standard((1,i), (1,j))
            assert signs_table[i, j] == h_sign, f"Sign mismatch at {i},{j}"
            assert indices_table[i, j] == h_idx, f"Index mismatch at {i},{j}"

def test_split_cross_validation():
    """Proves the Split O(n) multiplier perfectly matches the Split table builder."""
    n = 5  
    signs_table, indices_table = CDTableBuilder.split(n)
    
    holo = HolographicMultiplier()
    
    dim = 1 << n
    for i in range(dim):
        for j in range(dim):
            h_sign, h_idx = holo.split((1,i), (1,j), dim=n)
            assert signs_table[i, j] == h_sign, f"Split sign mismatch at {i},{j}"
            assert indices_table[i, j] == h_idx, f"Split index mismatch at {i},{j}"
def further_validation():
    holo = HolographicMultiplier()

    # --- standard: quaternions ---
    assert holo.standard((1, 1), (1, 2)) == (1, 3)   # e1·e2 = +e3
    assert holo.standard((1, 2), (1, 1)) == (-1, 3)  # e2·e1 = -e3
    assert holo.standard((1, 1), (1, 1)) == (-1, 0)  # e1² = -e0

    # --- split: split-complex ---
    assert holo.split((1, 1), (1, 1), dim=1) == (1, 0)   # e1² = +e0

    # --- split: split-octonions diagonal ---
    assert holo.split((1, 4), (1, 4), dim=3) == (1, 0)   # new element squares to +1
    assert holo.split((1, 1), (1, 1), dim=3) == (-1, 0)  # old element squares to -1

    # --- dual: nilpotency ---
    assert holo.dual((1, 1), (1, 1), dim=0) == (0, 0, 1)  # ε·ε = 0

    # --- dual: base product preserved ---
    s, idx, eps = holo.dual((1, 0), (1, 1), dim=1)
    assert eps == 0  # no ε involved

    # --- dual: one ε factor ---
    s, idx, eps = holo.dual((1, 0), (1, 3), dim=1)  # e0 · (εe1)
    assert eps == 1

    print("All HolographicMultiplier tuple tests passed.")