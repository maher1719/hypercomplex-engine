# tests/test_algebra.py
import pytest
import numpy as np
from hypercomplex import CDTableBuilder

def test_quaternion_anticommutativity():
    signs, indices = CDTableBuilder.standard(2)
    # i * j = k  => signs[1,2] == 1, indices[1,2] == 3
    assert signs[1, 2] == 1
    assert indices[1, 2] == 3
    
    # j * i = -k => signs[2,1] == -1, indices[2,1] == 3
    assert signs[2, 1] == -1
    assert indices[2, 1] == 3

def test_split_complex_square():
    signs, indices = CDTableBuilder.split(1)
    # In split-complex, e1^2 = +e0
    assert signs[1, 1] == 1
    assert indices[1, 1] == 0

def test_dual_nilpotency():
    signs, indices, eps = CDTableBuilder.dual(0)
    # epsilon^2 = 0
    assert signs[1, 1] == 0