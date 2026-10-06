import numpy as np


def index_dtype(dim: int):
    """
    Choose the smallest unsigned dtype capable of storing indices 0 ... dim-1.
    """
    if dim <= 1 << 8:
        return np.uint8
    if dim <= 1 << 16:
        return np.uint16
    if dim <= 1 << 32:
        return np.uint32
    return np.uint64


def real_table():
    """
    The real numbers A_0.
    """
    signs = np.array([[1]], dtype=np.int8)
    indices = np.array([[0]], dtype=np.uint8)
    return signs, indices