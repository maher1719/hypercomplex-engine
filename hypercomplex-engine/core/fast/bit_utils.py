def nu2(x: int) -> int:
    """
    2-adic valuation: index of the lowest set bit.

    nu2(1) = 0
    nu2(2) = 1
    nu2(4) = 2
    nu2(12) = 2
    """
    if x == 0:
        raise ValueError("nu2(0) is undefined in the O(1) structural descent")

    return (x & -x).bit_length() - 1


def popcount(x: int) -> int:
    """
    Number of set bits.
    """
    return int(x).bit_count()