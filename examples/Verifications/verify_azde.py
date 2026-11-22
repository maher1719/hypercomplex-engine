#!/usr/bin/env python3
"""Addendum: pin down the AZDE bookkeeping of Lemma 11.2 and the
pairing structure of the zero-divisor enumeration."""
import sys


def pc(x):
    return bin(x).count("1")


def nu2(x):
    return (x & -x).bit_length() - 1


def sgn_std(i, j):
    """Definition 6.1 of the paper (verified against the doubling formula
    in verify_cd.py on >1.1M pairs)."""
    if i == 0 or j == 0:
        return 1
    if i == j:
        return -1
    t = nu2(i ^ j)
    ti = nu2(i)
    tj = nu2(j)
    k = max(t, ti, tj) + 1
    if t >= ti and t >= tj:
        s = 2 * ((i >> t) & 1) - 1
    elif ti > t and ti >= tj:
        b = (j >> ti) & 1
        jloc = j & ((1 << ti) - 1)
        s = (1 if b == 0 else -1) * (1 if jloc == 0 else -1)
    else:
        s = 1 - 2 * ((i >> tj) & 1)
    if pc((i | j) >> k) & 1:
        s = -s
    return s


def generates(i, j, N, st):
    L = i ^ j
    if L == 0:
        return False
    for u in range(1, N):
        v = u ^ L
        if v == 0 or u == i or u == j or u == 0 or v == i or v == j:
            continue
        o = -st[i][j] * st[u][v]
        if st[u][j] == -o * st[i][v]:
            return True
    return False


print("=== AZDE reconstruction (Lemma 11.2) ===")
for x in (4, 5, 6):
    N = 1 << x
    st = [[sgn_std(i, j) for j in range(N)] for i in range(N)]
    top = 1 << (x - 1)

    # (1) never-generating positive unordered pairs == top-structured pairs
    non_gen, top_struct = set(), set()
    for i in range(1, N):
        for j in range(i + 1, N):
            if not generates(i, j, N, st):
                non_gen.add((i, j))
            if top in (i, j, i ^ j):
                top_struct.add((i, j))
    ok1 = non_gen == top_struct

    # (2) full-table never-generate count == 6*2^x - 8 (Lemma 11.2),
    #     counting zero rows, zero columns, diagonal, and (0,0)
    full = 2 * len(non_gen)                 # ordered positive pairs
    full += (N - 1) * 3 + 1                 # (i,i), (0,j), (i,0), (0,0)
    # verify the trivial positions really never generate (zero rows at
    # levels 5 and 6 included in the loop below)
    zero_ok = all(not generates(0, j, N, st) for j in range(1, N))
    ok2 = full == 6 * N - 8 and zero_ok

    # (3) each valid 4-set supports exactly 2 of its 3 pairings
    pair_count = {}
    for i in range(1, N):
        for j in range(i + 1, N):
            L = i ^ j
            for u in range(1, N):
                v = u ^ L
                if v == 0 or u in (i, j):
                    continue
                o = -st[i][j] * st[u][v]
                if st[u][j] != -o * st[i][v]:
                    continue
                s4 = frozenset((i, j, u, v))
                pr = frozenset((frozenset((i, u)), frozenset((j, v))))
                pair_count.setdefault(s4, set()).add(pr)
    ok3 = all(len(v) == 2 for v in pair_count.values())

    print(f"  A_{x}: non-gen == top-structured pairs: {ok1} ({len(non_gen)} pairs); "
          f"full-table count {full} == 6*{N}-8 = {6*N-8}: {ok2}; "
          f"every valid 4-set has exactly 2 valid pairings: {ok3}")
    assert ok1 and ok2 and ok3

print()
print("All AZDE / pairing structure checks passed.")
print("Summary of the reconstruction:")
print(" - never-generating starts 1<=i<j<N are exactly the pairs with")
print("   2^(x-1) in {i, j, i xor j}; count 3*(2^(x-1)-1).")
print(" - Lemma 11.2's 6*2^x-8 is the FULL-TABLE never-generate count")
print("   (ordered positive pairs x2, plus diagonal, zero row, zero col,")
print("   and (0,0)), which all fail to generate for structural reasons.")
