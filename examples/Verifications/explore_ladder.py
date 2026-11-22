#!/usr/bin/env python3
"""Where else does the compression ladder live?  Companion computations.

1.  CD alphabet ladder for A_n, n = 5..11: extends the paper's verified
    pattern [2, 6, 12^(n-4), 9, 4, 1] one level beyond the earlier check.
2.  The six canonical 2x2 sign-block types (printed for flavor).
3.  The SAME ladder run on the SPLIT family (Corollary 6.9 evaluator):
    does the plateau move when the construction is split?
4.  Distinct 8x8 canonical blocks in the split family (vs 12 for ordinary).
5.  Thue-Morse factor complexity p(L), L = 1..16: the 1-D engine GROWS;
    the 2-D CD table SATURATES.  Quantifies the contrast.
6.  Zero-divisor bookkeeping: which statistic equals (2^x-2)(2^x-4)(2^x-8)/16.
7.  Lemma 11.2 as projective geometry over F_2: never-generating pairs
    should be exactly 3 pairs in each 2-dim subspace through the top
    basis vector, giving 3 * (2^(x-1) - 1) = 3 * #planes.
"""
import sys
from collections import defaultdict


def pc(x):
    return bin(x).count("1")


def nu2(x):
    return (x & -x).bit_length() - 1


def sgn_std(i, j):
    """Definition 6.1 (verified against the doubling formula on >1.1M pairs)."""
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


def sgn_split(i, j, m):
    """Corollary 6.9: split at the final doubling of A_m (verified, 87,296 pairs)."""
    if i == 0 or j == 0:
        return 1
    h = 1 << (m - 1)
    if i == j:
        return 1 if i >= h else -1
    if i >= h and j >= h:
        iloc, jloc = i - h, j - h
        if iloc == jloc:
            return 1
        if iloc == 0:
            return -1
        if jloc == 0:
            return 1
        return sgn_std(iloc, jloc)
    return sgn_std(i, j)


def ladder(f, n):
    """Alphabet sizes under iterated dyadic 2x2 block-type compression."""
    N = 1 << n
    M = [[f(i, j) for j in range(N)] for i in range(N)]
    sizes = []
    while True:
        vals = set()
        for row in M:
            vals.update(row)
        sizes.append(len(vals))
        if len(M) == 1:
            break
        h = len(M) // 2
        rank = {}
        C = [[0] * h for _ in range(h)]
        for a in range(h):
            Ma, Mb = M[2 * a], M[2 * a + 1]
            for b in range(h):
                tup = (Ma[2 * b], Ma[2 * b + 1], Mb[2 * b], Mb[2 * b + 1])
                r = rank.get(tup)
                if r is None:
                    r = len(rank) + 1
                    rank[tup] = r
                C[a][b] = r
        M = C
    return sizes


print("=== 1. ordinary-family ladder, n = 5..11 ===", flush=True)
for n in range(5, 12):
    L = ladder(sgn_std, n)
    tag = "OK" if L == [2, 6] + [12] * (n - 4) + [9, 4, 1] else "DEVIATES"
    print(f"  A_{n}: {L}   [{tag}]", flush=True)

print("=== 2. the six canonical 2x2 block types (ordinary) ===", flush=True)
types = set()
for a in range(64):
    for b in range(64):
        i, j = 2 * a, 2 * b
        types.add((sgn_std(i, j), sgn_std(i, j + 1),
                   sgn_std(i + 1, j), sgn_std(i + 1, j + 1)))
for t in sorted(types):
    print(f"    {t}", flush=True)
print(f"  -> {len(types)} types (16 quadruples have product -1; only these "
      f"occur: top row of each block restricted by Statement 8.6)", flush=True)

print("=== 3. SPLIT-family ladder, m = 5..8 ===", flush=True)
for m in range(5, 9):
    L = ladder(lambda i, j: sgn_split(i, j, m), m)
    print(f"  split A_{m}: {L}", flush=True)

print("=== 4. split family: distinct 8x8 blocks (ordinary has 12) ===", flush=True)
for m in (8, 10):
    N = 1 << m
    blocks = set()
    for a in range(N // 8):
        for b in range(N // 8):
            blocks.add(tuple(sgn_split(8 * a + r, 8 * b + c, m)
                             for r in range(8) for c in range(8)))
    print(f"  split A_{m}: distinct 8x8 blocks = {len(blocks)}", flush=True)

print("=== 5. Thue-Morse factor complexity p(L) (1-D engine, for contrast) ===",
      flush=True)
N = 1 << 18
w = bytearray([0])
while len(w) < N:
    w += bytes(1 - x for x in w)
w = bytes(w[:N])
vals = []
for L in range(1, 17):
    s = set(w[p:p + L] for p in range(N - L + 1))
    vals.append(len(s))
print(f"  p(1..16) = {vals}", flush=True)

print("=== 6. zero-divisor statistics vs formula (2^x-2)(2^x-4)(2^x-8)/16 ===",
      flush=True)


def zd_survey(x):
    N = 1 << x
    st = [[sgn_std(i, j) for j in range(N)] for i in range(N)]
    raw, sets, pairings, starts = set(), set(), set(), set()
    for i in range(1, N):
        for j in range(i + 1, N):
            L = i ^ j
            for u in range(1, N):
                v = u ^ L
                if v == 0 or u == i or u == j:
                    continue
                o = -st[i][j] * st[u][v]
                if st[u][j] != -o * st[i][v]:
                    continue
                raw.add((i, j, u, v, o))
                sets.add(frozenset((i, j, u, v)))
                pairings.add(frozenset((frozenset((i, u)), frozenset((j, v)))))
                starts.add((i, j))
    return len(raw), len(sets), len(pairings), len(starts)


for x in (4, 5, 6):
    formula = (2 ** x - 2) * (2 ** x - 4) * (2 ** x - 8) // 16
    raw, sets, pairings, starts = zd_survey(x)
    which = [nm for nm, v in (("raw", raw), ("4-sets", sets),
                              ("pairings", pairings), ("starts", starts))
             if v == formula]
    print(f"  x={x}: formula={formula}  raw={raw}  4-sets={sets}  "
          f"pairings={pairings}  generating_starts={starts}  "
          f"-> formula counts: {which}", flush=True)

print("=== 7. Lemma 11.2 as projective geometry over F_2 ===", flush=True)
for x in (4, 5, 6):
    N = 1 << x
    top = 1 << (x - 1)
    st = [[sgn_std(i, j) for j in range(N)] for i in range(N)]
    non_gen = set()
    for i in range(1, N):
        for j in range(i + 1, N):
            L = i ^ j
            gen = False
            for u in range(1, N):
                v = u ^ L
                if v == 0 or u == i or u == j:
                    continue
                o = -st[i][j] * st[u][v]
                if st[u][j] == -o * st[i][v]:
                    gen = True
                    break
            if not gen:
                non_gen.add((i, j))
    by_span = defaultdict(set)
    for i, j in non_gen:
        by_span[frozenset((0, i, j, i ^ j))].add((i, j))
    planes = {frozenset((0, top, w, top ^ w))
              for w in range(1, N) if w != top}
    ok = (set(by_span) == planes
          and all(len(v) == 3 for v in by_span.values())
          and len(non_gen) == 3 * (2 ** (x - 1) - 1)
          and len(planes) == 2 ** (x - 1) - 1)
    print(f"  x={x}: {len(non_gen)} never-generating pairs; "
          f"{len(planes)} F_2-planes through e_top; "
          f"exactly 3 pairs per plane, spans = planes: {ok}", flush=True)

print("done.", flush=True)
