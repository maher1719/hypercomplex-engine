#!/usr/bin/env python3
"""
Compression ladders across the signature families of the Cayley-Dickson
framework, the block-d dissection of the split family, and bidirectional
(generative) ladder walks.  Companion to verify_cd.py / verify_azde.py.

Families (per-level coefficient of the d*b term in the doubling formula):
  ORD        sig = {}          (all -1)             -- paper Def 2.1
  SPLIT_TOP  sig = {n: +1}     (final doubling)     -- paper Thm 4.1 / Cor 6.9
  SPLIT_BOT  sig = {1: +1}     (split-complex root) -- classical split family
  FULL_SPLIT sig = all +1      (every doubling)
  DUAL       A_n[e]/(e^2)      (paper Section 9, block composition)
  TENSOR     Kronecker sign composition (paper Section 10)

Symbols of the compression are canonical nested tuples, so alphabets,
dictionaries and cross-n transfers are plain set comparisons.
"""
import sys

FAIL = []


def check(name, ok, detail=""):
    print(("PASS" if ok else "FAIL"), "-", name, detail)
    if not ok:
        FAIL.append(name)


def pc(x):
    return bin(x).count("1")


def nu2(x):
    return (x & -x).bit_length() - 1


# ------------------------------------------------------------------
# ground truth: the doubling formula with a per-level signature
# ------------------------------------------------------------------
def conj_full(x, L):
    if L == 0:
        return dict(x)
    h = 1 << (L - 1)
    a = {k: v for k, v in x.items() if k < h}
    b = {k - h: v for k, v in x.items() if k >= h}
    r = conj_full(a, L - 1)
    for k, v in b.items():
        r[k + h] = r.get(k + h, 0) - v
    return r


def mult_gen(x, y, L, sig):
    """(a,b)(c,d) = (ac + s*d*b, da + bc*) with s = sig.get(L, -1)."""
    if not x or not y:
        return {}
    if L == 0:
        return {0: x.get(0, 0) * y.get(0, 0)}
    h = 1 << (L - 1)
    ax = {k: v for k, v in x.items() if k < h}
    bx = {k - h: v for k, v in x.items() if k >= h}
    ay = {k: v for k, v in y.items() if k < h}
    by = {k - h: v for k, v in y.items() if k >= h}
    s = sig.get(L, -1)
    c1 = mult_gen(ax, ay, L - 1, sig)
    for k, v in mult_gen(conj_full(by, L - 1), bx, L - 1, sig).items():
        c1[k] = c1.get(k, 0) + s * v
    c2 = mult_gen(by, ax, L - 1, sig)
    for k, v in mult_gen(bx, conj_full(ay, L - 1), L - 1, sig).items():
        c2[k] = c2.get(k, 0) + v
    for k, v in c2.items():
        c1[k + h] = c1.get(k + h, 0) + v
    return c1


def sign_gen(i, j, n, sig):
    r = mult_gen({i: 1}, {j: 1}, n, sig)
    assert len(r) == 1 and next(iter(r)) == (i ^ j), (i, j, sig, r)
    return 1 if next(iter(r.values())) > 0 else -1


# verified evaluators (from verify_cd.py, already checked vs the algebra)
def sgn_std(i, j):
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


# ------------------------------------------------------------------
# compression with canonical nested-tuple symbols
# ------------------------------------------------------------------
def compress(M):
    h = len(M) // 2
    return [[(M[2 * a][2 * b], M[2 * a][2 * b + 1],
              M[2 * a + 1][2 * b], M[2 * a + 1][2 * b + 1])
             for b in range(h)] for a in range(h)]


def levels_of(M):
    out = [M]
    while len(M) > 1:
        M = compress(M)
        out.append(M)
    return out


def ladder(M):
    return [len({v for row in L for v in row}) for L in levels_of(M)]


def sym(i, j, k, S):
    """The level-k compression symbol of the 2^k x 2^k block at (i, j),
    built recursively in the same quadrant order as compress()."""
    if k == 0:
        return S(i, j)
    h = 1 << (k - 1)
    return (sym(i, j, k - 1, S), sym(i, j + h, k - 1, S),
            sym(i + h, j, k - 1, S), sym(i + h, j + h, k - 1, S))


def fam_table(n, sig):
    return [[sign_gen(i, j, n, sig) for j in range(1 << n)]
            for i in range(1 << n)]


_ord_cache, _sp_cache = {}, {}


def ord_t(n):
    if n not in _ord_cache:
        _ord_cache[n] = [[sgn_std(i, j) for j in range(1 << n)]
                         for i in range(1 << n)]
    return _ord_cache[n]


def sp_t(n):
    if n not in _sp_cache:
        _sp_cache[n] = [[sgn_split(i, j, n) for j in range(1 << n)]
                        for i in range(1 << n)]
    return _sp_cache[n]


# ==================================================================
print("=== 1. Regressions ===")
ok = all(sign_gen(i, j, 5, {}) == sgn_std(i, j)
         for i in range(32) for j in range(32))
check("mult_gen(sig={}) == sgn_std, all pairs of A5", ok)
ok = all(sign_gen(i, j, 5, {5: 1}) == sgn_split(i, j, 5)
         for i in range(32) for j in range(32))
check("mult_gen(sig={5:+1}) == sgn_split, all pairs of A5", ok)
L6 = levels_of(ord_t(6))
ok = all(L6[3][a][b] == sym(8 * a, 8 * b, 3, sgn_std)
         for a in range(8) for b in range(8))
check("sym() == compressed-table entries at level 3 (A6)", ok)
ok = all(L6[1][a][b] == sym(2 * a, 2 * b, 1, sgn_std)
         for a in range(32) for b in range(32))
check("sym() == compressed-table entries at level 1 (A6)", ok)

# ==================================================================
print("=== 2. Family ladders ===")
BOT = {1: 1}
FULL = {L: 1 for L in range(1, 9)}
fam_cache = {}


def bot_t(n):
    if n not in fam_cache:
        fam_cache[("bot", n)] = fam_table(n, BOT)
    return fam_cache[("bot", n)]


def full_t(n):
    if ("full", n) not in fam_cache:
        fam_cache[("full", n)] = fam_table(n, FULL)
    return fam_cache[("full", n)]


for n in (4, 5, 6, 7):
    print(f"   A_{n}:  ordinary={ladder(ord_t(n))}  split-top={ladder(sp_t(n))}"
          f"  split-bottom={ladder(bot_t(n))}  full-split={ladder(full_t(n))}")
print(f"   A_8:  ordinary={ladder(ord_t(8))}  split-top={ladder(sp_t(8))}")
print(f"   full-split small n: n=1:{ladder(full_t(1))} n=2:{ladder(full_t(2))}"
      f" n=3:{ladder(full_t(3))}")

# ==================================================================
print("=== 3. ORDINARY: bidirectional walk (the generative ladder) ===")
cls = lambda a, b: sym(8 * a, 8 * b, 3, sgn_std)
sigma, ok = {}, True
for a in range(64):
    for b in range(64):
        c = cls(a, b)
        quad = (cls(2 * a, 2 * b), cls(2 * a, 2 * b + 1),
                cls(2 * a + 1, 2 * b), cls(2 * a + 1, 2 * b + 1))
        if sigma.setdefault(c, quad) != quad:
            ok = False
check("substitution sigma well-defined on the 12 classes "
      "(4096 parents, children < 128)", ok, f"|classes|={len(sigma)}")

m4, ok4 = {}, True
for a in range(24):
    for b in range(24):
        s4 = sym(16 * a, 16 * b, 4, sgn_std)
        if m4.setdefault(cls(a, b), s4) != s4:
            ok4 = False
check("class determines the level-4 symbol (positions < 24)", ok4)

c2s1, okc = {}, True
for a in range(128):
    for b in range(128):
        s1 = sym(2 * a, 2 * b, 1, sgn_std)
        if c2s1.setdefault(cls(a, b), s1) != s1:
            okc = False
check("class determines the level-1 symbol (12 -> 6, positions < 128)", okc,
      f"6 symbols: {len(set(c2s1.values()))}")

# the walk: regenerate A_8 (256x256) from the apex + rulebook
cur = {(0, 0): cls(0, 0)}
rungs = [1]
for _step in range(7):          # levels 8 -> 1
    nxt = {}
    for (a, b), c in cur.items():
        q = sigma[c]
        nxt[(2 * a, 2 * b)] = q[0]
        nxt[(2 * a, 2 * b + 1)] = q[1]
        nxt[(2 * a + 1, 2 * b)] = q[2]
        nxt[(2 * a + 1, 2 * b + 1)] = q[3]
    cur = nxt
    rungs.append(len(set(cur.values())))
T = [[0] * 256 for _ in range(256)]
for (a, b), c in cur.items():
    q = c2s1[c]
    T[2 * a][2 * b] = q[0]
    T[2 * a][2 * b + 1] = q[1]
    T[2 * a + 1][2 * b] = q[2]
    T[2 * a + 1][2 * b + 1] = q[3]
ok = all(T[i][j] == sgn_std(i, j) for i in range(256) for j in range(256))
check("A_8 sign table (65,536 cells) regenerated exactly from "
      "apex + sigma + class->sym1", ok)
print("   walk rungs (apex -> base):", rungs, "then",
      len(set(c2s1.values())), "symbols, 2 signs")

# ==================================================================
print("=== 4. SPLIT (top doubling): ladder and block-d dissection ===")
T7s, T7o = sp_t(7), ord_t(7)
L7s, L7o = levels_of(T7s), levels_of(T7o)
h = 64
ok = True
for i in range(64):
    for j in range(64):
        d = T7s[h + i][h + j]
        pred = 1 if i == j else (-1 if i == 0 else (1 if j == 0 else sgn_std(i, j)))
        if d != pred:
            ok = False
check("block d of split A_7 = std A_6 with diag +1, row0 -1, col0 +1, "
      "interior std (all 4,096 cells)", ok)

# where does the split table differ from the ordinary one, per level?
for k in (1, 2, 3, 4, 5):
    As, Ao = L7s[k], L7o[k]
    m = len(As)
    B = m // 2
    out_ok = all(As[a][b] == Ao[a][b]
                 for a in range(m) for b in range(m) if a < B or b < B)
    alpha_o = {v for row in Ao for v in row}
    new = {v for row in As for v in row} - alpha_o
    posns = [(a, b) for a in range(m) for b in range(m)
             if As[a][b] in new and a >= B and b >= B]
    onspine = sum(1 for a, b in posns
                  if a - B == b - B or a == B or b == B)
    print(f"   level {k}: outside-d identical to ordinary: {out_ok}; "
          f"new symbols: {len(new)}; occurrences in d: {len(posns)} "
          f"(of which on d-spine: {onspine})")

# the new level-1 tiles, spelled out
o1 = {v for row in L7o[1] for v in row}
n1 = sorted({v for row in L7s[1] for v in row} - o1)
print("   new level-1 tiles of the split ladder:", n1)

# product invariant of Statement 8.6(1)
viol = [(a, b) for a in range(64) for b in range(64)
        if T7s[2 * a][2 * b] * T7s[2 * a][2 * b + 1] *
        T7s[2 * a + 1][2 * b] * T7s[2 * a + 1][2 * b + 1] != -1]
vd = sum(1 for a, b in viol if a >= 32 and b >= 32)
print(f"   2x2 blocks violating product = -1: {len(viol)} "
      f"(inside d: {vd}); e.g. {viol[:6]}")

# region grammar: is the split symbol determined by (ordinary symbol, region)?
for k in (2, 3, 4, 5, 6):
    As, Ao = L7s[k], L7o[k]
    m = len(As)
    B = m // 2
    gram, okg = {}, True
    for a in range(m):
        for b in range(m):
            if a < B or b < B:
                reg = "out"
            else:
                aa, bb = a - B, b - B
                reg = ("diag" if aa == bb else
                       "row0" if aa == 0 else ("col0" if bb == 0 else "dint"))
            if gram.setdefault((Ao[a][b], reg), As[a][b]) != As[a][b]:
                okg = False
    check(f"level {k}: split symbol determined by (ordinary symbol, region)",
          okg, f"rules={len(gram)}")

# cross-n uniformity: shape + aligned sub-apex alphabets
for nn in range(4, 9):
    lz = ladder(sp_t(nn))
    check(f"split A_{nn} ladder has the uniform shape "
          "[2, 6, 12^(n-4), 10, 4, 1]",
          lz == [2, 6] + [12] * (nn - 4) + [10, 4, 1], f"{lz}")
T8s = sp_t(8)
L8s = levels_of(T8s)
for k in range(1, 5):
    s8 = {v for row in L8s[k] for v in row}
    s7 = {v for row in L7s[k] for v in row}
    check(f"split A_8 level-{k} alphabet == split A_7's "
          "(aligned plateau levels)", s8 == s7, f"{len(s8)} symbols")

# ==================================================================
print("=== 5. Native split families: coherence ===")
for name, tab in (("split-bottom", bot_t), ("full-split", full_t)):
    sig = BOT if name == "split-bottom" else FULL
    ok = all(sign_gen(i, j, nn, sig) == sign_gen(i, j, nn + 1, sig)
             for nn in (4, 5) for i in range(1 << nn) for j in range(1 << nn))
    check(f"{name}: signs independent of ambient dimension (n=4,5)", ok)

for name, T7x, sig in (("split-bottom", bot_t(7), BOT),
                       ("full-split", full_t(7), FULL)):
    S = lambda i, j, T=T7x: T[i][j]
    clsx = lambda a, b: sym(8 * a, 8 * b, 3, S)
    sig_map, ok = {}, True
    for a in range(8):
        for b in range(8):
            c = clsx(a, b)
            quad = (clsx(2 * a, 2 * b), clsx(2 * a, 2 * b + 1),
                    clsx(2 * a + 1, 2 * b), clsx(2 * a + 1, 2 * b + 1))
            if sig_map.setdefault(c, quad) != quad:
                ok = False
    check(f"{name}: substitution well-defined on level-3 classes "
          "(parents < 8)", ok, f"classes seen={len(sig_map)}")
    # walk at n = 5 (32x32), rulebook from sigma + class->sym1
    c2s1, okc = {}, True
    for a in range(16):
        for b in range(16):
            s1 = sym(2 * a, 2 * b, 1, S)
            if c2s1.setdefault(clsx(a, b), s1) != s1:
                okc = False
    cur = {(0, 0): clsx(0, 0)}
    rungs = [1]
    for _step in range(4):
        nxt = {}
        for (a, b), c in cur.items():
            q = sig_map[c]
            nxt[(2 * a, 2 * b)] = q[0]
            nxt[(2 * a, 2 * b + 1)] = q[1]
            nxt[(2 * a + 1, 2 * b)] = q[2]
            nxt[(2 * a + 1, 2 * b + 1)] = q[3]
        cur = nxt
        rungs.append(len(set(cur.values())))
    T = [[0] * 32 for _ in range(32)]
    for (a, b), c in cur.items():
        q = c2s1[c]
        T[2 * a][2 * b] = q[0]
        T[2 * a][2 * b + 1] = q[1]
        T[2 * a + 1][2 * b] = q[2]
        T[2 * a + 1][2 * b + 1] = q[3]
    Tt = fam_table(5, sig)
    ok = all(T[i][j] == Tt[i][j] for i in range(32) for j in range(32))
    check(f"{name}: 32x32 table regenerated from apex + rulebook", ok,
          f"walk rungs {rungs}")

# ==================================================================
print("=== 6. Dual family ===")


def dual_table(n):
    N = 1 << n
    A = ord_t(n)
    T = [[0] * (2 * N) for _ in range(2 * N)]
    for i in range(N):
        for j in range(N):
            T[i][j] = A[i][j]
            T[i][N + j] = A[i][j]
            T[N + i][j] = A[i][j]
    return T


for n in (3, 4, 5, 6):
    dl = ladder(dual_table(n))
    al = ladder(ord_t(n))
    exp = [3] + [x + 1 for x in al[1:]] + [1]
    check(f"dual over A_{n}: ladder = [3] + (ordinary + 1 per rung) + [1]",
          dl == exp, f"{dl}")

# ==================================================================
print("=== 7. Tensor products (Kronecker sign composition) ===")


def kron(B, A):
    P, N = len(B), len(A)
    T = [[0] * (P * N) for _ in range(P * N)]
    for p in range(P):
        for q in range(P):
            for i in range(N):
                for j in range(N):
                    T[p * N + i][q * N + j] = B[p][q] * A[i][j]
    return T


C2 = [[1, 1], [1, -1]]
D2 = [[1, 1], [1, 0]]
Cp = [[1, 1], [1, 1]]
H4 = ord_t(2)
O8 = ord_t(3)
check("dual over A_4 == kron(D, A_4) (Section 9 = Section 10 instance)",
      dual_table(4) == kron(D2, ord_t(4)))
for name, B, A in (("C x C", C2, C2), ("C x H", C2, H4), ("H x C", H4, C2),
                   ("C x O", C2, O8), ("O x C", O8, C2), ("D x C", D2, C2),
                   ("C' x C", Cp, C2), ("C' x O", Cp, O8),
                   ("H x H", H4, H4), ("C x A5", C2, ord_t(5))):
    print(f"   {name}: ladder = {ladder(kron(B, A))}")

print()
print("FAILURES:", FAIL if FAIL else "none")
sys.exit(1 if FAIL else 0)
