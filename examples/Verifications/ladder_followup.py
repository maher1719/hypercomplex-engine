#!/usr/bin/env python3
"""Follow-up v2: same dissection, with timing checkpoints and flushed output."""
import time


def pc(x):
    return bin(x).count("1")


def nu2(x):
    return (x & -x).bit_length() - 1


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


def sym(i, j, k, S):
    if k == 0:
        return S(i, j)
    h = 1 << (k - 1)
    return (sym(i, j, k - 1, S), sym(i, j + h, k - 1, S),
            sym(i + h, j, k - 1, S), sym(i + h, j + h, k - 1, S))


def tick(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


T0 = time.time()
tick("start")

n = 7
T7s = [[sgn_split(i, j, n) for j in range(1 << n)] for i in range(1 << n)]
tick("T7s built")
T7o = [[sgn_std(i, j) for j in range(1 << n)] for i in range(1 << n)]
tick("T7o built")
L7s, L7o = levels_of(T7s), levels_of(T7o)
tick("levels of A_7 computed")
L4s = levels_of([[sgn_split(i, j, 4) for j in range(16)] for i in range(16)])
L4o = levels_of([[sgn_std(i, j) for j in range(16)] for i in range(16)])
tick("levels of A_4 computed")

print("=== the 10-rung under the microscope ===", flush=True)
for tag, Lv_s, Lv_o in (("split A_4 level 2", L4s[2], L4o[2]),
                        ("split A_7 level 5", L7s[5], L7o[5])):
    alpha_o = {v for row in Lv_o for v in row}
    alpha_s = {v for row in Lv_s for v in row}
    new = alpha_s - alpha_o
    lost = alpha_o - alpha_s
    pos_new = [(a, b) for a in range(len(Lv_s)) for b in range(len(Lv_s))
               if Lv_s[a][b] in new]
    pos_lost = [(a, b) for a in range(len(Lv_o)) for b in range(len(Lv_o))
                if Lv_o[a][b] in lost]
    print(f"  {tag} ({len(Lv_s)}x{len(Lv_s)}): new={len(new)} at {pos_new}; "
          f"lost={len(lost)} (ordinary positions: {pos_lost})", flush=True)
tick("10-rung dissection done")

print("=== block-d 2x2 patterns vs ordinary tiles ===", flush=True)
blk = lambda T, i, j: (T[i][j], T[i][j + 1], T[i + 1][j], T[i + 1][j + 1])
print("  split d-diag 2x2 at (72,72):", blk(T7s, 72, 72), flush=True)
print("  split d-diag 2x2 at (80,80):", blk(T7s, 80, 80), flush=True)
print("  ordinary at same positions :", blk(T7o, 72, 72), blk(T7o, 80, 80),
      flush=True)
o1 = {v for row in L7o[1] for v in row}
for t in sorted(o1):
    cnt = sum(row.count(t) for row in L7o[1])
    print("   tile", t, f"x{cnt}", flush=True)
for pat in (blk(T7s, 72, 72), blk(T7s, 80, 80)):
    where = [(a, b) for a in range(64) for b in range(64)
             if L7o[1][a][b] == pat]
    print(f"  {pat}: {'IS' if pat in o1 else 'NOT'} ordinary; "
          f"{len(where)} occurrences, e.g. {where[:4]}", flush=True)
tick("tile census done")

print("=== the 4-rung (level n-1) content: A_7 ===", flush=True)
g_s, g_o = L7s[6], L7o[6]
same = sum(1 for a in range(2) for b in range(2) if g_s[a][b] == g_o[a][b])
print(f"  2x2 quadrant grid: {same}/4 entries identical to ordinary; "
      f"distinct symbols split={len(set(v for r in g_s for v in r))}, "
      f"ord={len(set(v for r in g_o for v in r))}", flush=True)
tick("4-rung done")

print("=== cross-n anchoring of d-region 8x8 classes ===", flush=True)
S7 = lambda i, j: sgn_split(i, j, 7)
S8 = lambda i, j: sgn_split(i, j, 8)
alpha_o3 = {v for row in L7o[3] for v in row}
for name, a7, b7, a8, b8 in (
        ("d-corner (h,h)", 64, 64, 128, 128),
        ("d-edge (h,3h/2)", 64, 96, 128, 192),
        ("d-interior (3h/2,3h/2)", 96, 96, 192, 192),
        ("d-interior (5h/4,5h/4)", 80, 80, 160, 160)):
    c7 = sym(a7, b7, 3, S7)
    c8 = sym(a8, b8, 3, S8)
    print(f"  {name}: class(n=7)==class(n=8): {c7 == c8}; "
          f"ordinary tile: {c7 in alpha_o3}", flush=True)
tick("anchoring done")

h = 128
ok = all(sgn_split(i, j, 8) == sgn_std(i, j)
         for i in range(0, 256, 7) for j in range(0, 256, 7)
         if i < h or j < h)
print("  split A_8 outside block d == ordinary (sparse sample):", ok,
      flush=True)
tick("done")
