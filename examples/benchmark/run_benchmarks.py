#!/usr/bin/env python
"""Reproducible benchmarks for hypercomplex-engine.

Usage:
    python run_benchmarks.py            # full run (about 1-2 minutes)
    python run_benchmarks.py --quick    # smaller sizes, fewer repeats

Every number is the MINIMUM over several repeats (the standard way to reduce
scheduler noise), measured after a warm-up. Results are printed as Markdown
tables so they can be pasted into BENCHMARKS.md.

What is measured
    raw engine   FastStandard / StandardHolographic .multiply_indices(i, j)
                 (calculation only, no input validation)
    facade       multiply("standard", (1, i), (1, j), engine=...)
                 (validation + normalisation + calculation: what users call)
    table        build_table("standard", n), then numpy indexing
"""
import argparse
import math
import platform
import random
import statistics
import subprocess
import sys
import time
from importlib.metadata import version

import numpy as np

from hypercomplex import FastStandard, StandardHolographic, build_table, multiply

sys.set_int_max_str_digits(0)

FAST, HOLO = FastStandard(), StandardHolographic()


# --------------------------------------------------------------------------- helpers

def best_us(fn, number, repeat):
    """Minimum over `repeat` runs of `number` calls, in microseconds per call."""
    fn()  # warm-up
    best = float("inf")
    for _ in range(repeat):
        t = time.perf_counter()
        for _ in range(number):
            fn()
        best = min(best, time.perf_counter() - t)
    return best / number * 1e6


def dense_index(bits, rng):
    """Random `bits`-bit odd integer with the top bit set (no structural shortcuts)."""
    return rng.getrandbits(bits) | 1 | (1 << (bits - 1))


def fmt_us(us):
    if us >= 1e6:
        return f"{us / 1e6:.2f} s"
    if us >= 1e3:
        return f"{us / 1e3:.2f} ms"
    return f"{us:.2f} us"


def cpu_name():
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.lower().startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def table(headers, rows):
    print("| " + " | ".join(headers) + " |")
    print("|" + "|".join("---" for _ in headers) + "|")
    for r in rows:
        print("| " + " | ".join(str(c) for c in r) + " |")
    print()


# --------------------------------------------------------------------------- sections

def environment(args):
    print("## Environment\n")
    table(["Item", "Value"], [
        ("hypercomplex-engine", version("hypercomplex-engine")),
        ("Python", platform.python_version()),
        ("NumPy", np.__version__),
        ("Platform", platform.platform()),
        ("CPU", cpu_name()),
        ("Mode", "quick" if args.quick else "full"),
    ])


COLD_TABLE = (
    "import time\nfrom hypercomplex import build_table\n"
    "t=time.perf_counter(); S,I=build_table('standard',{n}); a=int(S[5][10]); b=int(I[5][10]); "
    "print(time.perf_counter()-t)\n"
)
COLD_CALL = (
    "import time\nfrom hypercomplex import multiply\n"
    "t=time.perf_counter(); multiply('standard',(1,5),(1,10),engine='{engine}'); "
    "print(time.perf_counter()-t)\n"
)


def _fresh(code):
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    return float(out.stdout.strip())


def cold_start(n, runs):
    """One product from a fresh interpreter: each method in its own process (no order effects)."""
    print(f"## 1. One-off product from a cold interpreter (n = {n})\n")
    print(f"Each value is the median of {runs} fresh processes (range in brackets); "
          "the import itself is not timed.\n")
    rows = []
    for label, code in (
        (f"Build table n={n} + one lookup", COLD_TABLE.format(n=n)),
        ("First `multiply` call, fast engine", COLD_CALL.format(engine="fast")),
        ("First `multiply` call, holographic engine", COLD_CALL.format(engine="holographic")),
    ):
        xs = sorted(_fresh(code) for _ in range(runs))
        med = statistics.median(xs)
        rows.append((label, fmt_us(med * 1e6), f"[{fmt_us(xs[0] * 1e6)} - {fmt_us(xs[-1] * 1e6)}]"))
    table(["Method", "Median", "Range"], rows)


def per_call(n, quick):
    print(f"## 2. Warm per-call cost and break-even against a table (n = {n})\n")
    rng = random.Random(12345)
    number, repeat = (5000, 3) if quick else (20000, 7)

    t0 = time.perf_counter()
    signs, idx = build_table("standard", n)
    build_s = time.perf_counter() - t0

    cases = {
        "single-bit indices (e1024 x e2048)": (1024, 2048),
        f"dense random {n}-bit odd indices": (dense_index(n, rng), dense_index(n, rng)),
    }
    rows = []
    for name, (i, j) in cases.items():
        assert FAST.multiply_indices(i, j) == HOLO.multiply_indices(i, j)
        look = best_us(lambda: (int(signs[i, j]), int(idx[i, j])), number * 5, repeat)
        rf = best_us(lambda: FAST.multiply_indices(i, j), number * 5, repeat)
        rh = best_us(lambda: HOLO.multiply_indices(i, j), number, repeat)
        ff = best_us(lambda: multiply("standard", (1, i), (1, j), engine="fast"), number, repeat)
        fh = best_us(lambda: multiply("standard", (1, i), (1, j), engine="holographic"), number, repeat)
        rows.append((name, fmt_us(look), fmt_us(rf), fmt_us(rh), fmt_us(ff), fmt_us(fh)))
        last = (look, rf, ff)
    table(["Inputs", "Table lookup", "Raw fast", "Raw holo", "Facade fast", "Facade holo"], rows)

    look, rf, ff = last
    print(f"Table build for n = {n}: **{build_s * 1e3:.1f} ms** "
          f"({(signs.nbytes + idx.nbytes) / 2**20:.0f} MiB).\n")
    rows = []
    for label, per in (("raw fast engine", rf), ("facade `multiply` (fast)", ff)):
        d = per - look
        rows.append((label, f"~{build_s * 1e6 / d:,.0f} products" if d > 0 else "never (engine is cheaper)"))
    table(["Compared with", "Table (build + lookups) becomes cheaper after"], rows)


def scaling(fast_sizes, holo_sizes, quick):
    print("## 3. Scaling with index bit-length (dense random odd indices, raw engines)\n")
    print("`n` = number of bits in each index = number of doublings; the algebra has 2^n basis "
          "elements. 'Local slope' is log(t2/t1)/log(n2/n1) against the previous row: "
          "0 = constant, 1 = linear, 2 = quadratic. 'Equals fast?' compares the holographic "
          "result with the fast engine's result on the same indices.\n")
    rng = random.Random(2026)
    for name, engine, sizes in (("Fast engine", FAST, fast_sizes), ("Holographic engine", HOLO, holo_sizes)):
        rows, prev = [], None
        for bits in sizes:
            i, j = dense_index(bits, rng), dense_index(bits, rng)
            once = best_us(lambda: engine.multiply_indices(i, j), 1, 1)
            number = max(1, min(2000, int(2e5 / max(once, 1))))
            repeat = 2 if quick or once > 1e5 else 5
            us = best_us(lambda: engine.multiply_indices(i, j), number, repeat)
            slope = "-" if prev is None else f"{math.log(us / prev[1]) / math.log(bits / prev[0]):.2f}"
            row = [f"{bits:,}", fmt_us(us), slope]
            if engine is HOLO:
                row.append("yes" if FAST.multiply_indices(i, j) == HOLO.multiply_indices(i, j) else "NO")
            rows.append(tuple(row))
            prev = (bits, us)
        headers = ["n (bits)", "Time per product", "Local slope"] + (["Equals fast?"] if engine is HOLO else [])
        print(f"### {name}\n")
        table(headers, rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="smaller sizes and fewer repeats")
    ap.add_argument("--cold-runs", type=int, default=None)
    ap.add_argument("--n", type=int, default=12, help="dimension exponent for the table comparison")
    a = ap.parse_args()

    if a.quick:
        fast_sizes = [10**3, 10**4, 10**5, 10**6]
        holo_sizes = [10**3, 4 * 10**3, 16 * 10**3]
    else:
        fast_sizes = [10**3, 10**4, 10**5, 10**6, 3 * 10**6]
        holo_sizes = [10**3, 4 * 10**3, 16 * 10**3, 64 * 10**3, 256 * 10**3]

    print("# Benchmark run\n")
    environment(a)
    cold_start(a.n, a.cold_runs or (3 if a.quick else 9))
    per_call(a.n, a.quick)
    scaling(fast_sizes, holo_sizes, a.quick)


if __name__ == "__main__":
    main()
