# Benchmark run

## Environment

| Item | Value |
|---|---|
| hypercomplex-engine | 0.4.2.1 |
| Python | 3.13.5 |
| NumPy | 2.5.3 |
| Platform | Debian 13 |
| CPU | 12th Gen Intel(R) Core(TM) i5-12450H |
| Mode | full |
| File | [examples/benchmark/run_benchmarks.py](https://github.com/maher1719/hypercomplex-engine/blob/main/examples/benchmark/run_benchmarks.py) |

## 1. One-off product from a cold interpreter (n = 12)

Each value is the median of 9 fresh processes (range in brackets); the import itself is not timed.

| Method | Median | Range |
|---|---|---|
| Build table n=12 + one lookup | 87.66 ms | [81.69 ms - 124.59 ms] |
| First `multiply` call, fast engine | 27.15 us | [25.96 us - 28.86 us] |
| First `multiply` call, holographic engine | 28.03 us | [27.17 us - 28.86 us] |

## 2. Warm per-call cost and break-even against a table (n = 12)

| Inputs | Table lookup | Raw fast | Raw holo | Facade fast | Facade holo |
|---|---|---|---|---|---|
| single-bit indices (e1024 x e2048) | 0.14 us | 0.38 us | 0.34 us | 3.67 us | 3.65 us |
| dense random 12-bit odd indices | 0.15 us | 0.37 us | 1.92 us | 3.60 us | 5.33 us |

Table build for n = 12: **108.8 ms** (48 MiB).

| Compared with | Table (build + lookups) becomes cheaper after |
|---|---|
| raw fast engine | ~493,108 products |
| facade `multiply` (fast) | ~31,446 products |

## 3. Scaling with index bit-length (dense random odd indices, raw engines)

`n` = number of bits in each index = number of doublings; the algebra has 2^n basis elements. 'Local slope' is log(t2/t1)/log(n2/n1) against the previous row: 0 = constant, 1 = linear, 2 = quadratic. 'Equals fast?' compares the holographic result with the fast engine's result on the same indices.

### Fast engine

| n (bits) | Time per product | Local slope |
|---|---|---|
| 1,000 | 0.72 us | - |
| 10,000 | 2.81 us | 0.59 |
| 100,000 | 19.48 us | 0.84 |
| 1,000,000 | 189.25 us | 0.99 |
| 3,000,000 | 590.62 us | 1.04 |

### Holographic engine

| n (bits) | Time per product | Local slope | Equals fast? |
|---|---|---|---|
| 1,000 | 177.85 us | - | yes |
| 4,000 | 879.62 us | 1.15 | yes |
| 16,000 | 6.50 ms | 1.44 | yes |
| 64,000 | 65.62 ms | 1.67 | yes |
| 256,000 | 861.66 ms | 1.86 | yes |
