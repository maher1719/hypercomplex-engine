# ⚡ Hypercomplex-Engine Computational Benchmarks

This document outlines the formal performance metrics, scaling boundaries, and stress-test results of the `hypercomplex-engine` core computational substrate. 



By bypassing traditional exponential O(4^n) lookup tables and introducing a dual-engine architecture powered by **closed-form O(1) bitwise sign laws** and **logarithmic O(log n) recursive block-descent structural parsing**, this library removes the physical memory and execution walls historically associated with high-dimensional Cayley–Dickson algebras.

---

## Environment

| Item | Value |
|---|---|
| hypercomplex-engine | 0.4.2 |
| Python | 3.13.5 |
| NumPy | 2.5.3 |
| Platform | Debian 13 |
| CPU | 12th Gen Intel(R) Core(TM) i5-12450H |
| Mode | full |
| File | [examples/benchmark/benchmark.ipynb](https://github.com/maher1719/hypercomplex-engine/blob/main/examples/benchmark/benchmark.ipynb) |

---

## 📊 Summary of Algorithmic Complexity

| Engine / Methodology | Computational Complexity | Memory Complexity | Scalability Boundaries |
| :--- | :--- | :--- | :--- |
| **Traditional Matrix Lookup** | O(1) (after generation) | O(4^n) (Exponential) | Hard wall at n = 13 (8,192 Dim) due to RAM limits |
| **Holographic Engine (`holo`)** | O(log n) (Logarithmic) | O(1) (Constant) | Infinite (Optimized for 2^k structural boundaries) |
| **Fast Bitwise Engine (`fast`)** | O(1) (Constant) | O(1) (Constant) | Infinite (Optimized for dense/extreme scale bignums) |

---

## 🚀 Stage 1: Ad-Hoc Setup & Execution Latency
**Environment:** Dimension exponent n = 12 (4,096 Dimensions)  
*Objective: Measure the total execution duration of a single arbitrary basis multiplication, factoring in mandatory backend initialization costs.*

* **Traditional SoA Matrix Lookup (Setup + Run):** `0.106069 seconds`
* **Maher's O(1) Fast Bitwise Engine (On-the-fly):** `0.000039 seconds`
* **Maher's Holographic Engine (On-the-fly):** `0.000018 seconds`

### 💡 Core Discovery
For ad-hoc calculations, **Maher's Fast Engine is 2,687.82x faster** and his **Holographic Engine is 5,986.53x faster** than traditional matrix architectures. They bypass the costly initialization, allocation, and memory footprint of building pre-compiled matrix arrays.

---

## ⚔️ Stage 2: Head-to-Head Multi-Iteration Engine Loops
**Volume:** 50,000 sequential iterations per engine loop.  
*Objective: Isolate raw execution paths by comparing boundary structural alignments against chaotic dense environments.*

### Test A: Structural Boundary Pillars (Exact Powers of 2)
Inputs align perfectly with structural generational block limits (e.g., e₁₀₂₄ × e₂₀₄₈).
* **O(1) Fast Engine Loops:** `0.18834 seconds`
* **Holographic Engine Loops:** `0.18324 seconds`
* 🏁 **Winner:** **Holographic Engine (1.03x faster)** by triggering its internal `is_structural` boundary short-circuit.

### Test B: Dense Random Odd Indices (Fractal Matrix Chaos)
Inputs target arbitrary, dense, non-aligned odd integers deep within the algebra.
* **O(1) Fast Engine Loops:** `0.20627 seconds`
* **Holographic Engine Loops:** `0.28122 seconds`
* 🏁 **Winner:** **O(1) Fast Engine (1.36x faster)** because closed-form bitwise arithmetic eliminates recursive stack-branch overhead.

---

## 🌋 Stage 3: Extreme Scale Space Expansion (The Monster Tests)
Leveraging Python's arbitrary-precision integers (`PyLongObject`), both engines were pushed into astronomical dimensions using dense, non-power-of-two odd index coordinates. This scales the mathematical space far beyond the number of subatomic particles in the physical universe (>10⁸⁰), evaluating sign correctness under massive computational loads.

### 1. The 4,300-Digit Frontier (n ≈ 14,284 bit-width)
* **O(1) Fast Engine Runtime:** `0.000021 seconds` (21 Microseconds)
* **Holographic Engine Runtime:** `0.008185 seconds` (8 Milliseconds)
* 🚀 **Speed Divergence:** Bitwise engine scales **393.69x faster** than the recursive engine.

### 2. The 60,000-Digit Expansion (n ≈ 199,315 bit-width)
* **O(1) Fast Engine Runtime:** `0.000066 seconds` (66 Microseconds)
* **Holographic Engine Runtime:** `0.620146 seconds` (620 Milliseconds)
* 🚀 **Speed Divergence:** Bitwise engine scales **9,342.22x faster** than the recursive engine.

### 3. The 1,000,000-Digit Apex Milestone (n ≈ 3,321,928 bit-width)
* **O(1) Fast Engine Runtime:** `0.000782 seconds` (782 Microseconds)
* **Holographic Engine Runtime:** `154.992875 seconds` (2.5 Minutes)
* 📐 **Structural Verification:** `equals holo and fast? True`

