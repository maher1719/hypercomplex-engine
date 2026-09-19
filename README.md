# hypercomplex-engine

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://img.shields.io/badge/tests-62%20passed-success)](#)
[![GitHub Repo](https://img.shields.io/badge/GitHub-maher1719%2Fhypercomplex--engine-black?logo=github)](https://github.com/maher1719/hypercomplex-engine)
[![PyPI version](https://badge.fury.io/py/hypercomplex-engine.svg)](https://pypi.org/project/hypercomplex-engine/)
[![Downloads](https://static.pepy.tech/badge/hypercomplex-engine)](https://pepy.tech/project/hypercomplex-engine)


Fast, validated multiplication and table generation for Cayley–Dickson algebras.

This library provides the computational substrate for high-dimensional hypercomplex algebra, featuring:

- **Full multiplication table generation** for standard, split, and dual algebras.
- **O(n) holographic** table-free recursive descent multiplication.
- **O(1) fast bitwise** closed-form multiplication.
- **Integer, graded, and LaTeX** notation formatting.
- **CSV export** for tables (matrix and long formats).
- A **simple facade API** for everyday use, and direct low-level classes for advanced physics/math engines.

---

## 📄 Publications & Preprints

This library serves as the formal verification substrate and computational engine for the following mathematical preprints:

**1. The Sign Structure of Cayley–Dickson and Split Algebras By Blocks**
*Proves the OPMT (Ordered-Pair Multiplication Table) sign laws, block decomposition, and the O(1) closed-form sign evaluator implemented in the `fast` engine of this library.*
* 📊 **Figshare:** [10.6084/m9.figshare.33705022](https://doi.org/10.6084/m9.figshare.33705022)
* 📦 **Zenodo:** [10.5281/zenodo.22051873](https://doi.org/10.5281/zenodo.22051873)

---

## Supported algebras

| Kind | Description |
|---|---|
| `standard` | Ordinary Cayley–Dickson algebras: real, complex, quaternions, octonions, sedenions, ... |
| `split` | Split Cayley–Dickson algebras: standard parent plus one split doubling at the top |
| `dual` | Dual extension of a standard algebra, with ε² = 0 |
| `dual_split` | Dual extension of a split algebra, with ε² = 0 |

---

## Installation

Clone the repository:

```bash
git clone https://github.com/maher1719/hypercomplex-engine.git
cd hypercomplex-engine
```
Install with pip:

```bash
pip install hypercomplex-engine
```

Install in editable mode:

```bash
pip install -e
```

Run the tests:

```bash
pytest -v
```

---

# Basic Use

The simplest way to use the library is through the top-level facade API.

```python
from hypercomplex import (
    build_table,
    multiply,
    format_element,
    print_table,
    export_csv,
)
```

---

## Build a table

```python
table = build_table("standard", 3)
```

This builds the octonion multiplication table.

Dimensions:

```text
n = 0 -> real numbers,       dimension 1
n = 1 -> complex numbers,    dimension 2
n = 2 -> quaternions,        dimension 4
n = 3 -> octonions,          dimension 8
n = 4 -> sedenions,          dimension 16
```

---

## Print a table

```python
print_table(table, title="Octonions", mode="integer")
```

Example output style:

```text
      e0  e1  e2  e3  e4  e5  e6  e7
e0 | +e0 +e1 +e2 +e3 +e4 +e5 +e6 +e7
e1 | +e1 -e0 +e3 -e2 +e5 -e4 -e7 +e6
...
```

You can also use graded notation:

```python
print_table(table, title="Octonions", mode="graded")
```

Example:

```text
       1  o1  o2  o3  o4  o5  o6  o7
1   | +1 +o1 +o2 +o3 +o4 +o5 +o6 +o7
o1  | +o1 -1 ...
...
```

---

## Export a table to CSV

Matrix-style CSV:

```python
export_csv(
    "octonions_graded.csv",
    table,
    mode="graded",
    csv_mode="matrix",
)
```

Long-format CSV for data analysis:

```python
export_csv(
    "octonions_long.csv",
    table,
    mode="integer",
    csv_mode="long",
)
```

The long format produces rows like:

```csv
i,j,sign,index
0,0,1,0
0,1,1,1
1,0,1,1
1,1,-1,0
...
```

---

## Multiply two basis elements

```python
result = multiply("standard", (1, 1), (1, 2))

print(result)
# (1, 3)
```

This means:

```text
e1 * e2 = +e3
```

Format the result:

```python
print(format_element(result, mode="integer"))
# +e3

print(format_element(result, mode="graded"))
# +o12

print(format_element(result, mode="latex"))
# +e_{12}
```

Note:

```text
integer mode uses the basis index:
    e3

graded mode uses the generator decomposition:
    index 3 = binary 011 = generators 1 and 2 = o12
```

---

## Split multiplication

```python
result = multiply("split", (1, 1), (1, 1), dim=1)

print(result)
# (1, 0)
```

In split-complex numbers:

```text
e1² = +e0
```

---

## Dual multiplication

```python
# eps*e0 represented as local tuple: (sign, local_index, eps_flag)
eps_e0 = (1, 0, 1)

result = multiply("dual", (1, 0), eps_e0, dim=1)

print(result)
# (1, 0, 1)

print(format_element(result, mode="integer"))
# +eps
```

Nilpotency:

```python
result = multiply("dual", eps_e0, eps_e0, dim=1)

print(result)
# (0, 0, 0)
```

This means:

```text
ε² = 0
```

---

# Intermediate Use

The facade API is enough for most users.

For more control, you can choose the computation engine and work directly with tables or multipliers.

---

## Engines

The `multiply` function supports two engines:

```python
multiply(kind, a, b, dim=None, engine="fast")
```

| Engine | Complexity | Description |
|---|---:|---|
| `"fast"` | O(1) | Bitwise closed-form sign evaluator |
| `"holographic"` | O(n) | Recursive block descent |

Example:

```python
from hypercomplex import multiply

a = (1, 3)
b = (1, 5)

fast_result = multiply("standard", a, b, engine="fast")
holo_result = multiply("standard", a, b, engine="holographic")

assert fast_result == holo_result
```

---

## Algebra kinds

```python
multiply("standard", a, b)
multiply("split", a, b, dim=3)
multiply("dual", a, b, dim=3)
multiply("dual_split", a, b, dim=3)
```

For `standard`, `dim` is not needed.

For `split`, `dual` and `dual_split`, `dim` is required.

---

## Table builders directly

```python
from hypercomplex import (
    StandardTableBuilder,
    SplitTableBuilder,
    DualTableBuilder,
)

standard_builder = StandardTableBuilder()
split_builder = SplitTableBuilder()
dual_builder = DualTableBuilder()

signs, indices = standard_builder.build(3)
signs, indices = split_builder.build(3)
signs, indices, eps = dual_builder.build(2, split=False)
```

Return conventions:

```text
standard:
    signs, indices

split:
    signs, indices

dual:
    signs, indices, eps
```

For dual tables:

- `signs[i, j]` is the sign.
- `indices[i, j]` is the local base index.
- `eps[i, j]` is the epsilon flag.

---

## Multipliers directly

```python
from hypercomplex import (
    StandardHolographic,
    SplitHolographic,
    DualHolographic,
)

holo = StandardHolographic()
result = holo.multiply((1, 1), (1, 2))

print(result)
# (1, 3)
```

Split:

```python
split_holo = SplitHolographic()
result = split_holo.multiply((1, 2), (1, 2), dim=2)

print(result)
# (1, 0)
```

Dual:

```python
dual_holo = DualHolographic(split=False)

result = dual_holo.multiply((1, 0), (1, 2), dim=1)

print(result)
# (1, 0, 1)
```

---

## Fast O(1) multipliers directly

```python
from hypercomplex import (
    FastStandard,
    FastSplit,
    FastDual,
)

fast = FastStandard()

result = fast.multiply((1, 1), (1, 2))

print(result)
# (1, 3)
```

Split:

```python
fast_split = FastSplit()

result = fast_split.multiply((1, 2), (1, 2), dim=2)

print(result)
# (1, 0)
```

Dual:

```python
fast_dual = FastDual(split=False)

result = fast_dual.multiply((1, 0), (1, 0, 1), dim=1)

print(result)
# (1, 0, 1)
```

---

## Formatting modes

| Mode | Example |
|---|---|
| `"integer"` | `+e5` |
| `"graded"` | `+o13` |
| `"latex"` | `+e_{13}` |
| `"latex_integer"` | `+e_{5}` |
| `"latex_graded"` | `+o_{13}` |

Example:

```python
from hypercomplex import format_element

element = (-1, 5)

print(format_element(element, mode="integer"))
# -e5

print(format_element(element, mode="graded"))
# -o13

print(format_element(element, mode="latex"))
# -e_{5}

print(format_element(element, mode="latex_integer"))
# -e_{5}

print(format_element(element, mode="latex_graded"))
# -o_{13}
```

---

# Advanced Use

This section is for contributors, benchmarking, physics engines, and symbolic pipelines.

---

## Direct low-level imports

If you prefer explicit imports:

```python
from hypercomplex.core.table_builder import (
    StandardTableBuilder,
    SplitTableBuilder,
    DualTableBuilder,
)

from hypercomplex.core.holographic import (
    StandardHolographic,
    SplitHolographic,
    DualHolographic,
)

from hypercomplex.core.fast import (
    FastStandard,
    FastSplit,
    FastDual,
)

from hypercomplex.printer import (
    CDFormat,
    CDTablePrinter,
)
```

---

## Cross-validating O(1) against the full table

```python
from hypercomplex import StandardTableBuilder, FastStandard

builder = StandardTableBuilder()
fast = FastStandard()

n = 4
signs, indices = builder.build(n)

dim = 1 << n

for i in range(dim):
    for j in range(dim):
        fast_sign, fast_idx = fast.multiply_indices(i, j)

        assert int(signs[i, j]) == fast_sign
        assert int(indices[i, j]) == fast_idx
```

This proves that the O(1) evaluator agrees with the O(4^n) table builder.

---

## Cross-validating split O(1) against the split table

```python
from hypercomplex import SplitTableBuilder, FastSplit

builder = SplitTableBuilder()
fast = FastSplit()

n = 4
signs, indices = builder.build(n)

dim = 1 << n

for i in range(dim):
    for j in range(dim):
        fast_sign, fast_idx = fast.multiply_indices(i, j, dim=n)

        assert int(signs[i, j]) == fast_sign
        assert int(indices[i, j]) == fast_idx
```

---

## Dual local and global indices

For dual multiplication, the total dimension is:

```text
2^(dim + 1)
```

The epsilon bit is bit `dim`.

Example for `dim=1`:

```text
lower half: 0, 1        base elements
upper half: 2, 3        epsilon elements
```

The dual multipliers accept both:

```python
# global index tuple
(1, 2)

# local tuple with epsilon flag
(1, 0, 1)
```

Both represent ε·e₀ when `dim=1`.

The output convention is:

```python
(sign, local_index, eps_flag)
```

This makes formatting easy:

```python
from hypercomplex import format_element

result = (1, 0, 1)

print(format_element(result, mode="integer"))
# +eps

print(format_element(result, mode="latex"))
# +\epsilon
```

---

## Using the fast engine in a physics loop

For simulations, avoid building large tables. Use the fast engine directly.

```python
from hypercomplex import FastStandard

fast = FastStandard()

def basis_product(i: int, j: int):
    sign, index = fast.multiply((1, i), (1, j))
    return sign, index

sign, index = basis_product(1, 2)

print(sign, index)
# 1 3
```

For octonionic or higher-dimensional simulations, this avoids O(4^n) memory.

---

## Table size warning

Full table generation grows as:

```text
entries = 4^n
```

where `n` is the dimension exponent.

| n | Dimension | Entries |
|---:|---:|---:|
| 0 | 1 | 1 |
| 1 | 2 | 4 |
| 2 | 4 | 16 |
| 3 | 8 | 64 |
| 4 | 16 | 256 |
| 5 | 32 | 1,024 |
| 6 | 64 | 4,096 |
| 8 | 256 | 65,536 |
| 10 | 1,024 | 1,048,576 |
| 12 | 4,096 | 16,777,216 |

For large dimensions, prefer:

```python
engine="fast"
```

or:

```python
engine="holographic"
```

---

# API Reference

## Top-level functions

### `build_table(kind, n)`

Builds a multiplication table.

```python
table = build_table("standard", 3)
```

Returns:

```text
standard:
    (signs, indices)

split:
    (signs, indices)

dual:
    (signs, indices, eps)

dual_split:
    (signs, indices, eps)
```

---

### `multiply(kind, a, b, dim=None, engine="fast")`

Multiplies two basis elements.

```python
result = multiply("standard", (1, 1), (1, 2))
```

Returns:

```text
standard:
    (sign, index)

split:
    (sign, index)

dual:
    (sign, local_index, eps_flag)

dual_split:
    (sign, local_index, eps_flag)
```

---

### `format_element(element, mode="integer")`

Formats a basis element tuple.

```python
format_element((1, 3), mode="integer")
# "+e3"

format_element((1, 3), mode="graded")
# "+o12"
```

---

### `print_table(table, title=None, limit=None, mode="integer")`

Prints a table.

```python
table = build_table("standard", 2)
print_table(table, mode="graded")
```

---

### `export_csv(path, table, mode="integer", csv_mode="matrix")`

Exports a table to CSV.

```python
table = build_table("standard", 3)

export_csv(
    "octonions.csv",
    table,
    mode="graded",
    csv_mode="matrix",
)
```

CSV modes:

| `csv_mode` | Output |
|---|---|
| `"matrix"` | Spreadsheet-style grid |
| `"long"` | One row per product |

---

## Algebra kinds

| Kind | Meaning |
|---|---|
| `"standard"` | Ordinary Cayley–Dickson |
| `"split"` | Split Cayley–Dickson |
| `"dual"` | Dual extension of standard algebra |
| `"dual_split"` | Dual extension of split algebra |

Aliases:

```text
standard: "std", "ordinary", "o"
split:    "s"
dual:     "d", "dual_standard"
dual_split: "split_dual", "ds"
```

---

## Engines

| Engine | Aliases | Complexity |
|---|---|---:|
| `"fast"` | `"o1"`, `"bitwise"`, `"constant"` | O(1) Word-RAM |
| `"holographic"` | `"on"`, `"descent"` | O(n) |

---

# Mathematical Background

## Basis product rule

For standard and split Cayley–Dickson algebras:

```text
e_i * e_j = sign(i, j) * e_{i XOR j}
```

The index is always:

```text
i XOR j
```

The sign is determined by the OPMT block laws.

---

## Standard doubling formula

```text
(a, b)(c, d) = (ac - d* b, da + b c*)
```

with conjugation:

```text
e0* = e0
ek* = -ek for k > 0
```

---

## Split doubling formula

```text
(a, b)(c, d) = (ac + d* b, da + b c*)
```

The only difference from the standard construction is the sign of the `d* b` term.

This causes Block d signs to invert relative to the standard algebra.

---

## Block decomposition

Each multiplication table splits into four blocks:

```text
[ a  b ]
[ c  d ]
```

where:

```text
Block a: e_i * e_j
Block b: e_i * (e_j ℓ)
Block c: (e_i ℓ) * e_j
Block d: (e_i ℓ) * (e_j ℓ)
```

For standard algebras:

```text
Block d interior sign = -σ_a
```

For split algebras:

```text
Block d interior sign = +σ_a
```

---

## Dual numbers

Dual algebras adjoin ε such that:

```text
ε² = 0
```

Multiplication rules:

```text
e_i * e_j       = parent product
e_i * (ε e_j)   = ε (e_i e_j)
(ε e_i) * e_j   = ε (e_i e_j)
(ε e_i) * (ε e_j) = 0
```

---

# Complexity

| Operation | Complexity | Memory |
|---|---:|---:|
| Full table generation | O(4^n) | O(4^n) |
| Holographic multiplication | O(n) | O(1) |
| Fast bitwise multiplication | O(1) Word-RAM | O(1) |

For arbitrary-precision integers, the fast evaluator uses O(n) bit operations, where:

```text
n = ceil(log2(max(i, j) + 1))
```

---

# Testing

Run all tests:

```bash
pytest -v
```

Run specific test files:

```bash
pytest tests/test_mega_mother.py -v
pytest tests/test_fast_mode.py -v
```

The test suite validates:

- Basis notation conversion.
- Input validation.
- Standard table generation.
- Split table generation.
- Dual table generation.
- Holographic O(n) multiplication.
- Fast O(1) multiplication.
- Cross-validation between tables and multipliers.
- Facade API behavior.
- CSV export.

---

# Repository Structure

```text
hypercomplex-engine/
├── examples/
│   ├── direct_implementation/
│   │   └── full_table_builder_simple.py
│   └── uses/
│       ├── outputs/
│       └── use.ipynb
├── hypercomplex/
│   ├── core/
│   │   ├── basis_element.py
│   │   ├── basis_notation.py
│   │   ├── validation.py
│   │   ├── table_builder/
│   │   │   ├── common.py
│   │   │   ├── standard.py
│   │   │   ├── split.py
│   │   │   └── dual.py
│   │   ├── holographic/
│   │   │   ├── standard.py
│   │   │   ├── split.py
│   │   │   └── dual.py
│   │   └── fast/
│   │       ├── bit_utils.py
│   │       ├── fast_standard.py
│   │       ├── fast_split.py
│   │       └── fast_dual.py
│   ├── printer/
│   │   ├── cd_format.py
│   │   └── cd_table_printer.py
│   ├── facade.py
│   └── __init__.py
├── tests/
│   ├── test_algebra.py
│   ├── test_fast_mode.py
│   ├── test_holographic_vs_table.py
│   └── test_mega_mother.py
├── LICENSE
├── README.md
└── pyproject.toml
```

---

# Citation

If you use this engine in your research, physics simulations, or geometric deep learning models, please cite the underlying theoretical preprints:

```bibtex
@article{ben abdessalem2026,
author = "maher ben abdessalem",
title = "{A Proven Sign Law for Cayley-Dickson Algebras: Ordinary, Split, dual Constructions and their computational proofs and implementations}",
year = "2026",
month = "9",
url = "https://figshare.com/articles/preprint/A_Proven_Sign_Law_for_Cayley-Dickson_Algebras_Ordinary_and_Split_Constructions/33705022",
doi = "10.6084/m9.figshare.33705022.v5"
}

```

---

# License

MIT License.

See [`LICENSE`](LICENSE) for details.

Copyright (c) 2026 Maher Ben Abdessalem
