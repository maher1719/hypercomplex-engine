# hypercomplex-engine


Fast, validated multiplication and table generation for Cayley–Dickson algebras.

This library provides:

- Full multiplication table generation for standard, split, and dual algebras.
- O(n) holographic table-free multiplication.
- O(1) fast bitwise multiplication.
- Integer, graded, and LaTeX notation formatting.
- CSV export for tables.
- A simple facade API for everyday use.
- Direct low-level classes for advanced use.

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

Install in editable mode:

```bash
pip install -e .
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
# +o_{12}
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

For `split`, `dim` is optional and can often be inferred from the indices.

For `dual` and `dual_split`, `dim` is required.

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

Copyright (c) 2026 Maher Ben Abdessalem
