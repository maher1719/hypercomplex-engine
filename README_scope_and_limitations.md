## Scope, strengths and limitations

> **Status: 0.4.x, pre-1.0.** The API can still change. Read this section before building on the package.

### What this package is

It computes the product of two **signed basis elements** of a Cayley–Dickson algebra: `±e_i × ±e_j → ±e_k` (or zero in dual algebras). It is a low-level component, a sign and index calculator, not a number type. There are no coefficients, no sums, and no vectors.

### Advantages

- **One job, small surface.** About 2,000 lines, one dependency (NumPy), and multipliers that keep no per-call state.
- **Two interchangeable engines.**
  - `fast` is a closed-form bitwise evaluator. It does constant work per call for word-sized indices, and cost grows only with bit length for arbitrary-precision indices.
  - `holographic` is an O(n) recursive descent, useful as a cross-check.
- **Works where tables cannot.** Only the two indices are needed, so indices hundreds of bits long are fine. Full tables need 4^n entries.
- **Cross-validated.** The test suite checks the table builders, the `fast` engine and the `holographic` engine against each other.
- **Strict about types.** `bool`, `float`, `str` and `list` inputs are rejected. Results are plain Python `int`s.
- **Tested on Python 3.10–3.14** and with NumPy 1.22 and newer.

### Input rules (0.4.x behavior)

| Kind | Element form | `dim` | Index range |
|---|---|---|---|
| `standard` | `(sign, index)` | not used | any integer ≥ 0 (no upper bound, since there is no `dim` to check against) |
| `split` | `(sign, index)` | **required** | `0 ≤ index < 2**dim` |
| `dual`, `dual_split` | `(sign, global_index)` or `(sign, local_index, eps)` | **required** | global: `0 ≤ i < 2**(dim+1)`; local: `0 ≤ i < 2**dim` |

- Elements must be **tuples**. Integers only, and NumPy integers are accepted.
- `sign` should be `-1` or `+1`. Note: `multiply` currently also accepts `0` as a formal zero and returns zero. **Do not rely on this**; it may be removed.
- Dual elements have two spellings of the same thing. At `dim=3`, `(1, 9)` and `(1, 1, 1)` both mean `+ε·e1` (global index `8 + 1`).
- **Known looseness:** at `dim=3`, `(1, 9, 1)` and even `(1, 9, 0)` are also accepted and read as `+ε·e1`. The second one contradicts itself, so don't rely on it.
- The dual result is always `(sign, local_index, eps)`, and `ε·ε` gives `(0, 0, 0)`.

### What the package cannot detect (your responsibility)

Elements are plain tuples, so they do not remember which algebra produced them.

- Feeding a result from one `dim` into an algebra of a different `dim` is not detected if the index happens to be in range.
- Mixing `standard` and `split` elements is not detected. The same tuple means different things: `e2·e2` is `-1` in the standard algebra, `+1` in `split` with `dim=2`, and `-1` in `split` with `dim>=3`.
- For `standard`, nothing tells you an index is "too big for the octonions".

### Not supported

- **No chain multiplication and no expressions.** `multiply` is strictly binary. You can pass a result into the next call yourself, but you choose the bracketing. From `n=3` (octonions) upward multiplication is **not associative**, so `(ab)c` and `a(bc)` can differ in sign.
- **No coefficients, sums, vectors, arrays, or parsing** of expressions like `e1+e2+e12`. This is planned for a separate package built on top of this one.
- **No addition, conjugation, norm, inverse, or division.**
- **`split` means one split doubling on top of a standard parent.** Other sign patterns are not offered.
- **`dual` and `dual_split`** extend a standard or split parent with a central `ε`, where `ε² = 0`.

### Algebra facts to keep in mind

| n | Standard | Split |
|---|---|---|
| 1 | complex, commutative | split-complex, commutative, has zero divisors |
| 2 | quaternions, associative, **not commutative** | split-quaternions, associative, not commutative |
| 3 | octonions, alternative, **not associative** | split-octonions, alternative, not associative, has zero divisors |
| ≥ 4 | sedenions and beyond: **zero divisors**, no longer alternative | no longer alternative |

### Convention

Products follow the standard doubling `(a, b)(c, d) = (ac − d*b, da + bc*)` (split algebras flip the sign of the `d*b` term at the top level). Other conventions give isomorphic algebras with **different tables**, so comparing against another source may show sign differences after relabeling. For example, here `e1·e2 = e3` and `e1·e6 = −e7`.

### Direct low-level classes

`FastStandard`, `FastSplit`, `FastDual`, `StandardHolographic`, `SplitHolographic` and `DualHolographic` are exported, but they validate **differently** from each other and from `multiply`. For example, `FastSplit` accepts a zero element while the other three reject it. Their `multiply_indices` methods do only minimal checks. **Prefer `multiply`.** These classes are not a stable API and are expected to become internal.

### Tables and memory

Full tables have `4**n` entries. `build_table` refuses requests above a default memory budget (`max_bytes`, 256 MiB):

- `standard` and `split` are allowed up to `n = 13`.
- `dual` and `dual_split` are allowed up to `n = 12`, since they are one doubling larger.

Use `estimate_table_bytes(kind, n)` to check a size, and pass a larger `max_bytes` (or `None`) if you really mean it. In dual tables, `sign == 0` marks `ε·ε`.

### Stability and planned changes

- Pin your dependency, for example `hypercomplex-engine>=0.4.1,<0.5`.
- *Planned for 0.5.0, subject to change:* stricter inputs (no zero elements, one dual form), and low-level classes made internal. Breaking changes will be listed in the release notes.
- Do not treat the package as 1.0-stable until it has been used by futur a palnned package the vector-and-expression package.
