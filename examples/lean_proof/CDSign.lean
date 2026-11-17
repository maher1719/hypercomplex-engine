/-!
# Cayley-Dickson signs: closed form vs. recursive doubling (core Lean 4, no Mathlib)

Convention: (a,b)(c,d) = (ac - d*b, da + b c*),  conj negates every coordinate except the first.
For basis elements e_i, e_j we have e_i * e_j = s * e_(i XOR j), and we compare the sign s.
-/

namespace CDSign

/-- Sign of conjugation on basis element e_x: e_0* = e_0, e_x* = -e_x otherwise. -/
def conjSign (x : Nat) : Int := if x = 0 then 1 else -1

/-- Recursive doubling. `sgnRec n i j` is the sign of e_i * e_j in dimension 2^n. -/
def sgnRec : Nat → Nat → Nat → Int
  | 0, _, _ => 1
  | n+1, i, j =>
    let h := 2 ^ n
    if i < h then
      if j < h then sgnRec n i j                      -- (a,0)(c,0) = (ac, 0)
      else sgnRec n (j - h) i                         -- (a,0)(0,d) = (0, d a)
    else
      if j < h then sgnRec n (i - h) j * conjSign j   -- (0,b)(c,0) = (0, b c*)
      else -(conjSign (j - h) * sgnRec n (j - h) (i - h)) -- (0,b)(0,d) = (-d* b, 0)

/-- Lowest set bit index (fuel-based so the kernel can evaluate it). -/
def nu2Aux : Nat → Nat → Nat
  | 0, _ => 0
  | f+1, x => if x % 2 = 1 then 0 else if x = 0 then 0 else 1 + nu2Aux f (x / 2)

def nu2 (x : Nat) : Nat := nu2Aux x x

/-- Number of set bits (fuel-based). -/
def popAux : Nat → Nat → Nat
  | 0, _ => 0
  | f+1, x => if x = 0 then 0 else x % 2 + popAux f (x / 2)

def popcount (x : Nat) : Nat := popAux x x

/-- Bit t of x. -/
def bit (x t : Nat) : Nat := (x / 2 ^ t) % 2

/-- The closed-form rule (a direct transcription of `multiply_indices`). -/
def closedSign (i j : Nat) : Int :=
  if i = 0 ∨ j = 0 then 1
  else if i = j then -1
  else
    let t  := nu2 (i ^^^ j)
    let ti := nu2 i
    let tj := nu2 j
    let r : Nat × Int :=
      if t ≥ ti ∧ t ≥ tj then
        (t + 1, 2 * (bit i t : Int) - 1)
      else if ti ≥ tj then
        let b : Nat := bit j ti
        let jloc := j % 2 ^ ti
        let δ : Int := if jloc = 0 then 1 else -1
        (ti + 1, if b = 0 then δ else -δ)
      else
        (tj + 1, 1 - 2 * (bit i tj : Int))
    if popcount ((i ||| j) / 2 ^ r.1) % 2 = 1 then -r.2 else r.2

/-- Exhaustive check in dimension 2^n. -/
def checkAll (n : Nat) : Bool :=
  (List.range (2 ^ n)).all fun i =>
    (List.range (2 ^ n)).all fun j => closedSign i j == sgnRec n i j

/-- Sanity: the recursion reproduces the quaternion/octonion facts. -/
example : sgnRec 2 1 2 = 1 := by decide
example : sgnRec 2 2 1 = -1 := by decide
example : sgnRec 1 1 1 = -1 := by decide

-- Kernel-checked (no `native_decide`): every pair of basis elements, dimensions 2^1 .. 2^6.
theorem verified_dim_2 : checkAll 1 = true := by decide +kernel
theorem verified_dim_4 : checkAll 2 = true := by decide +kernel
theorem verified_dim_8 : checkAll 3 = true := by decide +kernel
theorem verified_dim_16 : checkAll 4 = true := by decide +kernel
theorem verified_dim_32 : checkAll 5 = true := by decide +kernel
theorem verified_dim_64 : checkAll 6 = true := by decide +kernel

-- Larger runs by evaluation (compiled code, not kernel proofs).
#eval (List.range 10).map fun n => (n + 1, checkAll (n + 1))


end CDSign

#print axioms CDSign.verified_dim_64
