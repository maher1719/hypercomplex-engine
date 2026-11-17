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


/-! ## Block rules for the closed form (R1 is trivial). h = 2^n, a b < h. -/

def r2Holds (n : Nat) : Bool :=
  (List.range (2^n)).all fun a => (List.range (2^n)).all fun b =>
    closedSign a (b + 2^n) == closedSign b a
def r3Holds (n : Nat) : Bool :=
  (List.range (2^n)).all fun a => (List.range (2^n)).all fun b =>
    closedSign (a + 2^n) b == closedSign a b * conjSign b
def r4Holds (n : Nat) : Bool :=
  (List.range (2^n)).all fun a => (List.range (2^n)).all fun b =>
    closedSign (a + 2^n) (b + 2^n) == -(conjSign b * closedSign b a)

#eval (List.range 10).map fun n => (n, r2Holds n, r3Holds n, r4Holds n)

theorem R2 (n a b : Nat) (ha : a < 2^n) (hb : b < 2^n) :
    closedSign a (b + 2^n) = closedSign b a := sorry
theorem R3 (n a b : Nat) (ha : a < 2^n) (hb : b < 2^n) :
    closedSign (a + 2^n) b = closedSign a b * conjSign b := sorry
theorem R4 (n a b : Nat) (ha : a < 2^n) (hb : b < 2^n) :
    closedSign (a + 2^n) (b + 2^n) = -(conjSign b * closedSign b a) := sorry

/-- Main theorem: the closed form equals recursive doubling in every dimension. -/
theorem closed_eq_rec : ∀ (n i j : Nat), i < 2^n → j < 2^n → closedSign i j = sgnRec n i j := by
  intro n
  induction n with
  | zero =>
    intro i j hi hj
    have hi0 : i = 0 := by simp at hi; omega
    have hj0 : j = 0 := by simp at hj; omega
    subst hi0; subst hj0
    decide
  | succ n ih =>
    intro i j hi hj
    have h2 : 2 ^ (n+1) = 2 * 2 ^ n := by rw [Nat.pow_succ]; omega
    rw [h2] at hi hj
    simp only [sgnRec]
    by_cases hi' : i < 2 ^ n
    · by_cases hj' : j < 2 ^ n
      · simp only [hi', hj', if_true]
        exact ih i j hi' hj'
      · simp only [hi', hj', if_true, if_false]
        obtain ⟨b, rfl⟩ : ∃ b, j = b + 2 ^ n := ⟨j - 2 ^ n, by omega⟩
        have hb : b < 2 ^ n := by omega
        rw [R2 n i b hi' hb, Nat.add_sub_cancel]
        exact ih b i hb hi'
    · by_cases hj' : j < 2 ^ n
      · simp only [hi', hj', if_true, if_false]
        obtain ⟨a, rfl⟩ : ∃ a, i = a + 2 ^ n := ⟨i - 2 ^ n, by omega⟩
        have ha : a < 2 ^ n := by omega
        rw [R3 n a j ha hj', Nat.add_sub_cancel, ih a j ha hj']
      · simp only [hi', hj', if_false]
        obtain ⟨a, rfl⟩ : ∃ a, i = a + 2 ^ n := ⟨i - 2 ^ n, by omega⟩
        obtain ⟨b, rfl⟩ : ∃ b, j = b + 2 ^ n := ⟨j - 2 ^ n, by omega⟩
        have ha : a < 2 ^ n := by omega
        have hb : b < 2 ^ n := by omega
        rw [R4 n a b ha hb, Nat.add_sub_cancel, Nat.add_sub_cancel, ih b a hb ha]

end CDSign
#print axioms CDSign.closed_eq_rec
