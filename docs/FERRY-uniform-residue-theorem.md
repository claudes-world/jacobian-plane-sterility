# The Uniform Residue Theorem for split-edge cells (ferry document)

**From:** Claude (Fable), do-box jc lane, 2026-07-22 — **rev-2, 2026-07-23 ~03:15 ET** (audit repairs)
**To:** Liam, the originating desktop session, and the ChatGPT audit channel
**Status:** every displayed identity is sympy-verified exact (scripts named inline); the composed arguments are written here in full for audit. This document supersedes the Tail Exhaustion program (handoff Section 5) and closes the split-edge N ≥ 3 arc of the Line Theorem, conditional only on the audit of the three lemmas below.
**Rev-2.1 addendum (GPT-Pro audit, ~02:59 ET):** (4) the m < 0 terminal case of Theorem B, omitted in rev-1, is now written (one paragraph; the proof covers every m ∈ ℤ, making the claimed T4-independence real); (5) sharper kernel-exclusion at the germ grade via r | (N−1)k + fertility ⇒ r | k contradiction [GPT-Pro's arithmetic]; (6) Trace Lemma's unproven "τ polynomial" parenthetical deleted (x-freeness suffices); (7) two independently constructed exact fertile-cell pairs added as witnesses — (4,1,2) ratio −1/3 and (3,2,6) with germ in EXACT agreement with our verify_326_coset.py.
**Rev-2 changelog (from the 2026-07-23 cross-model audit round — GPT + chat-Claude, independently):** (1) the §5 bucket constants were displayed in the r = N specialization (true at the worked cells (3,1,3) and (3,2,3), false in general) — corrected to the general shift r − (N − 1); all conclusions survive because only bucket *differences* and the shift-invariant germ threshold v* = −(N−1)k enter. (2) Lemma 3's chain argument is refined by the ζ-coset decomposition, closing a genuine gap at cells with d = gcd(N−1, r) > 1 (e.g. the fertile (3,2,6)); new machine witness verify_326_coset.py 8/8. (3) Lemma 2's coordinate conventions (c-polygon vs u-polygon) clarified; the final admissibility inequalities are unchanged (they were machine-swept as plane-hull-equivalent all along).

## 0. Summary

Fix the split-edge chart of a cell (N, k, r): c = c_min + tail with c_min = Y + x^r Y^N (b₁ = β_N = 1 by torus normalization), pair equation J(w, c) = λ x^(−1−k) with λ ≠ 0, ζ := x^r Y^(N−1), weight v(x^a Y^b) = (N−1)a − rb. Chart solutions pull back to plane pairs through Y = x^σ y. A scaffold additionally carries P ∈ ℂ[W, C], D := deg_W P ≥ 1, with the master identity P_W(w, c) = μ x^(−m), μ = κ/λ ≠ 0, where m is the pair's exponent (m = σ − 1 − k in the realization at hand; no D-relation is assumed anywhere below).

Results, in logical order:

- **Theorem A (Shear Orbit).** At (3,1,3) the minimal obstruction ladder converges: its three forced shields are the binomial tail of the unipotent shear (x,Y) ↦ (x, Y+γx), the pair is exactly (w_germ, c_min) pulled through the shear, and it solves the pair equation to all orders. Lemma E-minimal is false: c_j = 0 for every j ≥ 4. Pair-level obstruction theory cannot close any fertile cell — the solution space is nonempty.
- **Lemma 1 (Fertility and germ uniqueness).** The cell has a germ iff r | Nk and r ∤ k; the germ is then the unique polynomial φ of degree Nk/r in ζ with w_germ = x^(−k) φ(ζ), and r ∤ k is precisely the condition that makes it unique on its line.
- **Lemma 2 (Hull rigidity).** For every polynomial tail and every plane realization σ, the fiber {c = t} (generic t) has the same x → 0 puncture structure as the cone: one main puncture (ζ → 0) and gcd(N−1, r) companion punctures (ζ → −1), with cone-determined ramifications and leading Puiseux data.
- **Lemma 3 (Residue localization).** After kernel normalization, w restricted to a generic fiber has pole order exactly k in x at every puncture, with x^(−k)-coefficient φ(0) at the main puncture and φ(−1) at every companion puncture — independent of the tail.
- **Theorem B (Uniform Residue Theorem).** φ(0) = −λ/k and φ(−1) = +λ/((N−1)k), so the residue ratio is −1/(N−1) for every k and r. For N ≥ 3 the master identity then forces μ = 0 at every degree D ≥ 1. **No scaffold exists over any fertile split-edge cell with N ≥ 3, for any polynomial tail.** At N = 2 the ratio is −1 and the argument degenerates exactly on the odd-D stratum — which is T7's jewel territory, already closed by Parity-Pole.

Combined with T5 (dominant edges), T6/No-Conspiracy (staircases, modulo the recursion mop-up), T7 (N = 2), and T8/Lemma 1 (which cells exist at all): every split-edge cell is closed. The Line Theorem's remaining debt is the T6 conspiracy recursion only.

## 1. Theorem A: the Shear Orbit

**Statement.** At (N,k,r) = (3,1,3) let c = γx + Y + x³(Y+γx)³ and w = −λ(x^(−1) + (3/2)x²(Y+γx)²). Then J(w, c) = λx^(−2) identically; in plane coordinates (Y = x³y) the pair is c = xq(1 + x⁵q²), w = −λ(1/x + (3/2)x⁴q²) with q := γ + x²y, and J_{x,y}(w, c) = λx with u(0,y) = γ ≠ 0 and w = a/x, x ∤ a.

**Proof.** The shear S: (x, Y) ↦ (x, Y + γx) has Jacobian determinant 1 and fixes x, so J(f∘S, g∘S) = J(f, g)∘S for any f, g, and (λx^(−2))∘S = λx^(−2). Apply to the germ identity J(w_germ, c_min) = λx^(−2) (artifact A6). Expanding c_min∘S gives the tail monomials γx, 3γx⁴Y², 3γ²x⁵Y, γ³x⁶ — exactly the shields the minimal ladder forces (machine: marcher2.py --minimal --no-e1; shield values 3γ, 3γ², γ³). Plane form and legitimacy: sympy-exact (verify_residue.py checks 1–2). ∎

**Consequences.** (i) The clean minimal march is quiet from bucket 12 through 46 (machine), and quiet forever by the closed form: Lemma E-minimal ("c_j ≠ 0 for all j") is **false**; c_j = 0 for all j ≥ 4, and the handoff's |c₄| = 30 falsification target fired. (ii) Since the shear orbit consists of exact pairs with nonconstant u, no argument at the pair level (obstruction ladders, tail exhaustion, terminal ideals) can close the cell. The closure must use the P-layer, and does, below. (iii) The shear is one member of the gauge pseudo-group (x,Y) ↦ (x, Y + h(x)), h(0) = 0: these all preserve the pair equation. The minimal ladder is the γx-orbit; the ladder's "shields" were the orbit's Taylor coefficients all along.

## 2. The line calculus (tool for everything below)

For a single weight-line based at x^a Y^b, write the line as x^a Y^b φ(ζ). The bracket against c_min is (machine: fertility_residue_sweep.py preamble; derivation is two lines from the monomial bracket):

    J(x^a Y^b φ(ζ), Y + x^r Y^N) = x^(a−1) Y^b [ (a + (Na − rb)ζ) φ + rζ(1+ζ) φ′ ]     (†)

Three structural facts, used repeatedly:

1. On monomials, (†) is bidiagonal: L[ζ^n] = (a + rn) ζ^n + (Na − rb + rn) ζ^(n+1).
2. The φ′-coefficient rζ(1+ζ) vanishes at ζ = 0 and ζ = −1 — one zero per component of the cone fiber {c_min = 0} = {Y = 0} ∪ {1 + ζ = 0}. Any identity L[φ] = ρ therefore **evaluates itself** at these two points: a·φ(0) = ρ(0) and (a − (Na − rb))·φ(−1) = ρ(−1).
3. The kernel of J(·, c_min) on a single line is nonzero only on the lines through (0, β), where it is exactly ℂ·c_min^β (proof: with a = 0, b = β the recursion from fact 1 gives φ_n = C(β, n)φ₀, i.e. φ = φ₀(1+ζ)^β, and Y^β(1+ζ)^β = c_min^β; on a ≠ 0 lines the ζ^0-equation forces φ₀ = 0 and the recursion kills everything unless the diagonal a + rn vanishes at some n ≥ 1, which relocates the line to a (0, β) line).

## 3. Lemma 1: fertility and germ uniqueness

**Statement.** The germ equation J(x^(−k) φ(ζ), c_min) = λ x^(−1−k), λ ≠ 0, has a polynomial solution iff r | Nk and r ∤ k. When it does, φ is unique, of degree Nk/r, given by φ₀ = −λ/k and (a + rn)φ_n = −(Na + r(n−1))φ_(n−1) with a = −k.

**Proof.** By (†) with a = −k, b = 0 the equation is the recursion above with source λ at ζ^0. Termination: the super-diagonal coefficient Na + r(n−1) vanishes iff n = 1 + Nk/r, an integer iff r | Nk; past that index all φ_n = 0, giving degree Nk/r. Obstruction: the diagonal a + rn vanishes iff n = k/r; if r | k this happens at a positive index where the super-diagonal coefficient equals −(N−1)k − r ≠ 0, forcing φ_(n−1) = 0 and, chasing back down, φ₀ = 0, hence λ = −kφ₀ = 0, contradiction. If r ∤ k the diagonal never vanishes and the recursion determines every φ_n from φ₀: existence and uniqueness. Uniqueness against kernel contamination on the germ line: by fact 3 of §2 the (−k, 0)-line carries no kernel (it contains no (0, β) monomial precisely because r ∤ k). ∎

Machine: 15 fertile germs across N = 3, 4 verified exact; criterion matches the shipped cone_solver on a 19-cell grid and reproduces the N = 2 jewel-ladder law (fertility_residue_sweep.py, 17/17). Note this upgrades T9 from window-certified to proved, and shows fertile cells are infinite in number (all k with r | Nk, r ∤ k) — including N ≥ 4.

## 4. Lemma 2: hull rigidity

**Setting.** Admissible tails are the chart monomials x^a Y^b, 0 ≤ b ≤ N, other than the excluded quartet (0,0), (1,0)=γ-slot, (0,1), (r,N), that keep the Newton data of the cell: strictly above the first edge and strictly above the c_min edge. [Rev-2 coordinate clarification: the chart monomial x^a Y^b is the c-monomial x^{a+σb} y^b, equivalently the u-monomial x^{a+σb−1} y^b (c = xu); the fiber-hull vertex from Y is (σ, 1) in c-coordinates, (σ−1, 1) in u-coordinates. The displayed tests below are in u-coordinates; the resulting inequalities were machine-swept as exactly equivalent to the plane-hull conditions, so they are unaffected.] In plane coordinates (u-monomial x^(a+σb−1) y^b) the two conditions are:

- **Edge 1** (from the constant vertex to the (σ−1, 1) vertex): strictly above ⟺ b(σ−1) < a + σb − 1 ⟺ **a + b ≥ 2**. The only lattice points failing this are the excluded (0,0), (1,0), (0,1).
- **Edge 2** (the c_min edge, chart points (0,1)–(r,N)): the cross product of (a, b−1) against the edge direction (r, N−1) is (N−1)a − r(b−1) = **v + r**, so strictly above ⟺ v > −r. σ-free.

(Machine: exact equivalence of these inequalities with the plane-hull conditions swept over cells (3,3),(3,6),(3,9),(4,4),(4,8) × σ = 2..11, full lattice.)

**Statement.** For every admissible tail and generic t, the closure of the fiber {c = t} meets x = 0 in exactly 1 + gcd(N−1, r) punctures: the **main** puncture, unramified, on which ζ → 0 and y ~ t x^(−σ); and d := gcd(N−1, r) **companion** punctures, each of ramification e := (N−1)/d, on which ζ → −1. The leading Puiseux coefficients at each puncture solve the cone's edge equations and are tail-independent.

**Proof.** The x → 0 branches of {c = t} are governed by the lower hull of the Newton polygon of c − t. The polygon's lower hull is the two-edge cone hull: the constant vertex (0,0) (from −t), the vertex (σ−1, 1) (from Y), and the vertex of x^r Y^N; every tail monomial lies strictly above both edges by the two inequalities, so the hull — hence the branch data — is that of c_min + γx − t. Edge 1 has lattice width 1: one unramified branch, balance γx + x^σ y·(1) = t + o(, i.e. y ~ t x^(−σ); on it ζ ~ t^(N−1) x^r → 0. Edge 2's balance is Y + x^r Y^N = 0, i.e. 1 + ζ = 0: Y^(N−1) = −x^(−r), which splits into d Puiseux cycles of ramification e each; on every one, ζ → −1 identically (the balance equation is ζ = −1). Leading coefficients solve these two balances, which contain no tail data. ∎

**Remark (gap in T11 as stated).** Admissibility requires only v > −r, so negative-weight slots exist — at (3,1,3) e.g. (1,1) (weight −1, plane x³y), strictly above both edges and not in the shear gauge orbit. The handoff's knob lattice (weights 0..12) omitted every −r < v < 0 slot, so the T11 ladder and any exhaustion argument built on it was not fully general. Nothing below depends on T11; the present argument covers all v > −r tails.

## 5. Lemma 3: residue localization

**Statement.** Let (w, c) be any pair in the fertile cell (N, k, r) (any admissible polynomial tail). Then there is a polynomial g such that w − g(c) restricted to the generic fiber {c = t} has, at every x → 0 puncture, pole order exactly k in x with x^(−k)-coefficient

    main puncture: φ(0)      every companion puncture: φ(−1),

where φ is the germ polynomial of Lemma 1. Since g(c) = g(t) is constant on the fiber, the master identity's puncture matching (§6) is unaffected by g.

**Proof.** w has finite chart support (w = a/x^d, a polynomial). Grade by v. (i) *Kernel normalization.* [Rev-2: constants corrected to the general shift. The bracket of monomials adds weights up to the constant shift r − (N−1): v(J(M₁,M₂)) = v₁ + v₂ + r − (N−1). The rev-1 display used "+1", the r = N specialization.] Let v* be the minimal line of w. Admissible tails satisfy v > −r, hence (integer weights) v ≥ −r + 1, and the γ-slot has v = N−1 > 0; so for any line v₁ of w and any tail monomial v₂ ≥ −r + 1, the bracket contributes at v₁ + v₂ + r − (N−1) ≥ v₁ − (N−1) + 1 > v₁ − (N−1), while the c_min-bracket of the line v₁ (c_min has v = −r) contributes at exactly v₁ − (N−1). Hence the minimal bucket of J(w, c) is v* − (N−1), achieved only by J(w_{v*}, c_min). If v* < −(N−1)k, then v* − (N−1) < −(N−1)(k+1) = v(λx^{−1−k}) = v(RHS), so L[w_{v*}] = 0, and by §2 fact 3, w_{v*} = const·c_min^ν (v* = −rν).

*Coset refinement* [rev-2; repairs the audit gap at d := gcd(N−1, r) > 1]. On a fixed weight level, two monomials (a, b), (a′, b′) satisfy (N−1)(a − a′) = r(b − b′), so b-steps along a level are multiples of (N−1)/d — while multiplication by ζ = x^r Y^{N−1} steps b by N−1, reaching only one of the d congruence classes of b within the level lattice. Each level therefore decomposes into d ζ-chains ("cosets"), and the rev-1 argument implicitly treated the level as a single chain. The decomposition is preserved by the operator: bracketing against Y fixes b, and against x^r Y^N shifts b by N−1 ≡ 0 modulo the coset index — so J(·, c_min) maps each coset to itself, and the §2 line calculus runs coset-by-coset with the same bidiagonal structure. The source λx^{−1−k} and the germ line x^{−k}φ(ζ) lie in the b ≡ 0 coset; on every other coset the recursion is homogeneous, §2 fact 3's kernel analysis applies verbatim per coset (kernel only on the (0, β) lines, absorbed into g(c)), and the coset solves to zero at and below the germ level. Sharper kernel-exclusion at the germ grade itself [GPT-Pro's arithmetic, adopted rev-2.1]: a kernel element at the germ grade −(N−1)k would be a c_min-power c_min^ν with −rν = −(N−1)k, i.e. r | (N−1)k; fertility gives r | Nk, and together these force r | k — contradicting the fertility condition r ∤ k. So no coset, including the germ's own, carries kernel at the germ grade: the germ is the unique inhabitant of its grade outright. All subsequent steps ((ii) germ minimality, (iii) puncture evaluation) then concern the b ≡ 0 coset only, as before. Machine witness at the fertile d = 2 cell (3, 2, 6): verify_326_coset.py 8/8 — germ w = −(λ/2)x^{−2}(1 + (3/2)ζ) exact; coset preservation; an odd-coset-only windowed solve forces λ = 0 (its kernel is exactly the odd c_min-powers); Lemma 2's count gives 1 + d = 3 punctures with the d = 2 companions unramified (e = (N−1)/d = 1), each carrying φ(−1) = +λ/4 against φ(0) = −λ/2 — ratio −1/2 = −1/(N−1); trace identity 1·φ(0) + (N−1)·φ(−1) = 0; D = 2 puncture matching forces μ = 0. Replace w by w − const·c^ν: this cancels the line (c^ν = c_min^ν + strictly-higher-v terms) and changes w only by a polynomial in c. Iterate (finitely, by finite support): all lines below −(N−1)k are absorbed into g(c). (ii) *Germ minimality.* Now v* = −(N−1)k (the RHS bucket must be sourced) and the line equation is the germ equation: by Lemma 1's uniqueness, w_{v*} = x^(−k) φ(ζ). (iii) *Puncture evaluation.* At the main puncture, a chart monomial x^a Y^b sizes as t^b x^a; on the germ line only (−k, 0) has a = −k, and any other line's monomials have a ≥ −k+1 (a monomial (−k, b), b ≥ 1 has v = −(N−1)k − rb < v*, excluded). So the x^(−k)-coefficient is φ(0) exactly. At a companion puncture with local parameter s, x = s^e: a monomial on the line v sizes as s^(ev/(N−1)) with tail-independent leading data (Lemma 2); the germ line gives s^(−ek)·φ(−1) (ζ → −1), and every line v > −(N−1)k contributes only above s^(−ek). So the s^(−ek) = x^(−k)-coefficient is φ(−1) exactly. Pole orders are exactly k because φ(0), φ(−1) ≠ 0 (Theorem B computes them). ∎

## 6. Theorem B: the Uniform Residue Theorem

**Statement.** For every fertile split-edge cell with N ≥ 3, every admissible polynomial tail, and every D ≥ 1: no polynomial P with deg_W P = D satisfies the master identity with κ ≠ 0. Hence **no scaffold exists over any such cell**, and by the Line Dichotomy no counterexample to JC2 arises from split-edge cells with N ≥ 3.

**Proof.** By §2 fact 2 applied to the germ equation (source ρ = λ, constant in ζ): a·φ(0) = λ and (a − Na)·φ(−1) = λ with a = −k, i.e.

    φ(0) = −λ/k,     φ(−1) = +λ/((N−1)k),     ratio = −1/(N−1).

Write P_W(W, C) = Σ_{j≤D−1} (j+1) p_{j+1}(C) W^j, and let m be the pair's exponent (m = σ − 1 − k for the realization at hand; no relation between m and D is assumed — T4's quantization is not used). Restrict the master identity P_W(w, c) = μ x^(−m) to a generic fiber and normalize w by Lemma 3's g (which only reshuffles the p_j(t) by polynomials in t — absorb and rename). At a puncture of ramification e, w ~ ρ_p s^(−ke)(1 + o(1)) with ρ_p ∈ {φ(0), φ(−1)}, both nonzero, so the W^j term sits at order s^(−jke) with leading coefficient (j+1) p_{j+1}(t) ρ_p^j, and these orders are strictly separated in j; the right side is μ s^(−me).

Induct down from the top degree. Let J ≥ 1 be the largest index with p_{J+1} ≢ 0. The s^(−Jke)-coefficient of the identity at each puncture reads

    (J+1) p_{J+1}(t) ρ_p^J = (μ if Jk = m, else 0)   for generic t.

If Jk ≠ m this forces p_{J+1} ≡ 0 (ρ_p ≠ 0), contradiction. If Jk = m, comparing the main puncture with a companion puncture gives φ(0)^J = φ(−1)^J, i.e. (−1/(N−1))^J = 1 — impossible for N ≥ 3, J ≥ 1 since the modulus is (N−1)^(−J) < 1 — again forcing p_{J+1} ≡ 0. So every W-coefficient above degree 0 vanishes and P_W = p₁(C). Terminal cases in m [rev-2.1: the m < 0 case was omitted in rev-1 — GPT-Pro audit]: if m ≠ 0 (either sign): restricting p₁(c) = μx^(−m) to a generic fiber gives p₁(t) = μx^(−m) with the left side constant in x and x nonconstant on the fiber, so μ = 0. If m = 0: p₁(c) = μ forces p₁ ≡ μ constant, and f = P(w, c) = μw + p₀(c) is a polynomial only if μw is — but w has a pole of order k ≥ 1 on Γ, so μ = 0. In every case μ = 0, i.e. κ = 0: not a scaffold, for EVERY m ∈ ℤ — the claimed independence from T4's m ≥ 1 is now actually delivered by the written proof. (T3's independent D = 1 closure and T4's quantization are consistent corroboration, not inputs.) ∎

## 6.5 Theorem C: the Trace Lemma (found after the above; simplifies and generalizes the D = 2 case)

**Statement.** Let (w, c) be ANY line pair (poles of w exactly on Γ = {x = 0}, J(w, c) = λx^m, m ≥ 1) with deg_y c = N ≥ 2. Then the fiberwise trace Tr_c(w)(x, t) := Σ_{i} w(x, y_i(x, t)) over the N sheets of c = t is **independent of x**: Tr_c(w) = τ(t). Consequently no scaffold with deg_W P = 2 exists over any such pair.

**Proof.** Along a fiber, dy_i/dx = −c_x/c_y, so (d/dx) w(x, y_i) = (w_x c_y − w_y c_x)/c_y = J(w, c)/c_y = λx^m/c_y. Summing over sheets: (d/dx) Tr = λx^m · Σ_i 1/c_y(x, y_i). For fixed x, the y_i are the roots of p(y) := c(x, y) − t, and Σ_i 1/p′(y_i) is the sum of finite residues of dy/p(y), which vanishes because the residue at y = ∞ is zero for deg p = N ≥ 2. So Tr is locally constant in x on the generic fiber, hence equal to a function τ(t); since Tr is algebraic over ℂ(x, t) and x-free, τ ∈ ℂ(t). [Rev-2.1: the rev-1 parenthetical claiming τ polynomial "for our pairs" was not established and is not needed — the D = 2 contradiction below only uses x-freeness. Deleted per the GPT-Pro audit.] For D = 2 the master identity reads 2p₂(c)w + p₁(c) = μx^(−m); trace it: 2p₂(t)τ(t) + N·p₁(t) = N·μ·x^(−m). The left side is x-free and the right is not (m ≥ 1), so μ = 0. ∎

Machine: on the sheared (3,1,3) pair, Tr(w) = 0 identically (symmetric-function computation); Tr(w²) is NOT x-free, exactly as it must be — for j ≥ 2, Σ w^(j−1)/c_y picks up a residue at y = ∞ once (j−1)·deg_y(numerator of w) exceeds N − 2, so the trace mechanism stops at D = 2. That boundary is real: A3's genuine P_W = μ(W² − 4C) (a W-degree-2 P_W hitting x^(−2)) lives exactly on the other side.

**What the trace explains.** (i) The residue ratio: boundedness of Tr as x → 0 forces the ramification-weighted leading coefficients to cancel, 1·φ(0) + (N−1)·φ(−1) = 0, which IS the ratio −1/(N−1). The magic number is just "the trace is blind to x," i.e. Galois symmetrization over the fibration — the Galois-collapse guardian from T2 in yet another costume. (ii) Scope: Theorem C kills every D = 2, m ≥ 1 scaffold over EVERY line pair with deg_y c ≥ 2 — no cell structure, no hull rigidity, no Puiseux needed — subsuming the D = 2 case of Theorem B and conceptually reproving A3's "P_W = μ/x machine-empty." Higher D still needs Theorem B's puncture induction (on N ≥ 3 cells) or T7 (N = 2 odd D).

**Boundary remark (N = 2).** The ratio is −1/(N−1) = −1, and (−1)^(D−1) = 1 exactly when D is odd. So for N = 2 the residue clash kills all even D — reproving T7's parity half in one line — and is silent on odd D, where the jewel pairs genuinely live and T7's Parity-Pole + 1/(2I+1) arguments were required. The two mechanisms agree on their overlap and the boundary is exactly where it must be. This is a strong consistency check of both.

## 7. What this closes and what remains

Closed by this document (audit pending): the entire split-edge N ≥ 3 arc — all fertile cells (infinitely many, classified by Lemma 1), all polynomial tails, all P-degrees. Together with T5 (dominant edges), T7 (N = 2), T8 (which (N, r) exist), and T3/T4: every cell-type stratum of the Line Theorem.

Remaining for the Line Theorem: the T6 conspiracy recursion mop-up (≥ 3-edge staircases with all-multiple deeper-edge roots) — separate machinery, quantified deficit already proven, recursion analysis not started. Then assembly + writeup.

**Audit requests (specific).** (1) §5 step (i): the bucket bound uses v(tail) > −r; check the inequality chain once more against the γ-slot (v = N−1) and negative-weight tails. (2) §6: the leading-term matching at each puncture — confirm no lower-order term of P_W(w, c) or subleading term of w can reach the s^(−k(D−1)e) coefficient (orders are strictly separated; stated, should be checked). (3) §4: genericity of t (finitely many bad fibers; identities of polynomials in t). (4) §6 no longer uses T4's quantization (the induction handles every (D, m) via order mismatch or clash) — confirm the removal is sound, i.e. that the s^(−Jke)-coefficient extraction is exhaustive: lower W-powers contribute at −jke + i, i ≥ 0, and −jke + i = −Jke would need i = −(J−j)ke < 0, impossible. (5) Lemma 2's d = gcd(N−1, r) cycle count — cross-check at (3,1,3): d = 1, one companion puncture of ramification 2 (the ±i sheets are one puncture; the residue computation used both sheets and got equal values, consistent).

## 8. Machine witnesses

- verify_artifacts.py — 13/13 (A1–A6, T10).
- marcher_repro.py + compare_repro.py + shipped bucket_marcher.py — T11 triple-reproduction, canonical spine −3, −9, +18.
- marcher2.py --minimal --no-e1 — clean minimal ladder: quiet 12..46, shields = shear binomials.
- verify_residue.py — 22/22: exact pair (chart + plane), Puiseux branches and residues −λ / +λ/2, clash μ = 0, hull-rigidity sweep, w-support knob/E1-freeness, windowed P_W emptiness.
- fertility_residue_sweep.py — 17/17: 15 exact germs, criterion ⟺ cone_solver on the grid, N = 2 jewel-law control, hull-inequality equivalence sweep.
- verify_residue_23.py — 15/15: (2,3) explicit instance (chart + plane pair, Puiseux residues −λ/2 / +λ/4, windowed P_W empty, Tr(w) = 0, hulls to N = 6).
- verify_326_coset.py — 8/8 [rev-2]: the fertile gcd(N−1, r) = 2 cell (3,2,6) — germ exact, coset preservation, odd-coset emptiness, unramified companion pair with uniform residues, ratio −1/2, trace identity, D = 2 clash.
- GPT-Pro adversarial battery [rev-2.1, external, reviews/gpt-pro/]: two NEW exact fertile-cell pairs built independently — (4,1,2) at σ = 3 (germ φ = −1 − 4z − (8/3)z², φ(0) = −1, φ(−1) = +1/3, ratio −1/3 = −1/(N−1) ✓; reducible plane c; one companion cycle of ramification 3) and (3,2,6) at σ = 4 (germ φ = −1/2 − (3/4)z — EXACT agreement with verify_326_coset.py's independently computed germ); 60 exact rational span searches for a master identity on both: empty. Plus the staircase batteries: 320 targeted + 540 deterministic-random exact systems incl. m < 0, all barren.
- Robustness runs: full-lattice march (--min-knob-wt -2: new bucket −2 obstruction, γ-spine intact, no terminals — the T11 gap is material but closure-neutral); E₂ probe (E2·c_min² genuinely obstructs at bucket −1 — not a kernel element; honest kernel g(c) is inert by J(g(c), c) = 0; the E1/E2 absorbability asymmetry is recorded as a curiosity).

All in the jc lane workspace (line-theorem/); shipped bundle + this doc + theorem_ledger_v2.json in ~/claudes-world/jc/.
