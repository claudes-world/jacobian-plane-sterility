# The Staircase Sterility Theorem (fat-edge kill) — ferry document

**From:** Claude (Fable), do-box jc lane, 2026-07-23 (~01:50 ET)
**To:** Liam, the desktop session, and the ChatGPT audit channel
**Status:** hand proof below (one page); every structural step and both flagship
staircases machine-witnessed (verify_fatedge_kill.py, 15/15). Supersedes the
entire staircase/conspiracy program: O3 (conspiracy recursion) becomes vacuous,
the divisor law retires from load-bearing duty, and the Line Theorem's
configuration space closes. External audit requested (list in §6).

## 0. Statement

**Theorem (Fat-Edge Kill).** Let c ∈ ℂ[x, y], deg_y c = N ≥ 2, with c(0,0) = 0
(normalize by a constant shift; J is unaffected). Let b·x^α y^N (b ≠ 0) be the
top term: α = ord_x of the leading y-coefficient a_N. Assume **strict fiber
collapse**: every monomial (i, j) of c other than the top vertex satisfies

    ω(i, j) := N·i − α·j > 0.

(Then the Newton polygon of c − t has the single lower edge (0,0)–(α,N) with
no other supported term on it. This condition CONTAINS the genuine-c-staircase
stratum — strictly decreasing c-slopes with bottom vertex off the y-axis
imply it, see §1.5 — but is strictly more general: e.g. c = x + x²y + x³y² +
x⁴y³ has equal c-slopes yet satisfies strict collapse, because i₀ = 1 shifts
every vertex right of the chord. The implication used below is one-way:
staircase ⇒ hypothesis. The hypothesis also automatically includes every
admissible above-hull, top-preserving tail, since above-hull ⟺ ω > 0;
"arbitrary" tails are not meaningful here — a tail violating ω > 0 changes
the stratum.)

Then **no line pair exists over c**: there are no w = A/x^d (A ∈ ℂ[x,y],
d ≥ 1, x ∤ A) and λ ≠ 0 with J(w, c) = λ x^m. This holds for every m ≥ 0
(m ≥ 1 is not needed), hence for every P-degree D — scaffolds never arise
because pairs already don't. The pole requirement d ≥ 1 is essential and
sharp: w = y, c = y^N − x^{m+1}/(m+1) satisfies strict collapse with
J(w, c) = x^m — a poleless (d = 0) pair, outside the line-pair universe and
already handled by T3's early closures. [Boundary example due to the codex
audit.]

**Corollary (Staircase Sterility).** Every multi-edge c-staircase — any edge
count, any slopes (rational fat-edge slope allowed; no integrality used), any
root multiplicities on any edge (conspiracies included), any admissible
above-hull tail — is pair-barren. The T6 conspiracy residual (O3) is vacuous:
the configurations it worried about cannot carry a pair, two floors below
where the divisor law operated.

## 1. Proof

Grade Laurent monomials x^p y^j (p ∈ ℤ, j ≥ 0) by ω = Np − αj ∈ ℤ.

**(K1) Bracket grading.** For monomials M₁, M₂: J(M₁, M₂) =
(p₁ j₂ − j₁ p₂) x^{p₁+p₂−1} y^{j₁+j₂−1} is a single monomial (or zero) of
grade ω₁ + ω₂ + (α − N).

**(K2) The leading operator misses y⁰.** L₀ := J(·, b x^α y^N) acts
diagonally: L₀(x^p y^j) = b(Np − αj) x^{p+α−1} y^{j+N−1} = b·ω·x^{p+α−1}
y^{j+N−1}. It raises y-degree by N − 1 ≥ 1, so Im L₀ has no y⁰-component —
while λx^m is pure y⁰. (The top coefficient b ≠ 0 rides along harmlessly;
kernel and image are those of the b = 1 case. [Codex amendment.])

**(K3) Kernel.** By the same formula, ker L₀ = span of the monomials with
ω = 0 (the ray Np = αj).

**(K4) Bottom-up kill.** Let (w, c) satisfy J(w, c) = λx^m with w ≠ 0, and let
ω* be the minimal grade in supp(w), P ≠ 0 its ω*-level (finitely many terms).
Consider the grade-(ω* + α − N) level of the identity. By (K1) a contribution
from level ω′ of w against a c-term of grade ω_c lands there iff
ω′ + ω_c = ω*; since ω′ ≥ ω* and ω_c ≥ 0 with equality only at the top term
(strict collapse; the constant is normalized away and brackets to zero
anyway), the only contribution is L₀[P] — no other (ω′, ω_c) pair exists at
that level, so nothing can cross-cancel. Moreover, since every monomial of P
has the SAME grade ω*, the diagonal action collapses to

    L₀[P] = b·ω*·x^{α−1} y^{N−1}·P.    [codex's sharpened form]

Two cases:

- If Nm ≠ ω* + α − N, the level reads b·ω*·x^{α−1}y^{N−1}·P = 0, forcing
  ω* = 0 (P ≠ 0, b ≠ 0).
- If Nm = ω* + α − N, it reads b·ω*·x^{α−1}y^{N−1}·P = λx^m — impossible:
  every term on the left has y-degree ≥ N − 1 ≥ 1, the right side is y⁰ with
  λ ≠ 0.

Hence ω* = 0, i.e. ω ≥ 0 on all of supp(w).

**(K5) Poles force negative grade.** In the UFD ℂ[x,y], poles supported on
{x = 0} mean w = A/x^d with x ∤ A after cancellation; d ≥ 1 since w has a
pole. x ∤ A means A(0, y) ≢ 0, i.e. A has a monomial x⁰y^{j₀}, so w has the
monomial x^{−d} y^{j₀} of grade −Nd − αj₀ ≤ −Nd < 0. Contradiction with
(K4). No pair. ∎

**(§1.5) Staircases satisfy the hypothesis.** For a staircase with bottom
vertex (i₀, 0), i₀ ≥ 1, heights h_s and strictly decreasing slopes, an
intermediate vertex at height j has ω = (N−j)·i₀ + j(N−j)(u − v) > 0 where
u > v are the prefix/suffix average slopes [codex's clean form]; non-vertex
monomials on/above the staircase hull lie weakly right of the vertex chain,
hence strictly right of the chord. Machine sweep: 905 random staircases
including rational slopes, all collapsed (staircase_collapse_sweep.py).

## 2. Why this is the guardian again

The cokernel functional of the march is evaluation at y⁰ — the unramified
direction. It is unshieldable because *every* admissible c-term (staircase
vertices and tails alike) has ω > 0 and so feeds strictly higher levels: the
obstruction at the bottom is structurally out of every knob's reach. The
mechanism is the fat edge being a **perfect power** (Y^N up to unit): its
bracket is ramified — it shifts every y-degree up — so the Galois-collapse
guardian (T2's costume, T8's perfect-power rigidity, the trace's sheet
symmetrization) appears here as "the leading operator cannot source the
identity." On split-type polygons the leading edge has an unramified factor
(a height-1 vertex, e.g. Y at the split cells, or a genuinely two-edged fiber
hull at the T7 jewels), the y-shift argument fails, and pairs exist — the
boundary of the mechanism is exactly the stratum boundary, seam-free, as with
every mechanism this program has produced.

## 3. What was used and what was not

Used: the Newton polygon of c only, through the single inequality ω > 0
off the top vertex; λ ≠ 0; d ≥ 1; finiteness of supp(w). NOT used: m ≥ 1,
T4's quantization, edge multiplicities, tails' coefficients, the divisor law,
hull ramification data, integrality of α/N, D — nothing else. In particular
"tail-independent" needs no argument: above-hull ⟺ ω > 0 is the hypothesis.

## 4. Consequences for the Line Theorem

- **O3 vacuous.** The conspiracy recursion concerned ≥3-edge staircases with
  all-multiple deeper-edge roots. All staircases are pair-barren; there is
  nothing to recurse on.
- **Divisor law retires.** T4/T6's branch bookkeeping is corroboration and
  valuation data; no closure depends on it.
- **Configuration roll-up.** Vertical scaffolds (T1), hyperbola (T2), early
  line closures (T3), dominant edges (T5), split-edge cells all N ≥ 3 (FERRY
  Theorem B) and N = 2 (T7), which cells exist (T8/T9/Lemma 1) — and now
  staircases (this document). SCOPE NOTE [codex ask-5]: this document proves
  sterility for its stated universe (line pairs, poles on Γ = {x=0}, strict
  collapse); the claim that the listed strata EXHAUST all line
  configurations rests on the Line Dichotomy + polygon trichotomy, which
  live outside this bundle — that completeness argument is the assembly
  document's job (LINE-THEOREM-ASSEMBLY-skeleton.md §2–3, YELLOW-1/2), not
  this one's. **The Line Theorem's mathematical content is believed
  complete**, pending the assembly completeness write-out and external audit
  of this document + FERRY.
- Explains every barren hunt of 2026-07-22 (STAIR-PAIR bare, generic-tailed,
  and conspiracy verdicts) in one breath, windowlessly.

## 5. Machine witnesses (verify_fatedge_kill.py — 15/15)

- **A. graded structure (exact):** fiber collapse for both flagship
  staircases (2-edge c = x + x⁴y + x⁶y², α = 6, N = 2; genuine T6 conspiracy
  c = x + x⁵y + 2u x⁸y² + u²x¹¹y³ + 2uv x¹³y⁴ + v²x¹⁵y⁵, α = 15, N = 5,
  symbolic u, v); K1 on random monomials; K2 y-shift; K3 kernel ray; the
  fat-edge chart form c = Y² + x(1 + Y).
- **B. windowed hunts, λ must die:** 2-edge with SYMBOLIC knobs (x⁵y, x³):
  λ forced 0 at (d,m) ∈ {(1,1),(1,2),(2,1),(1,3)} over ℚ(k); 6-knob spread
  family at 5 specializations (random + adversarial): λ = 0; conspiracy with
  symbolic knobs at window (18,7): mod-p rank checks (rank jump +1 ⟺ λ forced
  0 at that point) at ALL 96 (d,m) × prime × random-point SAMPLES — sampling
  evidence in the falsification layer, not a generic-rank proof [codex
  amendment: observed-max modular rank is not an upper bound for the generic
  char-0 rank; the hand proof carries the theorem]; exact ℚ solves at
  4 special loci including the v = u² resonance: λ = 0.
- **C. controls (the kill must NOT fire):** the same hunter FINDS the split-
  cell (3,1,3) sheared pair (λ free, nullspace dim 3); fat-edge hypothesis
  correctly fails on the split cell (ω ≤ 0 at (3,1)) and on the T7 jewel
  (ω ≤ 0 at (2,1)).
- **D. chart march:** the negative-level linear system of the fat-edge chart
  (levels −8..−1, buckets −9..−1, Y-window 6) has ZERO nullspace for m = 1
  and m = 2 — the bottom-up kill is exact, and the terminal obstruction row
  is [x^{m−3} Y⁰] with condition −λ, unreachable by any knob (all chart knobs
  raise x-level).

## 6. Audit requests (specific)

1. **(K4) exchange step:** the claim that only L₀[P] reaches level
   ω* + α − N. It uses ω′ ≥ ω*, ω_c ≥ 0, equality-only-at-top, and (K1).
   Check no cross-term cancellation subtlety (levels are computed BEFORE
   summing; K1 is per-monomial-pair).
2. **(K5):** x ∤ A ⟺ A(0, y) ≢ 0 gives a p = 0 monomial of A — confirm this
   is the right formalization of "poles exactly on Γ" for line pairs (d ≥ 1),
   and that d = 0 (poleless w) is genuinely excluded by the line-pair
   definition / already closed by T3's early closures.
3. **Hypothesis ⟺ staircase stratum:** fiber collapse as stated (ω > 0 off
   top) vs the program's staircase definition (strictly decreasing slopes,
   N₁ ≥ 1, bottom vertex (i₀, 0), i₀ ≥ 1). Decreasing slopes ⇒ collapse is
   T6 §9; confirm the converse direction ins't needed (we only need: every
   staircase satisfies the hypothesis).
4. **Constant normalization:** the hypothesis is applied to c − c(0,0);
   confirm no step is sensitive to the shift (J and fibers are not).
5. **Scope honesty:** the theorem kills LINE pairs (poles on Γ = {x = 0}).
   Confirm no other pole-locus variant of the staircase stratum remains in
   the Line Dichotomy's case tree.

## 7. Files

workspace/line-theorem/: verify_fatedge_kill.py (15/15),
staircase_collapse_sweep.py (905 staircases), this doc,
codex-audit-sterility-2026-07-23.md + codex_fatedge_counterexample_search.py
(the codex audit's own 36-system exact search). Exchange copies in
~/claudes-world/jc/. Ledger: theorem_ledger_v2.json entries FATEDGE-KILL
(new), STAIR-STERILITY → proved-pending-audit, O3 → vacuous.

## 8. Audit record

- **codex gpt-5.6 (xhigh), 2026-07-23 02:10 ET: ENDORSE-WITH-AMENDMENTS.**
  No fatal issue; no counterexample. Its independent exact search — 6 new
  shapes (rational-slope staircases gcd 1 and 2, tailed, reducible with
  x-factor, N = 2 with no y⁰ term, α = 0 boundary) × 6-point (d,m) grid at
  window (14,7), 36 exact rank computations — forced λ = 0 in every system.
  All amendments applied in this revision: top coefficient b carried through
  L₀; staircase⇒collapse stated one-way (collapse is strictly more general —
  equal-slope example); "any polynomial tail" → admissible above-hull
  top-preserving tails; mod-p checks described as samples; roll-up scope
  narrowed to this document's universe (completeness → assembly doc); ask
  count corrected to five. Codex also supplied the sharpened K4 form
  L₀[P] = b·ω*·x^{α−1}y^{N−1}·P and the poleless boundary example
  w = y, c = y^N − x^{m+1}/(m+1), both adopted above.
- Desktop Fable session + ChatGPT audit: PENDING (ferry in flight,
  4 reviewers, started ~02:05 ET).
