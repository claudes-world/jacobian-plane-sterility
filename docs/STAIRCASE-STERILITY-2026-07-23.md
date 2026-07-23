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
the stratum. [Rev-4, orientation pinned by sizing rather than by the word
"above": on the fat-edge branches y ~ s·x^{−α/N}, a monomial x^i y^j sizes
as x^{ω/N}, so ω > 0 monomials vanish relative to the edge as x → 0 —
harmless tails, right of the chord — while an ω < 0 monomial dominates the
balance, creates a lower-hull vertex below the chord, and moves c to a
different stratum of the case tree. See reviews/gemini/REBUTTAL.md.])

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

[Rev-4, one-line form — due to the chat-Claude audit: no case split is
needed. The level identity equates b·ω*·x^{α−1}y^{N−1}·P, whose every term
has y-degree ≥ N − 1 ≥ 1, with the grade-matching piece of λx^m, which is
pure y⁰. Disjoint y-supports force both sides to vanish separately; since
the monomial shift (p, j) ↦ (p + α − 1, j + N − 1) is injective and P ≠ 0,
this gives ω* = 0 directly — λ plays no role until K5.]

Hence ω* = 0, i.e. ω ≥ 0 on all of supp(w). (The original two-case split —
grade mismatch vs y-degree clash — is preserved in rev-2/3 history; the
one-line argument subsumes both.) Note the minimal grade ω* exists because
supp(w) is FINITE — for α > 0 the grade is unbounded below on y-power
series, so any future formal-series variant must re-prove the kill, not
cite it [chat-Claude ding (ii)].

**(K5) Poles force negative grade.** In the UFD ℂ[x,y], poles supported on
{x = 0} mean w = A/x^d with x ∤ A after cancellation; d ≥ 1 since w has a
pole. x ∤ A means A(0, y) ≢ 0, i.e. A has a monomial x⁰y^{j₀}, so w has the
monomial x^{−d} y^{j₀} of grade −Nd − αj₀ ≤ −Nd < 0. Contradiction with
(K4). No pair. ∎

**(§1.6) Strengthening (rev-4 chat-Claude; extended rev-5 GPT-Pro).** The
one-line K4 never uses the shape of the right side beyond "pure y⁰," so the
same five steps prove more: for c satisfying strict collapse and ANY
w = A/x^d with d ≥ 1, x ∤ A — equivalently any finite Laurent polynomial
w ∈ ℂ[x, 1/x, y] with negative x-valuation —

    J(w, c) ∉ ℂ[x, 1/x] — zero included, every m ∈ ℤ included.

(Every graded piece of a y-free Laurent RHS is a pure-y⁰ monomial; K4's
disjoint-support argument applies verbatim to each. So: no λx^m for any
m ∈ ℤ [GPT-Pro's form], no x-polynomial RHS [chat-Claude's form], and no
zero: the pole-bearing centralizer of c is empty.) Machine: centralizer
hunts (RHS ≡ 0) forced x | A (refute_fatedge.py, chat-lane); GPT-Pro's
batteries include m = −3..−1, all barren. Scope boundary kept visible:
poleless pairs (d = 0) genuinely exist over collapse shapes — J(y, y² − x)
= 1, and the exact family w = x^a y, c = x + x^{aN}y^N with J(w, c) = −x^a
[GPT-Pro] — so the d ≥ 1 restriction is load-bearing; d = 0 configurations
are outside the line-pair universe (poles exactly on Γ, which is nonempty)
and belong to T3's early-closure territory. The
staircase stratum also has N ≥ 2 by construction (top vertex strictly
above a first edge of height ≥ 1), which K2 requires; stated here so the
assembly can cite it [chat-Claude ask-5 notes].

**(§1.5) Staircases satisfy the hypothesis.** For a staircase with bottom
vertex (i₀, 0), i₀ ≥ 1, heights h_s and strictly decreasing slopes, an
intermediate vertex at height j has ω = (N−j)·i₀ + j(N−j)(u − v) > 0 where
u > v are the prefix/suffix average slopes [codex's clean form]; non-vertex
monomials on/above the staircase hull lie weakly right of the vertex chain,
hence strictly right of the chord. Machine sweep: 905 random staircases
including rational slopes, all collapsed (staircase_collapse_sweep.py).

## 1.7 Theorem′ (Ramified-Edge Kill) — post-commit generalization, rev-3

**Theorem′.** Let c be as above (deg_y c = N ≥ 2, c(0,0) = 0) with fiber
hull whose FIRST (lowest) edge has top endpoint (i₁, N₁), and assume the
first edge is **bare-ramified**: N₁ ≥ 2 and the only c-monomial on that
edge is its top endpoint b₁·x^{i₁}y^{N₁} (b₁ ≠ 0). Then no line pair
exists over c — regardless of the deeper hull structure (convex, collapsed,
or mixed), root multiplicities, or admissible tails.

**Proof delta.** Regrade by ω₁(i, j) := N₁·i − i₁·j. Every c-monomial other
than the first-edge endpoint lies strictly above the first edge's supporting
line (hull convexity + bareness), i.e. ω₁ > 0, and the endpoint has ω₁ = 0.
K1–K5 now run verbatim with L₀ = J(·, b₁ x^{i₁} y^{N₁}): the bracket adds ω₁
up to the constant shift i₁ − N₁; L₀ is diagonal with eigenvalue b₁ω₁ and
raises y-degree by N₁ − 1 ≥ 1, so it cannot source λx^m; its kernel is the
ω₁ = 0 ray; the bottom level of w is forced into the kernel, so ω₁ ≥ 0 on
supp(w); and a pole monomial x^{−d}y^{j₀} has ω₁ = −N₁d − i₁j₀ < 0. ∎

Theorem (§0) is the special case "the hull is one edge." What Theorem′ adds
to the assembly: every configuration whose first fiber edge is bare-ramified
is dead — including convex multi-edge hulls with tall first edge (the
YELLOW-1 class of the assembly skeleton) and every mixed shape. The
completeness residue shrinks to first edges with MULTIPLE supported points
(genuine edge polynomials — unramified punctures exist ⇒ cell/divisor-law
territory) and height-1 first edges (the split-edge chambers, closed).

Witness (suite section E): convex two-edge shape c = x + x²y² + x⁴y³
(N₁ = 2, slopes 1 < 2, genuinely non-collapsed: the middle vertex (2,2) has
ω = 3·2 − 4·2 = −2 ≤ 0 against the top vertex (4,3), so the §0 Theorem does
NOT apply — but under ω₁ = 2i − 2j only the first-edge endpoint (2,2) sits
at level 0, so Theorem′ does), bare and tailed: λ forced 0 at all (d, m)
tried.

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
hull at the T7 jewels), the y-shift argument fails — and pairs are KNOWN to
exist on the far side at two witnessed points (the T7 jewels; the (3,1,3)
shear cell). [Rev-4 precision, chat-Claude ding (i): this seam statement is
one-directional. What is proved: the kill's hypothesis fails exactly at the
stratum boundary. What is NOT claimed: barren ⟺ collapse — existence beyond
the boundary is governed by FERRY Lemma 1's fertility arithmetic (r | Nk ∧
r ∤ k), and some non-collapse shapes are also barren (e.g. xy + x³y²,
machine-checked in the chat-lane audit).]

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

## 5. Machine witnesses (verify_fatedge_kill.py — 18/18 as of rev-3)

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
- **Ferry round 1 verdicts (in ~02:35–02:40 ET, on rev-2):**
  - **chat-Claude (refutation attempt): FAILED TO REFUTE — strongest
    review.** ~210 exact hunts, 17 shapes, from-scratch hunter that
    blind-rediscovered the (3,1,3) shear pair and the T7 jewel pair
    (validation). Strengthened the theorem (one-line K4; J(w,c) ∉ ℂ[x],
    empty pole-bearing centralizer — adopted in rev-4 §1.6). Four dings
    (biconditional wording, finite-support note, mod-p one-directionality,
    ship-the-scripts) — all accepted, first two folded in rev-4.
    Independently found FERRY's r = N constant slip.
  - **GPT: could not refute the theorem** (356 exact systems, rational
    slopes, dense supports, d ≤ 5, m ≤ 7). Found three REAL FERRY issues:
    the r = N shift constants, the ζ-coset gap at gcd(N−1, r) > 1, and a
    Lemma-2 coordinate slip — all repaired in FERRY rev-2 with the new
    (3,2,6) witness (verify_326_coset.py 8/8).
  - **Gemini: claimed a "fatal geometry inversion" (admissible tails have
    ω < 0) — REFUTED**, see reviews/gemini/REBUTTAL.md: it reversed the
    Newton-hull orientation; the x→0 sizing (a monomial scales as x^{ω/N}
    on the fat-edge branches) fixes the direction unambiguously, both
    other reviewers independently machine-verified the exchange step, and
    the tails it predicts to enable cancellation demonstrably don't
    (every symbolic-tail hunt forces λ = 0). Its K5/constant-shift/
    rational-slope confirmations concur with the other reviewers.
  - **GPT-Pro (in ~02:59 ET, deepest FERRY read): did NOT refute the
    Fat-Edge Kill** — "the staircase core is sound and stronger than
    advertised." 320 targeted + 540 random exact systems (incl. m < 0,
    reducible/repeated factors, gcd(α,N) > 1): zero survivors. Its three
    main FERRY defects coincide with the rev-2 repairs (shift constants,
    coordinates, coset gap — its (3,2,6) germ matches verify_326_coset.py
    EXACTLY, independent derivation); two NEW items adopted in rev-2.1:
    the omitted m < 0 terminal case of Theorem B, and the sharper
    kernel-exclusion arithmetic (r | (N−1)k + fertility ⇒ r | k ⊥). It
    also built a second exact fertile cell (4,1,2) (ratio −1/3 ✓) + 60
    span searches (empty), supplied the poleless family w = x^a y, and
    its "all m ∈ ℤ" strengthening is folded into §1.6. Its residual
    caveat stands and is correct: the global Line Theorem roll-up is
    conditional on the external case partition (the assembly document's
    job), not on anything inside this bundle.
- **Rev-3 delta (post-commit, post-codex): §1.7 Theorem′ (Ramified-Edge
  Kill) + suite section E (18/18).** Theorem′ is NOT yet covered by the
  codex audit or the ferry round — include it in the next audit pass. The
  frozen hash-commit copies (frozen-20260723-0217/) remain rev-2, byte-exact;
  the public reveal will carry a rev changelog.
