# Composite-jewel cells: minimal-chart kill, the shear orbit, and the P-layer sweep

**rev-2, 2026-07-23 ~14:15 ET (post-audit)** — the E4-sliver session (assembly §5,
route 1 executed, route 3 partially executed via S-conjugation).
Witnesses: `marcher_composite.py` (runs in `out/mc_*.json`),
`verify_composite_jewel.py`, `probe_composite_cone.py`,
`classify_stratum.py`. Status: codex xhigh delta-audit FOLDED (verdict file
codex-audit-composite-2026-07-23.md; checker codexchk_composite_audit.py);
section 3 rewritten per its B/C findings — the shear-locus kill is now a
theorem for ALL (sigma, D). GPT-side audit optional per Liam.

## 0. Setting

The one open stratum of the Line Theorem partition (assembly v1.1, row E4):
two-edge split chambers whose cone edge is COMPOSITE (interior lattice
support, d = gcd(N−1, r) > 1) with distinct roots, sitting ON the fertile
locus of the Composite-Cone Fertility Lemma. First family: (N, r) = (3, 6),
d = 2, ξ = x³Y, cone polynomial E(ξ) = 1 + t·ξ + ξ² at t = ±√3
(P₁(t) ∝ t² − 3), pole depths k ≡ 1 (mod 3) — the classifier's independent
germ solver confirms fertile k ∈ {1, 4, 7} on this cell (classify_stratum
zoo line). WLOG t = +√3 (x ↦ −x maps t ↦ −t, λ ↦ −λ-compatible). Chart
data: c = γx + c_min + tails, c_min = Y + √3x³Y² + x⁶Y³, pair equation
J(w, c) = λx^{−1−k}. Grading ω(a, b) = a − 3b: cone edge at ω = −3,
J(·, c_min) graded of degree −1, γx at ω = +1, admissible tails ω > −3.
Realizations: Y = x^σ y, any σ ≥ 1; plane pair exponent m = σ − 1 − k.

Known before tonight: residue ratios at the cone punctures are primitive
cube roots of unity (Composite Residue Law φ(ρ)/φ(0) = d/((N−1)ρE′(ρ)))
⇒ the residue clash kills every scaffold with 3 ∤ J; CHI-DIV kills D ≤ 3;
windowed plane hunts with 3 knobs pair-barren.

## 1. Theorem (Minimal-Chart Kill)

**For k = 1, 4, and 7 — all three machine-confirmed fertile family
members — the minimal chart (tails = γx only, γ ≠ 0) admits no pair: the
bucket march terminates with an unshieldable cokernel obstruction at the
germ's first γ-rung, with core (k−2)·t·γ·λ (observed exactly at every
tested point; closed form CONJECTURED for the family).**

Precisely (marcher_composite.py, minimal mode, `--no-e1`):
- k = 1: the march reproduces the germ in closed form at bucket −2 —
  **w_germ = −λ(x⁻¹ + √3x²Y)** — sails through buckets −1..1 with all
  cokernel conditions vanishing, and at bucket 2 (source weight of
  J(germ, γx)) hits cokernel condition **−√3·γ·λ = 0** with NO admissible
  fresh slot (ω = 3: (3,0), (6,1), (9,2), (12,3)) entering linearly.
  With γ ≠ 0 (bottom-vertex normalization) and λ ≠ 0: contradiction.
- k = 4: germ = x⁻⁴(−λ/4 + 2√3λξ + (27λ/2)ξ² + 9√3λξ³ + (27λ/4)ξ⁴)
  found in closed form at bucket −5; terminal at bucket −1 (again the
  first γ-rung), core **+2√3·γ·λ**.
- k = 7: terminal at bucket −4 (first γ-rung again), core **+5√3·γ·λ**
  (out/mc_k7_minimal.json). Coefficient sequence −1, 2, 5 at k = 1, 4, 7:
  core = (k−2)√3γλ.
- t = −√3 twin: terminal at bucket 2, core **+√3γλ** = (k−2)·t·γ·λ at
  t = −√3 — the WLOG symmetry claim is now machine-witnessed, and the
  (k−2)·t·γ·λ form covers both signs.
- Robustness: identical spine with the E1·c_min kernel seed injected
  (no gauge leak — conditions E1-free). Deep-pole handling [wording per
  audit amendment 3]: the pre-bucket kernel recordings show the only
  in-window gauge directions are c_min-multiples, so w is NORMALIZED
  modulo ℂ[c] to pole depth exactly k in-window; the records certify
  the normalization, not a literal global vanishing statement.

Remark (structure): the terminal core is proportional to √3 = t, the
interior cone coefficient. On binomial cells (t = 0) the corresponding
γ-rung is shieldable — that is how jewels and the shear orbit exist in
E1. The interior monomial itself is the obstruction.

## 2. The forced cascade and the shear orbit

With ALL admissible tails symbolic (full-knob mode, weights −2..11,
buckets through 11): the march never terminates. Instead every bucket
β ≥ 2 fires exactly two cokernel conditions, absorbed by FORCED
impositions (unique-branch: linear coefficients are nonzero monomials
−(β−1)√3λ and −β√3λ) on the two b ∈ {0,1} tail slots of weight β−1;
the b ∈ {2,3} families stay free (23 free deformation knobs in-window
against 25 forced impositions; run mc_k1_full_deep2.json).

**The cascade IS the shear, value-exactly** (check_cascade_is_shear.py):
fixpoint-resolving all 25 forced imposition values and setting the free
knobs to 0 yields k4_1 = 2√3γ, k5_0 = √3γ², k7_2 = 3γ, k8_1 = 3γ²,
k9_0 = γ³ and ZERO on every other imposed slot — exactly the tail
monomials of c_min∘S below. Theorem (Shear Pair, exact —
verify_composite_jewel A1–A4): under S: (x, Y) ↦ (x, Y + γx),

  c_min∘S = c_min + γx + 2√3γx⁴Y + √3γ²x⁵ + 3γx⁷Y² + 3γ²x⁸Y + γ³x⁹

and (w_germ∘S, c_min∘S) is a genuine pair with the same λ, for symbolic
γ — a legitimate chart c (bottom vertex (1,0)). Verified in chart and in
the σ = 1 plane realization. So the composite-jewel cell DOES carry
pairs — the shear orbit of the germ, the exact composite analog of
FERRY's shear lemma — and the forced cascade is the shear orbit's shadow
plus residual tail freedom (the b ∈ {2,3} directions, in-window).

## 3. The Shear-Locus Kill (all σ, all D) — corrected and upgraded

**Correction (codex audit, finding B).** rev-1.1 claimed both scaffold
conditions transfer to the bare germ pair under S. That is HALF right:
the master identity transfers (S fixes x, Jacobian 1, so
P_W(w∘S, c∘S) = P_W(w, c)∘S = μx^{−m} ⟺ P_W(w, c) = μx^{−m}) — but
f-polynomiality does NOT: pole-cancellation must be imposed on the
SHEARED f∘S, and f∘S polynomial ⇎ f polynomial. Exact counterexample
(codex; reverified): P = W + W²C at σ = 2 — the bare f = P(w_germ, c_min)
carries an x⁻¹ term while the sheared f is polynomial. rev-1.1's sweep
used bare-f rows and was therefore not, as written, a certificate for
the shear locus (its numerics happened to agree; see witnesses below).

**Theorem (Shear-Locus Kill — every σ ≥ 1, every D ≥ 1).** No scaffold
exists on the composite-jewel shear orbit. Proof structure (codex's
target-weight + leading-pole argument, completed with the independence
step):
1. *Weight decoupling.* w_germ and c_min are ω-homogeneous (weights −1,
   −3), so P_W(w, c) = Σ_{(J,s)} J·a_{J,s} c^s w^{J−1} decomposes by the
   lines (J−1) + 3s = const, and the master identity constrains each
   line separately.
2. *Non-target lines vanish identically.* A non-target line must satisfy
   P_ℓ,W(w_germ, c_min) ≡ 0 as a function; since J(w_germ, c_min) =
   λx⁻² ≠ 0, the pair is functionally independent, so P_ℓ,W = 0 as a
   POLYNOMIAL, i.e. a_{J,s} = 0 for all J ≥ 1 off the target line. P
   reduces to the target line (J−1) + 3s = σ−2 plus p₀(C).
3. *Target count.* μ ≠ 0 forces the s = 0 target coefficient
   a_{σ−1,0} ≠ 0 (the only term reaching x^{2−σ}Y⁰); σ = 1 has no such
   term (μ = 0 immediately), σ = 2 is the m = 0 pole argument.
4. *Leading pole.* On the sheared pair at Y = 0: w_S = −λx⁻¹ + O(x³),
   c_S = γx + O(x⁵). The s = 0 target term contributes
   a_{σ−1,0}(−λ)^{σ−1}x^{−(σ−1)} to f∘S; target terms with s ≥ 1 sit at
   strictly higher order x^{−(σ−1)+4s}; p₀(c_S) is regular. The pole
   cannot cancel — contradicting f∘S ∈ ℂ[x, y]. Hence μ = 0. ∎

The bounded sweeps are now corroboration, not the proof:
- bare-rows full-grid (historical, rev-1.1): μ = 0, σ ≤ 12, D ≤ 8;
- corrected SHEARED-rows, line-restricted (codex checker): μ = 0,
  σ ≤ 8, D ≤ 6;
- corrected SHEARED-rows, FULL grid (this repo,
  verify_composite_sheared.py): μ = 0, σ ≤ 8, D ≤ 6 — and the machine
  independently rediscovers step 2 (all off-line coefficients forced 0).

**Verdict for the sliver.** The composite-jewel stratum E4 is now:
(i) pair-dead on minimal charts (k = 1, 4, 7; both t-signs) — §1;
(ii) pair-ALIVE on the shear locus — §2 — but scaffold-dead there for
     EVERY (σ, D) — §3 Theorem. The shear branch is CLOSED.
(iii) open ONLY for non-shear tail deformations (the b ∈ {2,3} cascade
     freedom; codex's exact window candidate below shows the direction
     is genuinely nonempty at pair level through weight 5), confined by
     the residue clash and CHI-DIV to D ≥ 4 with 3 | J. Closing tool:
     the ξ-ψ-calculus transplant (route 2) — the one remaining closure
     task of the Line Theorem.

**Codex's non-shear window candidate (audit finding E, adopted as
evidence).** An exact tail assignment with λ = γ = 1, free slot
k_{7,3} = 1 and a nonzero b = 2 slot (4,2) = √3/2 (full table in
codexchk_composite_audit.py) satisfies the pair equation exactly through
weight 5, first residual at weight 6, and lies outside even the extended
vertical-shear orbit Y ↦ Y + h(x). Window-partial (not an all-orders
pair) — but it certifies the b ∈ {2,3} freedom is real, exactly as the
23-free-knob cascade picture predicts.

## 4. Honest limits

- The marches are windowed (bucket ≤ 11, weights ≤ 11, b ≤ 18); the
  cascade's regularity (coefficient pattern linear in β) is strong
  structure but not yet a closed form with proof.
- §1 is exact and unconditional for its charts; it is a theorem.
- §3 is now a theorem for all (σ, D) (audit-contributed argument +
  independence step); the sweeps corroborate it machine-side.
- k = 10+ germ family members are fertility-implied but not yet marched;
  the k = 1, 4, 7 spines are identical in mechanism (core (k−2)tγλ) and
  we expect the general-k statement from the same route-2 calculus.

## 5. Witness ledger

- `marcher_composite.py` — clean-room composite marcher: rref-of-[M|I]
  architecture (no symbolic pivoting), resolved-sub imposition chains,
  cancel-normalized conditions, pre-bucket kernel-support recording,
  window asserts. Runs: `out/mc_k1_minimal.json` (terminal b2, core
  −√3γλ), `out/mc_k1_min_e1.json` (same, E1-seeded), `out/mc_k4_minimal.json`
  (terminal b−1, core 2√3γλ), `out/mc_k1_full_deep.json` (buckets → 11,
  forced cascade, no terminal; deep2 rerun with shields_all recording).
- `check_cascade_is_shear.py` — 25 impositions fixpoint-resolved == shear
  at free→0 (value-level identification).
- `verify_composite_jewel.py` — 5/5: A1–A4 shear pair exact; B-sweep
  P-layer μ = 0 forced, all σ ≤ 8, D ≤ 6 (245s, out/verify_composite_jewel.log).
- `probe_composite_cone.py` — fertility/ratio groundwork (Q1–Q3).
- `classify_stratum.py` — E4 recognized as its own stratum; fertile
  k ∈ {1, 4, 7} independent confirmation; 2,100-shape sweep green.
