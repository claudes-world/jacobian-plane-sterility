# The Line Theorem — assembly, v1.3 (2026-07-23 ~14:20 ET, post-delta-audit)

**Status:** the case tree is complete and every stratum is closed except ONE
named sliver (§5: composite-jewel cells), now reduced (COMPOSITE-JEWEL
rev-2, post codex delta-audit) to NON-SHEAR TAIL DEFORMATIONS ONLY: minimal
charts are pair-dead (k = 1, 4, 7, both t-signs); the entire shear locus is
scaffold-dead for EVERY (σ, D) by the Shear-Locus Kill theorem
(target-weight + functional-independence + leading-pole); survivors are
confined to D ≥ 4 with 3 | J by the residue clash and CHI-DIV.
The partition itself is now EXECUTABLE and fuzz-verified (obligation 8
closed: classify_stratum.py, 2,100 shapes + 15-shape zoo, every shape in
exactly one stratum). This document states the theorem once, gives the
partition, and cites the closing result per stratum. Supersedes
LINE-THEOREM-ASSEMBLY-skeleton.md; v1.1 sliver text superseded by §5.

## 0. Statement

**Line Theorem (target).** Every scaffold over Γ = {x = 0} — w ∈ ℂ(x,y) with
poles exactly on Γ, c ∈ ℂ[x,y], P ∈ ℂ[W,C] with f := P(w,c) ∈ ℂ[x,y] and
det J(f,c) = κ ∈ ℂ* — is an affine triangular automorphism.

**Corollary (exclusion form).** No JC2 counterexample arises from triangular
pole hiding over a line: by the Line Dichotomy a scaffold with nonconstant u
would BE one, and none exists. With T2 (hyperbola): the mechanism that
falsifies JC in dimension 3 — Vitushkin's 1999 plane prototype, promoted to
ℂ³ by the July 2026 counterexample — cannot operate over a line or the
standard hyperbola in ℂ². NOT claimed: other pole curves, non-triangular H,
JC2 itself.

## 1. Reductions (all closed)

1. Pair equation J(w,c) = λxᵐ + master identity P_W(w,c) = μx⁻ᵐ. [T2 unit
   machinery] 2. ū = γ ≠ 0; Dichotomy: u const = tame side; else
   counterexample-candidate. [Dichotomy] 3. D = 1, w ∈ ℂ(x), u ∈ ℂ[x],
   u y-affine, fiber ≅ ℂ or ℂ*: dead. [T3 + Abhyankar–Moh–Suzuki,
   Zaidenberg–Lin, Kaliman cited] 4. Poleless w (d = 0) excluded by the
   scaffold definition (poles exactly on Γ ≠ ∅); genuine poleless pairs
   exist (J(y, y²−x) = 1; w = x^a y family) and are tame/T3 territory.
   5. N := deg_y c ≥ 2 after (3); m ∈ ℤ unconstrained (kills below never
   use m ≥ 1; T4's m ≥ 1 + quantization used ONLY where cited).

## 2. The partition (fiber hull of c − t at x → 0)

Lower Newton hull from (0,0) to top vertex (α, N); edges convex (increasing
slopes). Let e₁ = first edge, endpoint (i₁, N₁). Complete trichotomy by
(#edges, structure of e₁, support on edges):

| # | Stratum | Closing result | Status |
|---|---------|----------------|--------|
| A | Single edge, some c-support on the open edge ("dominant") | T5 Discriminant-Ring (evaluate at critical point; GGV Thm 1.5(3) lineage) | CLOSED (scaffold-level, uses T4) |
| B | Single edge, no other support ("strict collapse": ω > 0 off top) | FATEDGE-KILL rev-5 (pair-level, all m ∈ ℤ, J(w,c) ∉ ℂ[x,1/x]) | CLOSED, 5-audited |
| C | ≥2 edges, e₁ bare-ramified (only endpoint supported, N₁ ≥ 2) | Theorem′ Ramified-Edge Kill (ω₁-grading, verbatim proof) | CLOSED (suite E) |
| D | ≥2 edges, N₁ ≥ 2, e₁ multi-supported | T6 divisor-law kill at simple e₁-roots; repeated e₁-roots: first-edge chart → recursed configuration (each recursion lands back in this table one level down; T6's undershoot bounds terminate it) | CLOSED mod §6-item-2 citation archaeology |
| E | ≥2 edges, N₁ = 1 (the split chambers): chart at e₁ → cells (N,k,r), cone edge = second hull edge | subdivided below | — |
| E1 | cone edge binomial (no interior lattice support; automatic when gcd(N−1,r) = 1) | N = 2: T7 (ψ-calculus, jewels, Parity-Pole). N ≥ 3: FERRY Theorem B rev-2.1 (fertility Lemma 1; Uniform Residue −1/(N−1); coset-refined). Which cells exist: T8/T9. | CLOSED, audited |
| E2 | cone edge composite, repeated root | No-Conspiracy telescoping (T6): live 2-edge configurations have simple second-edge roots | CLOSED (+ probe barren) |
| E3 | cone edge composite, distinct roots, generic coefficients | **Composite-Cone Fertility Lemma (NEW, tonight):** the germ recursion gains the termination condition P_k(t) = 0 — an explicit polynomial per k (e.g. (3,6): P₁ ∝ t²−3, P₂ ∝ t(2t²−9)); off its zero locus the cone is germ-barren, hence pair-barren | CLOSED (probe_composite_cone.py; needs prose + audit) |
| E4 | cone edge composite, distinct roots, ON the fertile locus P_k = 0 | **Composite Residue Law (NEW, tonight):** per-root ratio φ(ρ)/φ(0) = d/((N−1)ρE′(ρ)); trace identity holds by sum-of-residues; the clash kills every J with (ratio)^J ≠ 1 | PARTIALLY CLOSED → §5 |
| E0 | ≥3 fiber edges, N₁ = 1, SIMPLE deeper-edge roots (chart hull has middle vertices — not a clean 2-vertex cell) | T6 third-plus-edge undershoot kill (σ_s + I_s undershoots σ₁ by (N_{s−1}−1)(σ_s − σ_{s−1}) > 0 for s ≥ 3) + 2-edge level = E-machinery in the recursed chart | CLOSED mod §6-item-2 citation archaeology [row added 04:25 — was implicit in v1, surfaced by self-audit] |
| F | ≥3 edges, deeper-edge repeated roots (old "conspiracy" shapes) | concave/collapsed: subsumed by B (fat edge) or C. Convex: No-Conspiracy at 2-edge level; D3-TRACE/CHI-DIV for D ≤ 3 / D ≤ N′; recursion of D/E0 edge-by-edge | CLOSED mod §6-item-2 |

**Standing glue facts (used across rows, to state once in the paper):**
(i) c = xu with u(0,0) = γ ≠ 0 forces (1,0) ∈ supp(c): the hull's bottom
vertex is ALWAYS (1,0) — no vertical first edge, no bottom-vertex case
split. (ii) Fractional first-edge slopes: the covering trick x = χ^{q₁}
(handoff P6) transports every ramified-slope configuration to an
integer-slope one with the pole line preserved; the transport lemma
(pairs ↔ pairs, strata ↔ strata) must be stated and cited once — added to
§6 as obligation 7. (iii) Legitimacy/unit constraints only shrink strata:
every kill is a nonexistence statement, so omitting a constraint can only
make our task harder, never create a missed live stratum.

**Completeness:** every c with N ≥ 2 and nonconstant u has a fiber hull; the
hull has one edge (A/B split by open-edge support) or ≥ 2 (C/D/E split by
N₁ and e₁-support; E subdivides by the cone-edge polynomial's factorization:
binomial / repeated / distinct-generic / distinct-fertile). The table's
conditions are exclusive and exhaustive by construction. [Two-paragraph
formal proof to write for the paper; the trichotomy inequalities are
machine-swept per instance.]

## 3. What changed tonight (audit trail)

- Round-1 five-reviewer audit: theorem untouched; FERRY repaired (rev-2.1:
  shift constants, ζ-cosets + (3,2,6) witness, coordinates, m < 0, kernel
  exclusion); STERILITY strengthened (rev-5: J(w,c) ∉ ℂ[x,1/x]).
- Theorem′ (Ramified-Edge Kill) — closes stratum C wholesale.
- Composite-cone strata E2–E4 discovered (the cell machinery had silently
  assumed binomial cones), then 3/4 closed same-session.

## 4. Witness index

verify_artifacts 13/13 · marcher2 controls · verify_residue 22/22 ·
verify_residue_23 15/15 · fertility_residue_sweep 17/17 · trace_hierarchy
9/9 · verify_fatedge_kill 18/18 · verify_326_coset 8/8 ·
staircase_collapse_sweep 905/905 · probe_composite_cone (Q1 barren-generic,
Q2 hunts barren, P_k table) · marcher_composite (k=1/k=4 minimal TERMINAL
cores ∝ t·γ·λ; E1-seeded identical; full-knob deep march → bucket 11, forced
cascade, no terminal) · verify_composite_jewel 5/5 (shear pair exact; P-layer
μ=0 all σ≤8, D≤6) · classify_stratum 2,100 shapes + 15-zoo exactly-one-stratum
· five external reviewer batteries (~1800 systems).

## 5. The one named open sliver: composite-jewel cells [updated 05:10, route 1+3 executed]

At the fertile loci of composite cones (first family: (N,r) = (3,6),
k ≡ 1 mod 3 — machine-confirmed k ∈ {1,4,7} by the classifier's germ
solver — t = ±√3) the germ exists in closed form — at (3,6,1), t = √3:
**w_germ = −λ(x⁻¹ + √3·x²Y)** — and the residue ratios are **primitive
cube roots of unity**: the clash kills all scaffolds with 3 ∤ J.

**Session results 07-23 (COMPOSITE-JEWEL-2026-07-23.md, witnessed):**
1. **Minimal-Chart Kill (theorem):** minimal charts (tails = γx) are
   pair-dead at k = 1 and k = 4 — marcher terminal obstruction at the
   germ's first γ-rung, core ∝ t·γ·λ (−√3γλ at k=1, +2√3γλ at k=4);
   E1-robust, deep-pole-clean. The interior cone coefficient t is itself
   the obstruction (binomial cells shield this rung; composite cannot).
2. **Shear orbit (exact pairs):** (w_germ∘S, c_min∘S) under
   S: Y ↦ Y + γx is a genuine pair family — E4 carries pairs, the
   composite analog of the E1 shear lemma. The full-knob march's forced
   imposition cascade (2 conditions/bucket, absorbed at the b ∈ {0,1}
   lag-2 slots, coefficients −(β−1)√3λ, −β√3λ; no terminal through
   bucket 11) is the shear orbit's shadow + residual b ∈ {2,3} freedom.
3. **P-layer kill on the shear locus (S-conjugation):** scaffolds on
   shear pairs pull back to the bare germ pair; exact sweep kills
   μ ≠ 0 for all σ ≤ 8, D ≤ 6 (σ = 1 by Y-degree counting, σ = 2 by the
   m = 0 pole argument, σ ≥ 3 by the full linear system with
   pole-cancellation rows).

**Remaining open (the true sliver):** non-shear tail deformations of the
composite-jewel pairs ONLY (the shear branch is closed for all (σ, D) —
COMPOSITE-JEWEL rev-2 §3; codex's exact window candidate certifies the
non-shear direction is nonempty at pair level through weight 5). Confined
to D ≥ 4, 3 | J. Closing tool: **ξ-ψ-calculus transplant** (route 2) —
the cone ODE is first-order in ξ; the 𝒟-eigenvalue analog should deliver
the uniform-in-(D, tails, k) statement. This is the one remaining closure
task of the Line Theorem.

## 6. Remaining obligations (complete list)

1. §5 sliver closure (routes above).
2. Citation archaeology for D/F (T6's exact undershoot statement + recursion
   termination — the desktop session holds the derivations; extraction list
   ready for Liam's one-visit dossier).
3. FERRY asks 2–3 genericity wording (GPT-Pro §3.5 discharge — half done).
4. Composite-cone lemmas (E3/E4): prose + external audit round.
5. Completeness proof formal write-out (§2 note).
6. Then: 3a priority note → 3b two-pillar paper (Vitushkin-framed opening).
7. Covering-trick transport lemma: state once, cite handoff P6, verify the
   pair/strata correspondence under x = χ^{q₁} (glue fact ii).
8. ~~Partition fuzzer~~ **CLOSED 07-23 ~04:35** (classify_stratum.py):
   case tree EXECUTABLE as ten independent predicates; 2,100-shape sweep +
   15-shape zoo, every shape in exactly one row, all strata populated
   (out/classify_sweep_full.log). The machine answer to "did we miss a
   stratum" is: no.
## 7. Citation notes from lit pass-3 (2026-07-23 ~04:15)

- Vitushkin 1999 (Math. Notes 66:2, 245–249) VERIFIED: rational ℂ² map,
  constant Jacobian −2, pole on {y=0} — a scaffold over a line in our exact
  sense, killed by its own pole. Opens the paper. ACTION: pull the Springer
  PDF for the explicit formula before printing it (currently known only via
  paraphrase).
- Kraus 1884 VERIFIED (Rodríguez Díaz, CRMath 2025 / arXiv 2512.23614);
  Kraus's flawed step was "controlling ramification at infinity" — quote it.
- Orevkov: cite 1986/87 + 1996 (analytic realization of Vitushkin's
  topological covering). DO NOT cite the "1999 exceptional-set-line ⇒
  singularity" claim — unverifiable in pass-3.
- Gwoździewicz 1993 (alg-geom/9305008) VERIFIED: Keller map injective on
  one line ⇒ automorphism. Foil for the Dichotomy.
- Masuda–Miyanishi 2001 (Math. Ann.) VERIFIED-ADJACENT: étale endos of
  ℂ*-fibered surfaces; cite near the jewels; no scaffold content.
- ⚠ COMPETITOR: Shaska, "Graded Keller maps," arXiv 2607.20210 (22 Jul
  2026): every 𝔾ₘ-equivariant Keller map of ℂ² is an automorphism. Shares
  the grading tool, not the pole/tail machinery. Cite prominently in 3a;
  state the class boundary (equivariant vs our non-graded tails).
