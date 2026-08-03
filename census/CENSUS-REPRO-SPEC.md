# Census reproducibility spec (rev-3.11 witness layer)

The GPT-Pro referee (finding 12, 2026-07-31) flagged that "the 156k
enumeration lacks domain / exclusions / generator-version / expected-count /
digest / mapping." There are in fact **two distinct censuses** behind the
machine boundary. Both are specified here.

---

## Census A — the executable partition census

**Generator:** `classify_stratum.py` (in this repo; the same generator shipped
in the pre-refresh public repo).

**What it certifies.** Every candidate `c ∈ x·ℂ[x,y]` with the `γx` slot and
`N = deg_y c ≥ 2` lands in **exactly one** of the ten strata
`{A, B, C, D, E0, E1, E2, E3, E4, F}`, each carrying its closing-theorem
citation. Per-stratum predicates are evaluated *independently* and the fuzzer
asserts exactly one fires (stronger than a decision tree — catches overlap or
fall-through in the partition itself).

**Domain (deterministic).**
- Invocation: `python3 classify_stratum.py` (default `--seed 20260723`).
- A named zoo of 15 hand structures (flagship cells, jewels, staircases,
  ramified edges, composite cones, …), each with an expected label.
- 2,100 pseudo-random shapes from seven generators with fixed counts:
  `gen_staircase 500`, `gen_ramified 400`, `gen_cell 350`, `gen_cell(comp)
  150`, `gen_fertile_jewel 100`, `gen_multiedge 300`, `gen_random_cloud 300`.
- Germ-ansatz window `GERM_DEG = 16`; fertility search bound `KMAX = 8`.

**Exclusions.** `N < 2`, absent `γx` slot, and non-`x·ℂ[x,y]` inputs are out of
scope (they are not scaffold-legal `c`). The random generators never emit them.

**Expected count / result.** zoo 15/15 PASS; 2,100 shapes ALL GREEN (exactly
one stratum each). Per-stratum tally (seed 20260723):

    A 149   B 561   C 394   D 79   E0 275   E1 365   E2 36   E3 115   E4 100   F 26   (= 2100)

**Digest.** Captured run `classify_stratum_full.out`:

    sha256  0de84e8b1b10fb45c443cc32891db983de36a4562664022eb7c2a1a27b25fa83
    runtime ~68 s (this box, sympy + gmpy2 backend), exit 0

**Honest scope (GPT-Pro finding 11 / triage B4).** The E3/E4 split is a
**bounded diagnostic**, not exact geometric semantics: `KMAX = 8` means an E3
label reads "no fertility among tested candidates through k = 8." A k=10-fertile
cone can be labelled E3 by this surrogate. This does **not** reopen Theorem B
(E3 and E4 are both killed by Thm 7.9 independently of the label), but the
census must be cited as a bounded diagnostic. The rev-4 manuscript separates
the exact geometric partition from this diagnostic; until then the map above is
the honest claim.

---

## Census B — the reduced-slope termination enumeration ("the 156k")

**Provenance of the figure.** paper-3b rev-3.11's lexicographic-termination
proof (`thm:leafkill(1)`) cites *"a separate enumeration [that] checked
**156,354** reduced-slope support configurations through d=25, including gapped
residue-class supports"* — `codex-audit-round13-20260726.md` line 138, results
table line 450. **This enumeration was run auditor-side (codex, round-13) and
was never committed as a standalone script** anywhere in the program repo. That
is precisely the gap the referee named: there is no shippable generator, domain
file, or digest for the exact 156,354 count.

**The mathematical object it verified (round-13 §2.2, exact).** Along a
rational covering refinement, write the edge slope in lowest terms with
denominator `q ≥ 2`. Equal integral x-orders force every supported height into
one residue class mod q, so `e(Z) = Z^{j0}·ẽ(Z^q)`. For a nonzero root ρ the
map `Z ↦ Z^q` is étale, hence

    mult_ρ(e) = mult_{ρ^q}(ẽ) ≤ deg(ẽ) − ord₀(ẽ) ≤ (j1 − j0)/q < d.

"Gapped" support = interior residue-class heights may be absent; deleting an
interior coefficient only lowers `deg(ẽ) − ord₀(ẽ)` and the multiplicity of any
fixed nonzero root. So the child cluster degree strictly drops, supplying the
first-coordinate descent in the lexicographic pair `(d, A−v)`.

**Reconstruction shipped here:** `termination_slope_census.py`.
- Because gapping is *monotone* (it can only lower the bound), the exponential
  family of gapped subsets collapses to a polynomial-size **worst-case
  certificate**: the dense realizer `e(Z) = (Z^q − 1)^m` per triple `(d,q,m)`,
  `m = (j1−j0)/q`.
- Declared domain: `d ∈ [2,25]`, `q ∈ [2,d]`, `m ∈ [1,d−1]`.
- Invocation: `python3 termination_slope_census.py --emit termination_slope_domain.tsv`
- Result: **GREEN**, `RECON_TRIPLE_COUNT = 4900`, max observed nonzero-root
  multiplicity `= 24`, every triple satisfies `mult = m < d`.
- Domain-dump digest:
  `sha256 f87f36b83342a1f08147610c5fec1edd85a6fc5b63162a7392e10566e985168c`
  (`termination_slope_domain.tsv`).

**Reconciliation status — PENDING.** The reconstruction certifies the *same
mathematical bound* the auditor figure attested, over an explicitly declared
domain, and ships GREEN with a digest. It does **not** reproduce the exact
number **156,354**: that count is a function of the auditor's exact
subset-enumeration domain (which subsets, which q- and span-bounds), which is
not recorded in any committed script. The exact 156,354 is therefore marked
**PENDING recovery/reconstruction of the round-13 auditor script**; the
canonical citation remains `codex-audit-round13-20260726.md`. A rev-4 task is to
either recover that script or re-derive the exact count from a committed,
digest-pinned enumerator and update `thm:leafkill(1)`'s parenthetical to cite
the committed generator instead of an auditor-side run.

---

## Summary table

| Census | Generator in repo? | Domain fixed? | Expected count | Digest | Status |
|---|---|---|---|---|---|
| A — partition | yes (`classify_stratum.py`) | yes (seed 20260723) | zoo 15/15 + 2,100 green | `0de84e8b…` | GREEN; bounded-diagnostic scope (B4) |
| B — 156k slope enum | **no** (auditor-side) | recon domain only | 156,354 (auditor) / 4,900 (recon) | recon `f87f36b8…` | mechanism GREEN; exact 156,354 **PENDING** |
