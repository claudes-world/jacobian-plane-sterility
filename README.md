# Line Theorem — witness repository

Machine certificates for the Line/Hyperbola sterility program.

Exact arithmetic throughout (sympy; no floats, no Monte Carlo in any
proof-bearing step).

**Manuscript state.** The current manuscript is **paper-3b, rev-3.11**
(supersedes note-3a). Its line-theorem **completeness claim is SUSPENDED**
pending the program's external gates; the machine layer in this repository is
published as an immutable reproducibility record at that revision, not as a
completeness certificate. The manuscript source itself is an internal working
copy and is not distributed here until rev-4 (see *Manuscript source*, below).

**Cite:** Claude (Fable 5, Anthropic) and Liam Corrigan, *Sterility of
triangular pole-hiding over lines and hyperbolas in ℂ²*, Zenodo,
2026. [doi:10.5281/zenodo.21518214](https://doi.org/10.5281/zenodo.21518214)

## Priority commitment

`frozen/` contains BYTE-EXACT copies of the two core theorem documents whose
SHA-256 hashes were posted publicly on X on 2026-07-23 (~02:26 ET), before any
content reveal:

- `a40f311c2ee8aaabd3dada15b05dca81ef9dda540dcda4bb18c1b98352548dfc  STAIRCASE-STERILITY-2026-07-23.md`
- `58633ec91872c9746f3d6192d38dce0c106f88d94d0637dcf5f5815f3b06dc78  FERRY-uniform-residue-theorem.md`

Verify: `sha256sum frozen/*.md`. Post-audit successors of these documents live
in the program's working record.

## Adversarial audit history

Fourteen hostile audit rounds are folded into rev-3.11:

- **Rounds 4–8** — cross-vendor (codex + Gemini + GPT + GPT-Pro), note-3a → rev-3.x arc.
- **Rounds 9–13** — codex full-scope and delta legs; rebuilt the leaf-side
  architecture (Cluster Discriminant-Ring Kill, peel lemma, anchoredness),
  the Puiseux product formula, and the **lexicographic termination proof**
  (round-13: *"no third mathematical termination break exists"*).
- **Round 14** — text-fidelity delta (no mathematics changed).
- **GPT-Pro double-double (2026-07-31)** — two independent paper-3b reviews:
  *3b-review-1* (NO-GO: repaired architecture survives — termination core,
  anchoredness/peel/CDRT, FERRY, fat-edge, Axis Kill, Theorem C, I1 all
  VERIFIED; one exact false quantitative sentence in Thm 6.11; the executable
  partition certifies a bounded diagnostic, not the exact ten-stratum
  semantics; this machine boundary refreshed by the present tag);
  *3b-review-2* (BROKEN: a σ₁ ≥ 1 proof gap in Theorem 3.1). Both are folded
  into the rev-4 repair wave; **completeness remains SUSPENDED.**

Each round's findings are changelogged in the program's working documents. No
scaffold counterexample to the sterility machinery has been found across any
round.

## Machine gates (self-contained, run from the repository root)

| Gate | Command | Expected |
|---|---|---|
| Theorem ledger | `python3 validate_ledger.py` | GREEN — 76 entries, L1–L4 |
| Typed theorem graph v3 | `python3 check_theorem_graph.py` | GREEN — 35 nodes, G5–G8 (needs manuscript source; see below) |
| Executable partition | `python3 classify_stratum.py` | zoo 15/15 + 2,100 shapes, exactly one stratum each |
| Lexicographic termination faces | `python3 ladder_termination_face.py` | GREEN — F1–F5 |
| Termination slope census (recon) | `python3 termination_slope_census.py` | GREEN — span-division bound |

The **expected-RED whitelist** (`whitelist/expected-red-whitelist.json`) lists
round-pinned auditor snapshots that exit nonzero *by design* — they assert a
prior revision's state and go red when it is repaired. They are never findings.
Everything else in the root must be green; any other red is a finding.

## Claim ↔ witness map (rev-3.11)

| Claim | Script | Expected |
|---|---|---|
| Exponent/slope/divisor (T4) laws, ledger sum | scripts/t4_ledger_sum_check.py | GREEN |
| Fat-Edge / Deep-Pole kill (strict collapse, repaired scope) | scripts/deep_pole_fatedge_kill.py, scripts/verify_fatedge_kill.py | GREEN; 18/18 |
| Cluster Discriminant-Ring Kill (anchored blocks) | scripts/cluster_drt_kill.py, scripts/cluster_drt_kill_S4.py | GREEN |
| Peel lemma + anchoredness classification (pass-through nodes) | scripts/cluster_peel_anchoredness.py | GREEN |
| Puiseux product formula → terminal anchoring (+1 constant) | scripts/puiseux_product_formula.py | GREEN |
| Lexicographic termination (valuation climb / degree drop ≤ N₁−1) | scripts/ladder_termination_face.py, scripts/codexchk_round11_repeated_termination.py, scripts/codexchk_round12_termination_multiedge.py | GREEN; F1–F5 + hostile reproducers |
| Reduced-slope span-division bound (round-13 enumeration, reconstruction) | scripts/termination_slope_census.py | GREEN; see census/CENSUS-REPRO-SPEC.md |
| Uniform Residue ratio −1/(N−1) (fertile binomial cells) | scripts/verify_residue.py, scripts/verify_residue_23.py | 22/22, 15/15 |
| Fertility criterion r∣Nk ∧ r∤k | scripts/fertility_residue_sweep.py | 17/17 |
| ζ-coset refinement witness (3,2,6) | scripts/verify_326_coset.py | 8/8 |
| Composite-cone fertility / rigidity groundwork | scripts/probe_composite_cone.py, scripts/verify_composite_jewel.py | Q1–Q3 / 5/5 |
| Trace identity hierarchy | scripts/trace_hierarchy.py | 9/9 |
| Bracket/unit machinery + basic artifacts | scripts/verify_artifacts.py | 13/13 |
| Executable partition, exactly-one-stratum | scripts/classify_stratum.py | zoo 15/15 + 2,100 green |

**Scope note on the partition (GPT-Pro finding 11 / triage B4).** The
executable partition certifies a **bounded diagnostic**: the E3/E4 fertility
split is decided by a finite fertility search (`KMAX = 8`), so an E3 label
means "no fertility among tested candidates through k = 8", not the exact
geometric E3/E4 semantics. It is a strong falsifiable surrogate, not the exact
ten-stratum certificate. The exact geometric partition is separated from this
bounded diagnostic in the rev-4 manuscript.

## Census reproducibility

`census/CENSUS-REPRO-SPEC.md` documents both enumerations behind the machine
boundary — the executable partition census (domain, generator, expected
counts, digest) and the round-13 reduced-slope termination enumeration
(mathematical object, span-division bound, expected-count provenance, and the
reconstruction shipped here). Read it before citing either census.

## Manuscript source

`check_theorem_graph.py` gates the typed theorem graph against the manuscript
(`paper-3b.tex`) in checks G7/G7b. That source is an internal working copy and
is **not** shipped in the rev-3.11-witness tag; run those two checks against
the manuscript distributed at rev-4. The ledger gate, partition, termination
faces, and census are fully self-contained and run without the manuscript.

## Requirements

python3 + sympy (any recent). Every script is standalone:
`python3 scripts/<name>.py`.

## License

Scripts: MIT. Documents (`docs/`, `frozen/`): © the authors, all rights
reserved pending journal submission; cite doi:10.5281/zenodo.21518214.
