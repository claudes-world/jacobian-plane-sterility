# Line Theorem — witness repository

Machine certificates for the Line/Hyperbola sterility program (see
`docs/note-3a.pdf` for the current preprint draft).

Exact arithmetic throughout (sympy; no floats in any proof-bearing step).

**Cite:** Claude (Fable 5, Anthropic) and Liam Corrigan, *Sterility of
triangular pole-hiding over lines and hyperbolas in ℂ²*, Zenodo preprint,
July 23 2026. [doi:10.5281/zenodo.21518214](https://doi.org/10.5281/zenodo.21518214)

## Priority commitment

`frozen/` contains BYTE-EXACT copies of the two core theorem documents whose
SHA-256 hashes were posted publicly on X on 2026-07-23 (~02:26 ET), before
any content reveal:

- `a40f311c2ee8aaabd3dada15b05dca81ef9dda540dcda4bb18c1b98352548dfc  STAIRCASE-STERILITY-2026-07-23.md`
- `58633ec91872c9746f3d6192d38dce0c106f88d94d0637dcf5f5815f3b06dc78  FERRY-uniform-residue-theorem.md`

Verify: `sha256sum frozen/*.md`. Current revisions (post-audit changelogged
successors) live in `docs/`.

## Claim ↔ witness map

| Claim | Script | Expected |
|---|---|---|
| Bracket/unit machinery + basic artifacts | scripts/verify_artifacts.py | 13/13 |
| Uniform Residue ratio −1/(N−1) (fertile binomial cells) | scripts/verify_residue.py, verify_residue_23.py | 22/22, 15/15 |
| Fertility criterion r∣Nk ∧ r∤k across (N,k,r) grid | scripts/fertility_residue_sweep.py | 17/17 |
| Trace identity hierarchy | scripts/trace_hierarchy.py | 9/9 |
| Fat-Edge Kill (strict collapse; J(w,c) ∉ ℂ[x,1/x]) + Theorem′ | scripts/verify_fatedge_kill.py | 18/18 |
| ζ-coset refinement witness (3,2,6) | scripts/verify_326_coset.py | 8/8 |
| Staircase ⇒ collapse (random sweep incl. rational slopes) | scripts/staircase_collapse_sweep.py | 905/905 |
| Composite-cone fertility P_k(t) + ratio groundwork | scripts/probe_composite_cone.py | Q1–Q3 |
| (3,1,3) bucket marcher (clean-room) | scripts/marcher2.py | spine repro |
| Composite-jewel Minimal-Chart Kill k=1,4,7 (+ t=−√3 twin) | scripts/marcher_composite.py | TERMINAL, core (k−2)·t·γ·λ |
| Shear pair exact + P-layer sweeps | scripts/verify_composite_jewel.py, verify_composite_sheared.py, check_cascade_is_shear.py | 5/5 + CLEAN + PASS |
| Partition executable, exactly-one-stratum | scripts/classify_stratum.py | zoo 15/15 + 2,100 green |

Run artifacts (JSON march records, sweep logs): `runs/`.

## Requirements

python3 + sympy (any recent). Every script is standalone:
`python3 scripts/<name>.py`. Marcher CLI examples in each header.

## Audit trail

Six adversarial model-audit rounds are folded into the documents in `docs/`
(changelogs inside). The delta-audit verdict + exact checker for the newest
results: `scripts/codexchk_composite_audit.py` (ALL EXACT CHECKS PASSED).

## License

Scripts: MIT. Documents (`docs/`, `frozen/`): (c) the authors, all rights
reserved pending journal submission; cite doi:10.5281/zenodo.21518214.
