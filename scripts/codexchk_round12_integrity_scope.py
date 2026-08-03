#!/usr/bin/env python3
"""Round-12 ledger-integrity and conditional-scope delta checker."""

from __future__ import annotations

from collections import Counter
import json
import re


def check(name: str, condition: bool) -> None:
    if not condition:
        raise AssertionError(name)
    print(f"PASS  {name}")


ledger = json.load(open("theorem_ledger_v2.json"))
entries = ledger["entries"]
by = {entry["id"]: entry for entry in entries}
paper = open("paper-3b.tex").read()


print("== R3. Ledger structural repairs ==")

check("76 unique ledger IDs", len(entries) == len(by) == 76)
dangling = [
    (entry["id"], dep)
    for entry in entries
    for dep in entry.get("depends_on", [])
    if dep not in by
]
check("zero dangling dependencies", not dangling)
check("12-class closed status_class vocabulary", len(ledger["status_classes"]) == 12)
check(
    "every entry has exactly one declared status_class",
    all(entry.get("status_class") in set(ledger["status_classes"]) for entry in entries),
)

expected_repoints = {
    "CLUSTER-DRT": ["LOCAL-REDUCTION", "WALL", "T4"],
    "LEAFKILL-CLASSIFICATION": [
        "CLUSTER-DRT",
        "UNDERSHOOT-LAW",
        "T4",
        "LOCAL-REDUCTION",
    ],
    "PUISEUX-PRODUCT-FORMULA": ["T4"],
    "T4-UNCOVERED-RESTATEMENT": ["T4", "TRANSPORT-LEMMA"],
}
check(
    "all nine dangling endpoints repointed to the intended live records",
    all(by[entry]["depends_on"] == deps for entry, deps in expected_repoints.items()),
)
semantic_tokens = {
    "LOCAL-REDUCTION": ("Local Reduction", "T1 chart identity"),
    "WALL": ("wall", "C[t]"),
    "T4": ("Exponent/slope/divisor laws", "m>=1"),
    "UNDERSHOOT-LAW": ("undershoot identity", "floor_s"),
    "TRANSPORT-LEMMA": ("transport lemma", "slopes multiply"),
}
check(
    "repoint targets match the intended mathematical titles",
    all(all(token.lower() in by[key]["title"].lower() for token in tokens) for key, tokens in semantic_tokens.items()),
)
check(
    "RECURSION-TERMINATION remains broken pending this round",
    by["RECURSION-TERMINATION"]["status_class"] == "broken",
)
check(
    "rev-3.5 LEAFKILL-CLASSIFICATION is superseded",
    by["LEAFKILL-CLASSIFICATION"]["status_class"] == "superseded",
)
check(
    "historical graph-v2 count corrected to 71",
    "71 math edges" in by["TYPED-THEOREM-GRAPH"]["verification"],
)

spot_ids = [
    "WALL",
    "T4",
    "T11",
    "O1",
    "O2",
    "RECURSION-TERMINATION",
    "LEAFKILL-CLASSIFICATION",
    "LINE-THEOREM-COMPLETE-R7",
    "PUISEUX-PRODUCT-FORMULA",
    "TYPED-THEOREM-GRAPH",
    "LOCAL-REDUCTION",
    "ANCHOREDNESS-CLASSIFICATION",
]
check("status-class spot check covers >=10 records", len(spot_ids) >= 10)
expected_classes = {
    "WALL": "proved",
    "T4": "proved",
    "T11": "window_certified",
    "O1": "superseded",
    "O2": "open",
    "RECURSION-TERMINATION": "broken",
    "LEAFKILL-CLASSIFICATION": "superseded",
    "LINE-THEOREM-COMPLETE-R7": "retracted",
    "PUISEUX-PRODUCT-FORMULA": "proved+witnessed",
    "TYPED-THEOREM-GRAPH": "machine-verified",
}
check(
    "ten unambiguous sampled classes are honest",
    all(by[key]["status_class"] == value for key, value in expected_classes.items()),
)

# Two additional sampled entries are not honest after the R1 counterroute.
check(
    "LOCAL-REDUCTION is still witnessed although it depends on broken termination and claims rows D/F closed",
    by["LOCAL-REDUCTION"]["status_class"] == "witnessed"
    and "RECURSION-TERMINATION" in by["LOCAL-REDUCTION"]["depends_on"]
    and "Rows D/F closed" in by["LOCAL-REDUCTION"]["title"],
)
check(
    "ANCHOREDNESS-CLASSIFICATION still claims universal ascent as proved+witnessed",
    by["ANCHOREDNESS-CLASSIFICATION"]["status_class"] == "proved+witnessed"
    and "refinement ascends to the anchor" in by["ANCHOREDNESS-CLASSIFICATION"]["title"],
)


print("\n== R4. Conditional-scope sweep ==")

remark_start = paper.index(r"\begin{remark}[Conditional scope")
remark_end = paper.index(r"\end{remark}", remark_start)
remark = paper[remark_start:remark_end]
check(
    "scope remark now names ferry/fertility and the composite-consumer interface",
    "uniform-residue mechanism and the binomial fertility classification" in remark
    and "rigidity theorem" in remark
    and "composite cones" in remark,
)
check(
    "scope remark explicitly demotes the N=2 jewel route under UCK",
    "jewel ladder, listed as a pillar in rev-3.8, is not" in remark
    and r"Theorem~\ref{thm:UCK} closes every" in remark,
)
check(
    "Table 1 E1 row names UCK closure and jewel corroboration",
    "E1 & two edges, binomial cone & Unconditional Cell Kill" in paper
    and "jewel/parity arc: independent corroboration" in paper,
)
check(
    "ferry and fertility titles carry [W; R] retagging",
    "Fertility and the shear orbit {[W; R:" in paper
    and "Uniform Residue Theorem {[W; R:" in paper,
)

# Exact contradictory remnants: these are findings because the remark says
# that any load-bearing record lean outside its list is an error.
check(
    "cells prose still says the N=2 arc is closed instead by the [R] jewel ladder",
    "the cell arc is\nclosed instead by the program's jewel ladder [R]" in paper,
)
check(
    "reader guide still says exactly six record-scoped pillars",
    "The [R] items are exactly the six record-scoped" in paper,
)
check(
    "reader guide still says round-11 full-scope audit is pending",
    "the round-11 full-scope audit is pending" in paper,
)

# UCK genuinely subsumes N=2: its statement quantifies over any split-edge
# cell, Lemma K permits arbitrary E (including repeated roots), and its proof
# has no N>=3, residue-ratio, parity, or jewel hypothesis.
uck_start = paper.index(r"\begin{theorem}[Theorem C --- Unconditional Cell Kill")
uck_end = paper.index(r"\end{proof}", uck_start)
uck = paper[uck_start:uck_end]
check(
    "UCK is genuinely N-uniform and does not consume the jewel route",
    "on any split-edge cell" in uck
    and "N \\ge 3" not in uck
    and "jewel" not in uck.lower()
    and "parity" not in uck.lower()
    and r"Lemma~\ref{lem:K}" in uck,
)

print("\nROUND-12 R3/R4 AUDIT COMPLETE")
print("The dependency repoints and ten core status classes are correct.")
print("R3 still has two overclassified live/history records; R4 is directly")
print("contradicted by the unchanged N=2 closure paragraph and stale guide.")
