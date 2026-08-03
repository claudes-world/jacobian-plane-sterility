#!/usr/bin/env python3
"""Round-10 adversarial checks of the typed theorem graph."""

import json
import re


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    print(f"PASS  {name}")


paper = open("paper-3b.tex", encoding="utf-8").read()
graph = json.load(open("theorem_graph_v1.json", encoding="utf-8"))
nodes = {node["id"]: node for node in graph["nodes"]}
labels = set(re.findall(r"\\label\{([^}]+)\}", paper))

print("== G1. Identity and coverage checks omitted by the shipped linter ==")
misnamed = {"cor:laws-cov", "lem:puiseux-pf"} & set(nodes)
check("graph uses IDs that are not the corresponding TeX labels",
      misnamed == {"cor:laws-cov", "lem:puiseux-pf"}
      and "cor:lawscov" in labels and "lem:puiseux" in labels)
check("three live split-cell producers are absent from the graph",
      {"thm:fertility", "thm:fertile", "thm:rigidity"} <= (labels - set(nodes)))


print("\n== G2. Direct statement-consumption edges ==")
leaf_start = paper.index(r"\begin{theorem}[Node classification and kills]")
leaf_end = paper.index(r"\end{theorem}", leaf_start)
leaf_statement = paper[leaf_start:leaf_end]
leaf_refs = set(re.findall(r"\\ref\{([^}]+)\}", leaf_statement))
check("leafkill statement directly consumes Axis Kill and terminal anchoring",
      {"thm:axiskill", "cor:anchor"} <= leaf_refs)
check("typed graph omits both direct leafkill edges",
      not ({"thm:axiskill", "cor:anchor"} & set(nodes["thm:leafkill"]["math"])))
check("adding the omitted Axis edge exposes the written mutual dependency",
      "thm:leafkill" in nodes["thm:axiskill"]["math"])

# The cluster proof invokes T2's governance as a mathematical premise, but
# the graph files thm:T as frame-only so its realization-level lint will not
# inspect the pair<-quasi consumption.
check("cluster hides its T2 premise in frame_deps",
      "thm:T" in nodes["lem:cluster"]["frame_deps"]
      and "thm:T" not in nodes["lem:cluster"]["math"])


print("\n== G3. Sink/checker strength ==")
checker = open("check_theorem_graph.py", encoding="utf-8").read()
check("the 'every kill' assertion is only a hard-coded seven-node list",
      'kills = ["thm:fatedge"' in checker
      and '"thm:UCK"' not in checker.split("kills =", 1)[1].split("missing =", 1)[0])
check("graph sink remains explicitly suspended with blockers",
      nodes["thm:B"]["status"] == "suspended"
      and bool(nodes["thm:B"].get("blockers")))

old_checker = open("codexchk_round9_graph.py", encoding="utf-8").read()
check("the claimed pinned 66-entry round-9 snapshot is not actually pinned",
      'json.loads((ROOT / "theorem_ledger_v2.json").read_text())' in old_checker
      and 'len(entries) == 66' in old_checker)

print("\nROUND-10 GRAPH CHECKS PASSED (THE CHECKER PASSES A NON-FAITHFUL GRAPH)")
