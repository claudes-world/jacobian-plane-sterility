#!/usr/bin/env python3
"""Round-10b hostile audit of typed theorem graph/checker v2."""

import json
import re
import subprocess


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    print(f"PASS  {name}")


paper = open("paper-3b.tex", encoding="utf-8").read()
graph = json.load(open("theorem_graph_v1.json", encoding="utf-8"))
nodes = {node["id"]: node for node in graph["nodes"]}


print("== A. The shipped checker is green ==")
shipped = subprocess.run(
    ["python3", "check_theorem_graph.py"],
    check=False,
    capture_output=True,
    text=True,
)
check("shipped checker exits zero", shipped.returncode == 0)
check("shipped checker claims G7 clean", "G7 statement citations" in shipped.stdout)


print("\n== B. Exact G7 parser coverage counterexample ==")

# This is the regex copied from check_theorem_graph.py.
shipped_regex = re.compile(
    r"\\begin\{(theorem|lemma|corollary|proposition)\}"
    r"(\[[^\]]*\])?\s*"
    r"\\label\{([^}]+)\}(.*?)\\end\{\1\}",
    re.S,
)
shipped_matches = list(shipped_regex.finditer(paper))

# The title is allowed to contain nested status brackets such as {[W]}.
# Rather than attempting to parse the title, stop at the first label in
# the environment and then capture its statement body.
robust_regex = re.compile(
    r"\\begin\{(theorem|lemma|corollary|proposition)\}"
    r"(?:(?!\\end\{\1\}).)*?"
    r"\\label\{([^}]+)\}(.*?)\\end\{\1\}",
    re.S,
)
robust_matches = list(robust_regex.finditer(paper))

real_graph_ids = set(nodes) - {"thm:B", "thm:T2obj"}
shipped_linted_ids = {match.group(3) for match in shipped_matches}
robust_linted_ids = {match.group(2) for match in robust_matches}
missed = sorted(real_graph_ids - shipped_linted_ids)

check("robust scan finds all 32 TeX-backed graph statements", robust_linted_ids == real_graph_ids)
check("shipped G7 scans only 12 of those 32 statements", len(shipped_linted_ids) == 12)
check("shipped G7 silently misses 20 graph statements", len(missed) == 20)
check("the missed set includes the load-bearing thm:T", "thm:T" in missed)
check("the missed set includes Axis Kill and UCK", {"thm:axiskill", "thm:UCK"} <= set(missed))


print("\n== C. Undeclared live citations hidden by the parser gap ==")
undeclared = {}
for match in robust_matches:
    label = match.group(2)
    if label not in nodes:
        continue
    references = set(re.findall(r"\\ref\{([^}]+)\}", match.group(3)))
    node = nodes[label]
    declared = (
        set(node.get("math", []))
        | set(node.get("frame_deps", []))
        | set(node.get("expository", []))
    )
    bad = sorted(ref for ref in references if ref in nodes and ref not in declared)
    if bad:
        undeclared[label] = bad

expected_hidden = {
    "thm:laws": ["cor:lawscov", "thm:clusterdrt", "thm:undershoot"],
    "thm:T": ["lem:U", "thm:UCK", "thm:axiskill", "thm:leafkill"],
    "thm:rigidity": ["thm:ferry"],
}
check("robust lint finds the exact hidden citation set", undeclared == expected_hidden)
check(
    "thm:T's graph declaration omits every closure cited in its statement",
    set(nodes["thm:T"]["math"]) == {"lem:transport", "thm:T2obj"}
    and set(expected_hidden["thm:T"]).isdisjoint(nodes["thm:T"]["math"]),
)


print("\n== D. Restoring thm:T's statement edges exposes cycles and levels ==")

# Dependency orientation in the JSON/checker is dependency -> consumer.
all_edges = [
    (dependency, node["id"])
    for node in graph["nodes"]
    for dependency in node.get("math", []) + node.get("frame_deps", [])
]
restored_t_edges = [(dependency, "thm:T") for dependency in expected_hidden["thm:T"]]
augmented_edges = all_edges + restored_t_edges

adjacency = {node_id: [] for node_id in nodes}
for source, target in augmented_edges:
    adjacency[source].append(target)


def reaches(source, target):
    stack = [source]
    seen = set()
    while stack:
        current = stack.pop()
        if current == target:
            return True
        if current in seen:
            continue
        seen.add(current)
        stack.extend(adjacency[current])
    return False


# Existing leafkill <- T is T -> leafkill.  Restoring the statement's
# leafkill -> T creates an immediate 2-cycle.
check(
    "restored T<->leafkill statement consumption is cyclic",
    reaches("thm:T", "thm:leafkill") and reaches("thm:leafkill", "thm:T"),
)
# Existing leafkill -> axiskill, plus restored axiskill -> T, yields the
# longer T -> leafkill -> axiskill -> T cycle.
check(
    "restored T/leafkill/Axis closure is cyclic",
    "thm:leafkill" in adjacency["thm:T"]
    and "thm:axiskill" in adjacency["thm:leafkill"]
    and "thm:T" in adjacency["thm:axiskill"],
)

level = {name: index for index, name in enumerate(graph["levels_order"])}
level_violations = [
    dependency
    for dependency in expected_hidden["thm:T"]
    if level[nodes[dependency]["level"]] > level[nodes["thm:T"]["level"]]
]
check(
    "restored thm:T edges expose quasi<-scaffold level violations",
    set(level_violations) == {"thm:UCK", "thm:axiskill", "thm:leafkill"},
)


print("\n== E. Two advertised ID/type claims are not literal ==")
synthetic_ids = sorted(set(nodes) - set(re.findall(r"\\label\{([^}]+)\}", paper)))
check("the graph still contains exactly two non-TeX synthetic ids", synthetic_ids == ["thm:B", "thm:T2obj"])
check(
    "clusterdrt declares transport expository although its statement uses it for slope reduction",
    "lem:transport" in nodes["thm:clusterdrt"].get("expository", [])
    and "lem:transport" not in nodes["thm:clusterdrt"].get("frame_deps", []),
)

print("\nROUND-10b GRAPH COUNTEREXAMPLE CONFIRMED")
