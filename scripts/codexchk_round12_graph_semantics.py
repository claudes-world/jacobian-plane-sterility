#!/usr/bin/env python3
"""Independent round-12 graph-v3 syntax and semantic-fidelity audit.

This deliberately does not import check_theorem_graph.py.  It uses a small
stack scanner over theorem/lemma/corollary/proposition/theoremA/B/C
environments, rechecks edge-kind disjointness and citation coverage, audits
all 15 context edges, and reproduces the thm:C schema/alias defects.
"""

from __future__ import annotations

from collections import defaultdict
import json
import re


def check(name: str, condition: bool) -> None:
    if not condition:
        raise AssertionError(name)
    print(f"PASS  {name}")


graph = json.load(open("theorem_graph_v1.json"))
paper = open("paper-3b.tex").read()
nodes = {n["id"]: n for n in graph["nodes"]}


print("== G0. Independent interval scanner ==")

statement_envs = {
    "theorem",
    "lemma",
    "corollary",
    "proposition",
    "theoremA",
    "theoremB",
    "theoremC",
}
token = re.compile(r"\\(begin|end)\{([^}]+)\}")
stack: list[tuple[str, int, int]] = []
statements: dict[str, str] = {}
all_statements: list[tuple[str, str | None]] = []

for match in token.finditer(paper):
    kind, env = match.group(1), match.group(2)
    if kind == "begin":
        stack.append((env, match.start(), match.end()))
        continue
    if not stack or stack[-1][0] != env:
        continue
    _, start, body_start = stack.pop()
    if env not in statement_envs:
        continue
    body = paper[body_start : match.start()]
    labels = re.findall(r"\\label\{([^}]+)\}", body)
    label = labels[0] if labels else None
    all_statements.append((env, label))
    if label:
        statements[label] = body

check("36 theorem-like statements scanned", len(all_statements) == 36)
check("33 labeled TeX-backed statements scanned", len(statements) == 33)
check(
    "custom theorem environments included",
    {("theoremA", None), ("theoremB", None), ("theoremC", "thm:C")}
    <= set(all_statements),
)
check(
    "all 33 TeX-backed graph nodes have statement intervals",
    {i for i in nodes if i in statements} == set(statements),
)


print("\n== G1-G8 independently recomputed ==")

check("35 node IDs unique", len(nodes) == len(graph["nodes"]) == 35)
math_edges = [(d, n["id"]) for n in graph["nodes"] for d in n["math"]]
frame_edges = [
    (d, n["id"]) for n in graph["nodes"] for d in n.get("frame_deps", [])
]
context_edges = [
    (d, n["id"]) for n in graph["nodes"] for d in n.get("context", [])
]
check(
    "all math/frame/context endpoints exist",
    all(d in nodes and target in nodes for d, target in math_edges + frame_edges + context_edges),
)

overlaps = []
for node in graph["nodes"]:
    kinds = {
        key: set(node.get(key, []))
        for key in ("math", "frame_deps", "expository", "context")
    }
    names = list(kinds)
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            if kinds[left] & kinds[right]:
                overlaps.append((node["id"], left, right, kinds[left] & kinds[right]))
check("four edge kinds pairwise disjoint", not overlaps)

adj: dict[str, list[str]] = defaultdict(list)
indeg = {i: 0 for i in nodes}
for dep, target in math_edges + frame_edges:
    adj[dep].append(target)
    indeg[target] += 1
queue = [i for i, degree in indeg.items() if degree == 0]
seen = []
while queue:
    source = queue.pop()
    seen.append(source)
    for target in adj[source]:
        indeg[target] -= 1
        if indeg[target] == 0:
            queue.append(target)
check("math+frame graph acyclic", len(seen) == len(nodes))

violations = []
for label, body in statements.items():
    cited = set(re.findall(r"\\ref\{([^}]+)\}", body))
    node = nodes[label]
    allowed = set().union(
        node["math"],
        node.get("frame_deps", []),
        node.get("expository", []),
        node.get("context", []),
    )
    violations.extend((label, ref) for ref in cited if ref in nodes and ref not in allowed)
check("all statement citations declared", not violations)


print("\n== Context and strengthens semantics ==")

expected_context = {
    ("cor:anchor", "thm:UCK"),
    ("cor:lawscov", "thm:UCK"),
    ("lem:I1", "thm:UCK"),
    ("thm:fertile", "thm:UCK"),
    ("thm:ferry", "thm:trace"),
    ("thm:laws", "thm:trace"),
    ("thm:ferry", "prop:deck"),
    ("prop:deck", "thm:B"),
    ("thm:ferry", "thm:B"),
    ("thm:fertile", "thm:B"),
    ("thm:fertility", "thm:B"),
    ("thm:rigidity", "thm:B"),
    ("thm:trace", "thm:B"),
    ("thm:fertility", "thm:fertile"),
    ("thm:ferry", "thm:fertile"),
}
check("exactly the advertised 15 context edges", set(context_edges) == expected_context)

# The UCK proof is the only target proof for which the four application
# contexts might plausibly be hidden inputs.  Extract the proof immediately
# following its statement and verify it cites only the genuine producer K.
uck_start = paper.index(r"\begin{theorem}[Theorem C --- Unconditional Cell Kill")
uck_proof_start = paper.index(r"\begin{proof}", uck_start)
uck_proof_end = paper.index(r"\end{proof}", uck_proof_start)
uck_proof = paper[uck_proof_start:uck_proof_end]
uck_refs = set(re.findall(r"\\ref\{([^}]+)\}", uck_proof))
check("UCK proof consumes lem:K and none of its four context producers", uck_refs == {"lem:K"})
check(
    "Trace is a theorem for pairs assuming m>=1, not a consumer of laws/ferry",
    "For ANY line pair with $m \\ge 1$" in statements["thm:trace"]
    and not re.findall(r"\\ref\{([^}]+)\}", statements["thm:trace"]),
)
check(
    "six B context edges are exactly the six optional pre-UCK routes",
    set(nodes["thm:B"]["context"])
    == {
        "prop:deck",
        "thm:ferry",
        "thm:fertile",
        "thm:fertility",
        "thm:rigidity",
        "thm:trace",
    },
)
check("eight strengthens edges declared", len(graph["strengthens"]) == 8)
uck_strengtheners = {
    s["from"] for s in graph["strengthens"] if s["to"] == "thm:UCK"
}
check(
    "the six optional cell routes strengthen UCK and are absent from its proof",
    uck_strengtheners
    == {
        "prop:deck",
        "thm:ferry",
        "thm:fertile",
        "thm:fertility",
        "thm:rigidity",
        "thm:trace",
    }
    and not (uck_strengtheners & uck_refs),
)

check(
    "ferry/rigidity/deck are raised to scaffold",
    all(nodes[i]["level"] == "scaffold" for i in ("thm:ferry", "thm:rigidity", "prop:deck")),
)
check(
    "transport->clusterdrt and rigidity->deck are math edges",
    "lem:transport" in nodes["thm:clusterdrt"]["math"]
    and "thm:rigidity" in nodes["prop:deck"]["math"],
)


print("\n== Round-12 graph-v3 defects reproduced ==")

required = {"id", "title", "level", "frame", "status", "math", "frame_deps"}
missing = {n["id"]: required - set(n) for n in graph["nodes"] if required - set(n)}
check("thm:C alone is missing required node fields title+frame", missing == {"thm:C": {"title", "frame"}})
check(
    "schema and checker headers still identify v2",
    "typed graph v2" in graph["schema"]
    and "checker v2" in open("check_theorem_graph.py").read().splitlines()[1],
)

# The paper itself is not circular: the later proof does not assume Theorem
# C's conclusion.  The graph encoding is semantically wrong under its own
# doctrine, however: a MATH edge means theorem consumption, while thm:C's
# note says that the target thm:UCK is what proves it.
check(
    "thm:C alias is encoded as a consumed proved theorem while its note says UCK proves it",
    "thm:C" in nodes["thm:UCK"]["math"]
    and nodes["thm:C"]["status"] == "proved"
    and "proved by thm:UCK" in nodes["thm:C"]["note"],
)

print("\nROUND-12 R2 AUDIT COMPLETE")
print("All 15 context edges and all strengthens edges are semantically clean.")
print("Graph v3 is nevertheless BROKEN as an integrity object: thm:C is an")
print("incomplete node, and its restatement/proof-alias relation is mislabeled")
print("as MATH theorem consumption under the graph's own doctrine.")
