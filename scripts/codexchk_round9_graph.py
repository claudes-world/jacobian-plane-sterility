#!/usr/bin/env python3
"""Audit the literal theorem_ledger_v2.json dependency graph.

The script separates basic graph well-formedness from semantic repair-edge
coverage.  It exits zero when the recorded Round-9 findings reproduce.
"""

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ledger = json.loads((ROOT / "theorem_ledger_v2.json").read_text())
entries = ledger["entries"]
by_id = {entry["id"]: entry for entry in entries}


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    print(f"PASS  {name}")


print("== A. Literal graph shape ==")
check("66 entries are present", len(entries) == 66)
check("all entry ids are unique", len(by_id) == len(entries))
edges = [
    (entry["id"], dependency)
    for entry in entries
    for dependency in entry["depends_on"]
]
check("the ledger contains 100 literal depends_on edges", len(edges) == 100)
check("every dependency names an existing node", all(dep in by_id for _, dep in edges))


def has_cycle():
    state = {}

    def visit(node):
        state[node] = 1
        for dependency in by_id[node]["depends_on"]:
            if state.get(dependency) == 1:
                return True
            if state.get(dependency, 0) == 0 and visit(dependency):
                return True
        state[node] = 2
        return False

    return any(state.get(node, 0) == 0 and visit(node) for node in by_id)


check("the literal graph is acyclic", not has_cycle())


print("\n== B. New repair edges are absent or reversed ==")
axis_deps = set(by_id["AXIS-KILL"]["depends_on"])
local_deps = set(by_id["LOCAL-REDUCTION"]["depends_on"])
check("AXIS-KILL has no DEEP-POLE-KILL dependency", "DEEP-POLE-KILL" not in axis_deps)
check("LOCAL-REDUCTION has no DEEP-POLE-KILL dependency", "DEEP-POLE-KILL" not in local_deps)


def reachable(start, target):
    seen = set()
    stack = [start]
    while stack:
        node = stack.pop()
        if node == target:
            return True
        if node in seen:
            continue
        seen.add(node)
        stack.extend(by_id[node]["depends_on"])
    return False


check(
    "DEEP-POLE-KILL instead depends transitively on AXIS-KILL",
    reachable("DEEP-POLE-KILL", "AXIS-KILL"),
)
check(
    "DEEP-POLE-KILL instead depends transitively on LOCAL-REDUCTION",
    reachable("DEEP-POLE-KILL", "LOCAL-REDUCTION"),
)
deep_dependents = [node for node, dependency in edges if dependency == "DEEP-POLE-KILL"]
check("DEEP-POLE-KILL is a sink with no recorded consumer", deep_dependents == [])
check(
    "adding the paper's AXIS-KILL <- DEEP-POLE-KILL edge would create a cycle",
    reachable("DEEP-POLE-KILL", "AXIS-KILL"),
)


print("\n== C. Status schema and paper-consumer inventory ==")
declared_statuses = set(ledger["statuses"])
actual_statuses = {entry["status"] for entry in entries}
outside = sorted(actual_statuses - declared_statuses)
check("actual status values exceed the declared four-value schema", bool(outside))
print(f"INFO  undeclared status strings: {len(outside)}")

paper_lines = (ROOT / "paper-3b.tex").read_text().splitlines()
consumer_patterns = {
    "T5 quantization": "slope quantization",
    "undershoot value law": "The value law demands",
    "T3 budget transport": "value law reads in the chart",
    "Axis climb m>=1": r"\ord_x W_{\beta^*} \ge m \ge 1",
    "reach ROOT quantization": r"quantization (Theorem~\ref{thm:laws})",
    "all-fat C value law": "the transported value",
    "Puiseux V3 value law": "Taylor expansion along the tracked branch, with the value law",
    "Trace m>=1": r"For ANY line pair with $m \ge 1$",
}
found = {}
for name, needle in consumer_patterns.items():
    found[name] = [
        number
        for number, line in enumerate(paper_lines, 1)
        if needle in line
    ]
    check(f"paper contains direct T4 consumer: {name}", bool(found[name]))
print("INFO  direct semantic T4 consumers:", ", ".join(found))

print("\nALL ROUND-9 GRAPH CHECKS PASSED")
