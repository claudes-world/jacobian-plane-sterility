#!/usr/bin/env python3
"""Round-10c independent theorem-statement and typed-graph audit.

Unlike check_theorem_graph.py, this uses a begin/end environment stack
instead of a title regular expression.  It independently checks current
statement coverage, declared citations, graph cycles, realization levels,
the thm:T split, and the four requested textual confirmations.
"""

from collections import defaultdict, deque
import json
from pathlib import Path
import re


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    print(f"PASS  {name}")


def strip_tex_comments(text):
    cleaned = []
    for line in text.splitlines(keepends=True):
        cut = None
        for index, char in enumerate(line):
            if char != "%":
                continue
            backslashes = 0
            cursor = index - 1
            while cursor >= 0 and line[cursor] == "\\":
                backslashes += 1
                cursor -= 1
            if backslashes % 2 == 0:
                cut = index
                break
        if cut is None:
            cleaned.append(line)
        elif line.endswith("\n"):
            cleaned.append(line[:cut] + "\n")
        else:
            cleaned.append(line[:cut])
    return "".join(cleaned)


def parse_environments(text):
    token_re = re.compile(r"\\(begin|end)\{([^}]+)\}")
    stack = []
    parsed = []
    for match in token_re.finditer(text):
        kind, name = match.group(1), match.group(2)
        if kind == "begin":
            stack.append((name, match.start(), match.end()))
            continue
        if not stack:
            raise AssertionError(f"unmatched end of {name}")
        open_name, open_start, content_start = stack.pop()
        if open_name != name:
            raise AssertionError(
                f"environment mismatch: began {open_name}, ended {name}"
            )
        parsed.append(
            {
                "name": name,
                "start": open_start,
                "end": match.end(),
                "content": text[content_start : match.start()],
            }
        )
    check("TeX environment stack closes", not stack)
    return parsed


paper_raw = Path("paper-3b.tex").read_text()
paper = strip_tex_comments(paper_raw)
graph = json.loads(Path("theorem_graph_v1.json").read_text())
nodes = {node["id"]: node for node in graph["nodes"]}
levels = {name: index for index, name in enumerate(graph["levels_order"])}


print("== A. Stack-based statement coverage ==")
environments = parse_environments(paper)
statement_kinds = {"theorem", "lemma", "corollary", "proposition"}
statement_envs = [
    environment
    for environment in environments
    if environment["name"] in statement_kinds
]
statements = {}
unlabelled = []
for environment in statement_envs:
    labels = re.findall(r"\\label\{([^}]+)\}", environment["content"])
    if not labels:
        unlabelled.append(environment)
        continue
    label = labels[0]
    if label in statements:
        raise AssertionError(f"duplicate parsed statement label {label}")
    environment["label"] = label
    environment["citations"] = set(
        re.findall(
            r"\\(?:auto|eq)?ref\{([^}]+)\}",
            environment["content"],
        )
    )
    statements[label] = environment

all_tex_labels = set(re.findall(r"\\label\{([^}]+)\}", paper))
tex_backed_nodes = {node_id for node_id in nodes if node_id in all_tex_labels}
check(
    "only the intentionally graph-external hyperbola Theorem A is unlabeled",
    len(unlabelled) == 1
    and "Theorem A" in unlabelled[0]["content"]
    and "no scaffold exists" in unlabelled[0]["content"],
)
check("robust scan finds exactly 32 labeled statements", len(statements) == 32)
check("graph has exactly 32 TeX-backed nodes", len(tex_backed_nodes) == 32)
check(
    "robust scan and graph TeX-backed coverage agree exactly",
    set(statements) == tex_backed_nodes,
)


print("\n== B. Independent statement-citation lint ==")
undeclared = []
declared_pairs = []
for label, statement in statements.items():
    node = nodes[label]
    allowed = (
        set(node.get("math", []))
        | set(node.get("frame_deps", []))
        | set(node.get("expository", []))
    )
    for citation in sorted(statement["citations"] & set(nodes)):
        declared_pairs.append((label, citation))
        if citation not in allowed:
            undeclared.append((label, citation))
check("no graph-node statement citation is undeclared", not undeclared)

thm_t_citations = statements["thm:T"]["citations"] & set(nodes)
check(
    "thm:T graph citations are transport plus the leafkill pointer",
    thm_t_citations == {"lem:transport", "thm:leafkill"},
)
check(
    "thm:T consumes transport at quasi level",
    "lem:transport" in nodes["thm:T"]["math"],
)
check(
    "thm:T declares leafkill expository rather than consuming it",
    "thm:leafkill" in nodes["thm:T"].get("expository", [])
    and "thm:leafkill" not in nodes["thm:T"]["math"]
    and "thm:leafkill" not in nodes["thm:T"].get("frame_deps", []),
)


print("\n== C. Declared graph edges ==")
math_edges = [
    (dependency, node["id"])
    for node in graph["nodes"]
    for dependency in node.get("math", [])
]
frame_edges = [
    (dependency, node["id"])
    for node in graph["nodes"]
    for dependency in node.get("frame_deps", [])
]
check("declared edge counts are 71 math and 3 frame", (len(math_edges), len(frame_edges)) == (71, 3))
check(
    "every math/frame endpoint exists",
    all(left in nodes and right in nodes for left, right in math_edges + frame_edges),
)

adjacency = defaultdict(list)
indegree = {node_id: 0 for node_id in nodes}
for dependency, consumer in math_edges + frame_edges:
    adjacency[dependency].append(consumer)
    indegree[consumer] += 1
queue = deque(sorted(node_id for node_id, degree in indegree.items() if degree == 0))
topological_order = []
while queue:
    current = queue.popleft()
    topological_order.append(current)
    for consumer in sorted(adjacency[current]):
        indegree[consumer] -= 1
        if indegree[consumer] == 0:
            queue.append(consumer)
check(
    "math+frame graph is acyclic with all 34 nodes sorted",
    len(topological_order) == len(nodes) == 34,
)

level_violations = [
    (consumer, dependency)
    for dependency, consumer in math_edges
    if levels[nodes[dependency]["level"]] > levels[nodes[consumer]["level"]]
]
check("all math consumptions are realization-level clean", not level_violations)


print("\n== D. thm:T statement/Route-map separation ==")
thm_t_body = statements["thm:T"]["content"]
check(
    "thm:T contains T1-T3 and repeated-root degree descent",
    all(
        fragment in thm_t_body
        for fragment in (
            "(T1)",
            "(T2)",
            "(T3)",
            "each refinement of a repeated root lowers",
            "only finitely many repeated/ramified",
        )
    ),
)
check(
    "thm:T does not contain the old closure roll-call citations",
    not (
        {"lem:U", "thm:UCK", "thm:axiskill"}
        & statements["thm:T"]["citations"]
    ),
)

route_maps = [
    environment
    for environment in environments
    if environment["name"] == "remark"
    and "Route map --- expository only" in environment["content"]
]
check("there is exactly one Route map remark", len(route_maps) == 1)
route_map = route_maps[0]["content"]
check("Route map is unlabeled", r"\label{" not in route_map)
route_map_graph_citations = set(
    re.findall(r"\\(?:auto|eq)?ref\{([^}]+)\}", route_map)
) & set(nodes)
check(
    "Route map's only explicit graph citation is leafkill",
    route_map_graph_citations == {"thm:leafkill"},
)


print("\n== E. Requested textual confirmations ==")
title_start = paper.index(r"\title{")
title_end = paper.index(r"\author{")
title_banner = paper[title_start:title_end]
check(
    "rev-3.7 banner keeps completeness suspended pending round-10c",
    "Manuscript rev-3.7" in title_banner
    and "Completeness claims REMAIN SUSPENDED pending round-10c"
    in title_banner,
)
check(
    "banner chronology restores rev-3.3 then rev-3.4 then rev-3.5",
    title_banner.index("REPAIR STATE (rev-3.3")
    < title_banner.index("mechanical items fixed in rev-3.4")
    < title_banner.index("REV-3.5"),
)

leafkill = statements["thm:leafkill"]["content"]
check(
    "witness wording is narrowed to the two producer scopes",
    "classifier-v2 verdicts on" in leafkill
    and "the remaining classes and" in leafkill
    and "phrase ``per-class witnesses'' overstated" in leafkill,
)
check(
    "thm:leafkill(1) contains the discrete-ladder termination proof",
    all(
        fragment in leafkill
        for fragment in (
            "discrete-ladder argument",
            "strictly lower the",
            "uniformizer is fixed",
            "valuation is an integer",
            "bounded above by",
            "reaches an anchored node in finitely many",
        )
    ),
)

print(
    "\nALL ROUND-10c GRAPH CHECKS PASSED "
    f"({len(statements)} statements, {len(declared_pairs)} declared citation pairs, "
    f"{len(math_edges)} math edges, {len(frame_edges)} frame edges)"
)
