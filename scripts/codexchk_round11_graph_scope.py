#!/usr/bin/env python3
"""Round-11 independent graph, external-statement, and scope audit.

This intentionally does not reuse the round-10c environment-stack parser.
The TeX corpus is recovered by a non-nesting interval scan: theorem-class
environments cannot contain another theorem-class environment in this paper,
so each begin token is paired with the first same-kind end token.

The script exits green when the asserted audit observations (including
defects) are reproduced.
"""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent
PAPER = (ROOT / "paper-3b.tex").read_text()
GRAPH = json.loads((ROOT / "theorem_graph_v1.json").read_text())
LEDGER = json.loads((ROOT / "theorem_ledger_v2.json").read_text())
NODES = {node["id"]: node for node in GRAPH["nodes"]}


def check(name: str, condition: bool) -> None:
    if not condition:
        raise AssertionError(name)
    print(f"PASS  {name}")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strip_comments(text: str) -> str:
    """Strip TeX comments without parsing environments."""
    output = []
    for line in text.splitlines(keepends=True):
        cut = len(line)
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
        suffix = "\n" if line.endswith("\n") and cut < len(line) else ""
        output.append(line[:cut] + suffix)
    return "".join(output)


def interval_statements(text: str) -> list[dict[str, object]]:
    """Find non-nested theorem-class intervals without a begin/end stack."""
    kinds = ("theorem", "lemma", "corollary", "proposition")
    begin_re = re.compile(r"\\begin\{(" + "|".join(kinds) + r")\}")
    intervals = []
    for begin in begin_re.finditer(text):
        kind = begin.group(1)
        end = re.search(rf"\\end\{{{re.escape(kind)}\}}", text[begin.end() :])
        if end is None:
            raise AssertionError(f"unclosed {kind} at byte {begin.start()}")
        stop = begin.end() + end.start()
        content = text[begin.end() : stop]
        labels = re.findall(r"\\label\{([^}]+)\}", content)
        intervals.append(
            {
                "kind": kind,
                "start": begin.start(),
                "content": content,
                "labels": labels,
            }
        )
    return intervals


print("== A. Custody ==")
expected = {
    "paper-3b.tex": "703b7a080d94f37e77ba8d3a472b1fc00dbaa22e30d953d95b636c62672f6bc3",
    "paper-3b.pdf": "f520c00c7f8952e968aa77f8dd659d5d92385b5fe15ec373f2d3eb9e25bae58a",
    "theorem_ledger_v2.json": "9dce5bfeb1d680b0ce5591f52f78aff0dd0e31777c9981ebad26a8e87583f9fe",
    "theorem_graph_v1.json": "1aceae24d1e7626ccf40b6dcd8d59fd86bbf8294f994459c60d18eb60fd2ada5",
}
for filename, digest in expected.items():
    check(f"{filename} has the Round-11 digest", sha256(ROOT / filename) == digest)
    check(
        f"{filename} is byte-identical to the Round-11 freeze",
        (ROOT / filename).read_bytes()
        == (ROOT / "frozen-20260725-round11" / filename).read_bytes(),
    )


print("\n== B. Independent statement corpus ==")
paper = strip_comments(PAPER)
intervals = interval_statements(paper)
labeled = {}
unlabeled = []
for interval in intervals:
    labels = interval["labels"]
    if not labels:
        unlabeled.append(interval)
        continue
    check("theorem-class interval has one label", len(labels) == 1)
    label = labels[0]
    check(f"statement label {label} is unique", label not in labeled)
    labeled[label] = interval

tex_labels = set(re.findall(r"\\label\{([^}]+)\}", paper))
tex_backed_graph = set(NODES) & tex_labels
check("independent interval scan finds 32 labeled statements", len(labeled) == 32)
check("the graph has the same 32 TeX-backed statement IDs", set(labeled) == tex_backed_graph)
check(
    "the new conditional-scope remark is outside the theorem corpus",
    "rem:conditional-scope" not in labeled and "rem:conditional-scope" not in NODES,
)
check(
    "the standard theorem-class corpus has one unlabeled Hyperbola Theorem A",
    len(unlabeled) == 1
    and "Theorem A" in unlabeled[0]["content"]
    and "no scaffold exists" in unlabeled[0]["content"],
)

custom_envs = re.findall(r"\\begin\{(theoremA|theoremB|theoremC)\}(.*?)\\end\{\1\}", paper, re.S)
check("intro has three additional custom, unlabeled theorem environments", len(custom_envs) == 3)
check("custom theorem environments have no labels", all(r"\label{" not in body for _, body in custom_envs))

# Independent reproduction of G7/G7b.  Unlike check_theorem_graph.py's single
# tempered regex, this uses the already recovered non-nesting intervals and
# searches only the text following the statement's unique label, exactly the
# doctrine's "statement body" (proofs excluded).
citation_violations = []
linted = set()
for label, interval in labeled.items():
    if label not in NODES:
        continue
    linted.add(label)
    body = interval["content"]
    body_after_label = body[body.index(r"\label{" + label + "}") :]
    cited = set(re.findall(r"\\ref\{([^}]+)\}", body_after_label))
    node = NODES[label]
    allowed = (
        set(node.get("math", []))
        | set(node.get("frame_deps", []))
        | set(node.get("expository", []))
    )
    citation_violations.extend(
        (label, cited_label)
        for cited_label in cited
        if cited_label in NODES and cited_label not in allowed
    )
check("independent statement-citation lint has zero undeclared refs", citation_violations == [])
check("independent statement-citation lint covers all 32 TeX-backed graph nodes", linted == tex_backed_graph)


print("\n== C. Graph-external imports and route-map scope ==")
uck = labeled["thm:UCK"]["content"]
check(
    "labeled thm:UCK explicitly imports the intro statement",
    r"As stated in \S\ref{sec:intro}" in uck,
)
check(
    "the imported intro object is the custom unlabeled theoremC",
    any(kind == "theoremC" and "any split-edge cell" in body for kind, body in custom_envs),
)

route_match = re.search(
    r"\\begin\{remark\}\[Route map --- expository only,"
    r".*?\\end\{remark\}",
    paper,
    re.S,
)
check("exactly one Route-map remark is found", route_match is not None)
route_body = route_match.group(0)
check("Route-map remark is unlabeled", r"\label{" not in route_body)
check(
    "no text outside the Route-map remark names or references it",
    "Route map --- expository only" not in paper[: route_match.start()]
    and "Route map --- expository only" not in paper[route_match.end() :],
)
check(
    "no linted statement consumes Hyperbola Theorem A",
    all("Theorem A" not in interval["content"] for interval in labeled.values()),
)


print("\n== D. Declared graph structure and semantic type defects ==")
math_edges = [
    (dependency, node["id"])
    for node in GRAPH["nodes"]
    for dependency in node.get("math", [])
]
frame_edges = [
    (dependency, node["id"])
    for node in GRAPH["nodes"]
    for dependency in node.get("frame_deps", [])
]
check("declared graph remains 34 nodes / 71 math / 3 frame", (len(NODES), len(math_edges), len(frame_edges)) == (34, 71, 3))

duplicate_types = {}
for node in GRAPH["nodes"]:
    buckets = {
        "math": set(node.get("math", [])),
        "frame": set(node.get("frame_deps", [])),
        "expository": set(node.get("expository", [])),
    }
    overlap = (buckets["math"] & buckets["frame"]) | (buckets["math"] & buckets["expository"]) | (buckets["frame"] & buckets["expository"])
    if overlap:
        duplicate_types[node["id"]] = sorted(overlap)
check(
    "two consumers assign the same dependency incompatible edge types",
    duplicate_types
    == {
        "cor:anchor": ["lem:puiseux"],
        "thm:leafkill": ["thm:T"],
    },
)

check(
    "clusterdrt's slope-transport consumption is dressed as a frame edge",
    "lem:transport" in NODES["thm:clusterdrt"].get("frame_deps", [])
    and "lem:transport" not in NODES["thm:clusterdrt"].get("math", []),
)
check(
    "UCK declares four context edges absent from its displayed proof",
    set(NODES["thm:UCK"]["math"]) - {"lem:K"}
    == {"cor:anchor", "cor:lawscov", "lem:I1", "thm:fertile"},
)
check(
    "rigidity consumes deck in the graph while deck does not consume rigidity",
    "prop:deck" in NODES["thm:rigidity"]["math"]
    and "thm:rigidity" not in NODES["prop:deck"]["math"],
)
check(
    "paper order/wording makes rigidity the producer for the deck locus",
    PAPER.index(r"\label{thm:rigidity}") < PAPER.index(r"\label{prop:deck}")
    and "forcing the conjugate pair" in labeled["thm:rigidity"]["content"]
    and "is exactly the locus" in labeled["prop:deck"]["content"],
)
check(
    "ferry and rigidity are typed pair despite asserting scaffold kills",
    NODES["thm:ferry"]["level"] == "pair"
    and "no scaffold exists" in labeled["thm:ferry"]["content"]
    and NODES["thm:rigidity"]["level"] == "pair"
    and "dies by the residue clash" in labeled["thm:rigidity"]["content"],
)
check(
    "graph doctrine says IDs are TeX labels but checker-only synthetic IDs remain",
    {node_id for node_id in NODES if node_id not in tex_labels} == {"thm:B", "thm:T2obj"}
    and "node IDs ARE TeX labels" in GRAPH["doctrine"],
)
check(
    "Theorem-B sink omits a node for the load-bearing N=2 jewel ladder",
    "jewel" not in " ".join(NODES["thm:B"]["math"]).lower()
    and "$N = 2$" in paper[paper.index(r"\label{rem:conditional-scope}") : paper.index(r"\end{remark}", paper.index(r"\label{rem:conditional-scope}"))],
)


print("\n== E. Ledger integrity ==")
entries = LEDGER["entries"]
ids = {entry["id"] for entry in entries}
check("ledger has 76 unique entries", len(entries) == len(ids) == 76)
missing = sorted(
    (entry["id"], dependency)
    for entry in entries
    for dependency in entry.get("depends_on", [])
    if dependency not in ids
)
check(
    "ledger contains the exact nine dangling dependency endpoints",
    missing
    == [
        ("CLUSTER-DRT", "T-LOCAL-REDUCTION"),
        ("CLUSTER-DRT", "T4-LAWS"),
        ("CLUSTER-DRT", "WALL-LOCAL-FORM"),
        ("LEAFKILL-CLASSIFICATION", "T-LOCAL-REDUCTION"),
        ("LEAFKILL-CLASSIFICATION", "T4-LAWS"),
        ("LEAFKILL-CLASSIFICATION", "UNDERSHOOT-KILL"),
        ("PUISEUX-PRODUCT-FORMULA", "T4-LAWS"),
        ("T4-UNCOVERED-RESTATEMENT", "COVERING-TRANSPORT"),
        ("T4-UNCOVERED-RESTATEMENT", "T4-LAWS"),
    ],
)

declared_statuses = set(LEDGER["statuses"])
used_statuses = Counter(entry["status"] for entry in entries)
check("ledger's advertised status vocabulary does not contain most live statuses", len(set(used_statuses) - declared_statuses) > 20)
entry_by_id = {entry["id"]: entry for entry in entries}
check(
    "refuted rev-3.5 LEAFKILL entry is still marked proved+witnessed",
    entry_by_id["LEAFKILL-CLASSIFICATION"]["status"] == "proved+witnessed"
    and "mixed" in entry_by_id["LEAFKILL-CLASSIFICATION"]["title"].lower(),
)
check(
    "ledger still claims strict-N-descent termination",
    entry_by_id["RECURSION-TERMINATION"]["status"] == "proved"
    and "strict N-descent" in entry_by_id["RECURSION-TERMINATION"]["title"],
)
check(
    "typed-graph ledger entry records the stale 70-edge count",
    "34 nodes / 70 math edges" in entry_by_id["TYPED-THEOREM-GRAPH"]["verification"],
)
check(
    "graph sink blockers retain resolved rev-3.7/R-scope text",
    any("rev-3.7" in blocker for blocker in NODES["thm:B"]["blockers"])
    and any("[R]-handoff scope decision" in blocker for blocker in NODES["thm:B"]["blockers"]),
)


print("\n== F. Conditional-scope completeness ==")
scope_start = paper.index(r"\label{rem:conditional-scope}")
scope_end = paper.index(r"\end{remark}", scope_start)
scope = paper[scope_start:scope_end]
for token in ("prop:reductions", "thm:T5", "lem:I1", "T6/D3-TRACE", "$N = 2$", "Theorem A"):
    check(f"conditional-scope list names {token}", token in scope)
check(
    "FERRY proof explicitly delegates its full derivation outside the paper",
    "fully derived in the program's ferry document" in paper,
)
check(
    "FERRY is not one of the six listed record-scoped pillars",
    "ferry document" not in scope.lower() and "Uniform Residue" not in scope,
)
scope_graph_refs = {
    label for label in re.findall(r"\\ref\{([^}]+)\}", scope) if label in NODES
}
check(
    "conditional-scope remark's graph refs are the three explicit expository interfaces",
    scope_graph_refs == {"prop:reductions", "thm:T5", "lem:I1"},
)
check(
    "no theorem-class statement consumes the conditional-scope remark",
    all("rem:conditional-scope" not in interval["content"] for interval in labeled.values()),
)

print("\nROUND-11 GRAPH/SCOPE DEFECTS REPRODUCED")
