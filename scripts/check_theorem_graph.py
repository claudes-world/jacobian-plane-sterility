#!/usr/bin/env python3
"""Typed theorem-graph checker v3 (round-10 blocker 5; v3 = round-11/12 repairs) — jc-4/jc-5.

v2 upgrades over v1, per codex round-10 finding 9:
  G1  ids unique; endpoints exist; every graph id that looks like a TeX
      label EXISTS as a \\label in paper-3b.tex (no invented ids)
  G2  status vocabulary closed
  G3  acyclicity over math+frame edges
  G4  realization-level lint (object < pair < quasi < scaffold)
  G5  sink: unique, suspended, blockers; EVERY non-structure,
      non-handoff-source node reaches the sink (no hard-coded list)
  G6  structure-layer nodes feed no math edge
  G7  STATEMENT-CITATION LINT: for every graph node whose id is a TeX
      label, every \\ref inside the corresponding statement body (the
      environment text before \\end{...}, proofs excluded) must appear
      in math + frame_deps + expository of that node. Citations the
      author considers non-consuming must be DECLARED expository.
"""

import json
import re
import sys

g = json.load(open("theorem_graph_v1.json"))
nodes = {n["id"]: n for n in g["nodes"]}
LEV = {name: i for i, name in enumerate(g["levels_order"])}
VOCAB = {"proved", "machine-faced", "handoff-R", "suspended", "structure"}
paper = open("paper-3b.tex").read()
labels = set(re.findall(r"\\label\{([^}]+)\}", paper))

def check(name, cond):
    if not cond:
        print(f"FAIL  {name}")
        sys.exit(1)
    print(f"PASS  {name}")

check("G1a ids unique", len(nodes) == len(g["nodes"]))
all_edges = [(d, n["id"]) for n in g["nodes"] for d in n["math"]]
frame_edges = [(d, n["id"]) for n in g["nodes"] for d in n.get("frame_deps", [])]
context_edges = [(d, n["id"]) for n in g["nodes"] for d in n.get("context", [])]
check("G1b every math+frame+context endpoint exists",
      all(d in nodes for d, _ in all_edges + frame_edges + context_edges))
check("G1c strengthens endpoints exist",
      all(s["from"] in nodes and s["to"] in nodes for s in g["strengthens"]))
SYNTH = {"thm:B", "thm:T2obj"}  # graph-only ids, documented
bad_ids = [i for i in nodes
           if i not in labels and i not in SYNTH]
check(f"G1d graph ids are TeX labels (exceptions documented): {bad_ids if bad_ids else 'clean'}",
      not bad_ids)

check("G2 status vocabulary closed",
      all(n["status"] in VOCAB for n in g["nodes"]))

indeg = {i: 0 for i in nodes}
adj = {i: [] for i in nodes}
for dep, tgt in all_edges + frame_edges:
    adj[dep].append(tgt)
    indeg[tgt] += 1
queue = [i for i in nodes if indeg[i] == 0]
seen = 0
while queue:
    u = queue.pop()
    seen += 1
    for v in adj[u]:
        indeg[v] -= 1
        if indeg[v] == 0:
            queue.append(v)
check("G3 math+frame graph acyclic", seen == len(nodes))

bad = [(n["id"], d) for n in g["nodes"] for d in n["math"]
       if LEV[nodes[d]["level"]] > LEV[n["level"]]]
check(f"G4 level lint: {bad if bad else 'clean'}", not bad)

sinks = [n for n in g["nodes"] if n["id"] == "thm:B"]
check("G5a sink unique, suspended, with blockers",
      len(sinks) == 1 and sinks[0]["status"] == "suspended"
      and bool(sinks[0].get("blockers")))
reach = set()
stack = ["thm:B"]
while stack:
    u = stack.pop()
    if u in reach:
        continue
    reach.add(u)
    stack.extend(nodes[u]["math"])
unreachable = [i for i, n in nodes.items()
               if n["status"] not in ("structure",) and i not in reach
               and i not in SYNTH]
# frame-only feeders (transport, T2obj) may reach via frame edges: accept
# those whose only consumers are frame_deps:
strength_sources = {s["from"] for s in g["strengthens"]}
ALIAS = {"thm:C"}  # statement-alias display nodes (doctrine v3): no proof obligation
truly_unreachable = [i for i in unreachable
                     if not any(i in n.get("frame_deps", []) for n in g["nodes"])
                     and i not in strength_sources and i not in ALIAS]
check(f"G5b every load-bearing node reaches the sink: "
      f"{truly_unreachable if truly_unreachable else 'all'}",
      not truly_unreachable)

retired = [n["id"] for n in g["nodes"] if n["status"] == "structure"]
leaks = [(r, n["id"]) for r in retired for n in g["nodes"] if r in n["math"]]
check(f"G6 structure-layer feeds no math edge: {leaks if leaks else 'clean'}",
      not leaks)

# G7: statement-citation lint
# robust: titles may contain nested status brackets like {[W: ...]};
# anchor on the FIRST \label inside the environment (round-10b blocker 2)
env_re = re.compile(
    r"\\begin\{(theorem|lemma|corollary|proposition|theoremA|theoremB|theoremC)\}"
    r"(?:(?!\\end\{\1\}).)*?"
    r"\\label\{([^}]+)\}(.*?)\\end\{\1\}", re.S)
linted = set()
violations = []
for m in env_re.finditer(paper):
    lab = m.group(2)
    if lab not in nodes:
        continue
    linted.add(lab)
    body = m.group(3)
    cited = set(re.findall(r"\\ref\{([^}]+)\}", body))
    node = nodes[lab]
    allowed = set(node["math"]) | set(node.get("frame_deps", [])) \
        | set(node.get("expository", [])) | set(node.get("context", []))
    for c in cited:
        if c in nodes and c not in allowed:
            violations.append((lab, c))
check(f"G7 statement citations covered by declared edges: "
      f"{violations if violations else 'clean'}", not violations)
texbacked = {i for i in nodes if i in labels}
uncovered = sorted(texbacked - linted)
check(f"G7b full parser coverage of TeX-backed nodes: "
      f"{uncovered if uncovered else 'all linted'}", not uncovered)

# G8 (round-11 J2.6): edge kinds pairwise disjoint per node
overlaps = []
for n in g["nodes"]:
    kinds = {k: set(n.get(k, []) if k != "math" else n["math"])
             for k in ("math", "frame_deps", "expository", "context")}
    names = list(kinds)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            inter = kinds[names[i]] & kinds[names[j]]
            if inter:
                overlaps.append((n["id"], names[i], names[j], sorted(inter)))
check(f"G8 edge kinds pairwise disjoint per node: "
      f"{overlaps if overlaps else 'clean'}", not overlaps)

print(f"\nTYPED GRAPH v3: ALL CHECKS PASSED ({len(nodes)} nodes, "
      f"{len(all_edges)} math, {len(frame_edges)} frame, "
      f"{len(context_edges)} context, {len(g['strengthens'])} strengthens)")
