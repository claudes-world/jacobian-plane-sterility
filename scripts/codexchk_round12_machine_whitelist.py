#!/usr/bin/env python3
"""Round-12 structural and behavioral check of the canonical red whitelist."""

from __future__ import annotations

from pathlib import Path
import re
import subprocess
import sys


def check(name: str, condition: bool) -> None:
    if not condition:
        raise AssertionError(name)
    print(f"PASS  {name}")


expected = {
    "codexchk_round9_graph.py",
    "codexchk_round10_laws_puiseux.py",
    "codexchk_round10_graph.py",
    "codexchk_round10b_graph.py",
    "codexchk_round10c_graph.py",
    "peel_rig_standalone.py",
    "leaf_taxonomy_sweep.py",
    "codexchk_round11_graph_scope.py",
    "_cache_peel.py",
}

text = Path("POST-AUDIT-FOLD-CHECKLIST.md").read_text()
section = text.split("## EXPECTED-RED WHITELIST", 1)[1]
rows = re.findall(r"^\|\s*([^|]+?\.py)\s*\|\s*([^|]+)\|\s*([^|]+)\|$", section, re.M)
listed = {artifact.strip() for artifact, _, _ in rows}

check("canonical whitelist contains exactly nine Python artifacts", listed == expected)
check(
    "every whitelist row carries both a pin/reason and red-since disposition",
    all(pin.strip() and since.strip() for _, pin, since in rows),
)
check(
    "round-11 graph snapshot was added with rev-3.8 custody/pre-repair reason",
    any(
        artifact.strip() == "codexchk_round11_graph_scope.py"
        and "rev-3.8 custody digests" in pin
        and "rev-3.9" in since
        for artifact, pin, since in rows
    ),
)
check("all nine whitelist artifacts exist in the repo root", all(Path(name).is_file() for name in expected))

print("\n== behavioral sample: execute the nine whitelisted artifacts ==")
results = {}
for name in sorted(expected):
    proc = subprocess.run(
        [sys.executable, name],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=120,
    )
    results[name] = proc.returncode
    print(f"{name:<44} rc={proc.returncode}")

expected_nonzero = expected - {"leaf_taxonomy_sweep.py"}
check(
    "eight expected-red/retired artifacts exit nonzero",
    all(results[name] != 0 for name in expected_nonzero),
)
check(
    "known-buggy leaf taxonomy archive retains its documented exit 0",
    results["leaf_taxonomy_sweep.py"] == 0,
)

print("\nROUND-12 R5 WHITELIST CHECK COMPLETE")
