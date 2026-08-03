#!/usr/bin/env python3
"""Leaf-taxonomy quadchotomy sweep (A2 data) — jc-4, 2026-07-25.

For a batch of plane parents, chart every deeper-edge repeated root and
classify the child's LOCAL (cluster) hull:

    SPLIT        edges == 1, span == 1           -> Axis Kill stratum
    SINGLE-EDGE  edges == 1, span >= 2           -> Cluster-DRT dead
    RECURSE      edges >= 2, some deeper repeated root
    MIXED        edges >= 2, all deeper roots simple

The classification is BY CONSTRUCTION exhaustive on (edge count,
first-edge span, deeper-root multiplicity) — that tautology is the
point of the A2 restatement. This sweep provides witnesses per class,
flags interior support (A-analog) and the mu = N1 boundary, and
reports M-hat = m' + sigma + 1 signs for m_root in {1,2}.
"""

import sympy as sp
from fractions import Fraction

x, y, Y, t, z = sp.symbols("x y Y t z")

def support(c):
    pts = {}
    for term in sp.Add.make_args(sp.expand(c)):
        pd = term.as_powers_dict()
        b = int(pd.get(y, 0)); a = int(pd.get(x, 0))
        co = term / (x**a * y**b)
        pts[(a, b)] = sp.expand(pts.get((a, b), 0) + co)
    return {k: v for k, v in pts.items() if v != 0}

def lower_hull(pts_by_height):
    """anchored lower hull: heights ascending, slopes ascending."""
    hs = sorted(pts_by_height)
    cur = (pts_by_height[hs[0]], hs[0])
    hull = [cur]
    rest = [h for h in hs if h > cur[1]]
    while rest:
        best = None
        for h in rest:
            sl = Fraction(pts_by_height[h] - cur[0], h - cur[1])
            if best is None or sl < best[0] or (sl == best[0] and h > best[1]):
                best = (sl, h)
        sl, h = best
        cur = (pts_by_height[h], h)
        hull.append(cur)
        rest = [hh for hh in rest if hh > h]
    return hull  # [(a,b), ...] anchored at height hs[0]

def fiber_hull_edges(c):
    """plane fiber hull of c - t (anchor (0,0) from -t)."""
    pts = support(c)
    byb = {0: 0}  # the -t anchor
    for (a, b), co in pts.items():
        if b not in byb or a < byb[b]:
            byb[b] = min(a, byb.get(b, a))
    hull = lower_hull(byb)
    edges = []
    for i in range(len(hull) - 1):
        (a0, b0), (a1, b1) = hull[i], hull[i+1]
        edges.append((Fraction(a1 - a0, b1 - b0), (a0, b0), (a1, b1)))
    return pts, edges

def edge_poly(pts, edge):
    sl, (a0, b0), (a1, b1) = edge
    E = 0
    for (a, b), co in pts.items():
        if b0 <= b <= b1 and Fraction(a - a0) == sl * (b - b0):
            E += co * z**(b - b0)
    return sp.expand(E)

def chart_child(c, sigma, rho):
    """child of c - t under y = x^{-sigma}(rho + Y); returns support dict."""
    child = sp.expand(c.subs(y, x**(-sigma) * (rho + Y)))
    pts = {}
    for term in sp.Add.make_args(child):
        pd = term.as_powers_dict()
        b = int(pd.get(Y, 0)); a = pd.get(x, 0)
        co = term / (x**a * Y**b)
        key = (sp.Rational(a), b)
        pts[key] = sp.expand(pts.get(key, 0) + co)
    return {k: v for k, v in pts.items() if v != 0}

def classify_child(child_pts, mu):
    """local hull on heights 0..mu (anchor height 0 = the -t row at order 0)."""
    byb = {0: sp.Rational(0)}
    for (a, b), co in child_pts.items():
        if 1 <= b <= mu:
            if b not in byb or a < byb[b]:
                byb[b] = a
    hull = lower_hull(byb)
    edges = []
    for i in range(len(hull) - 1):
        (a0, b0), (a1, b1) = hull[i], hull[i+1]
        edges.append((Fraction(a1 - a0, b1 - b0), (a0, b0), (a1, b1)))
    n_edges = len(edges)
    span1 = edges[0][2][1] - edges[0][1][1] if edges else 0
    interior = False
    if n_edges == 1 and span1 >= 2:
        sl, (a0, b0), (a1, b1) = edges[0]
        for (a, b), co in child_pts.items():
            if b0 < b < b1 and Fraction(a - a0) == sl * (b - b0):
                interior = True
    deeper_repeated = False
    for e in edges[1:]:
        E = 0
        sl, (a0, b0), (a1, b1) = e
        for (a, b), co in child_pts.items():
            if b0 <= b <= b1 and Fraction(a - a0) == sl * (b - b0):
                E += co * z**(b - b0)
        E = sp.expand(E)
        if E != 0 and sp.degree(sp.gcd(sp.Poly(E, z), sp.Poly(sp.diff(E, z), z)), z) > 0:
            deeper_repeated = True
    if n_edges == 1 and span1 == 1:
        cls = "SPLIT"
    elif n_edges == 1:
        cls = "SINGLE-EDGE" + ("(interior=A-analog)" if interior else "(bare)")
    elif deeper_repeated:
        cls = "RECURSE"
    else:
        cls = "MIXED"
    return cls, n_edges, span1, edges

PARENTS = [
    ("codex-witness", x - x**4*y + x**11*y**2 + x**8*y**2 - 2*x**15*y**3 + x**22*y**4),
    ("codex-overlap-E", x - x**4*y + x**8*y**2 - 2*x**15*y**3 + x**22*y**4 + x**30*y**5),
    ("partB-regular", x + x**5*y + x**10*y**2 + x**110*y**3),
    ("partC-tuned", x + x**4*y + 3*x**7*y**2 + 3*x**10*y**3 + x**13*y**4),
    ("deep-split-mix", x + x**3*y + x**8*y**2 + 2*x**13*y**3 + x**18*y**4),
    ("triple-root", x + x**4*y + 3*x**9*y**2 + 3*x**14*y**3 + x**19*y**4 + x**40*y**5),
]

print(f"{'parent':>16} | {'edge':>5} | {'rho':>4} mu | class                    | child edges")
print("-" * 100)
found = {"SPLIT": 0, "SINGLE": 0, "RECURSE": 0, "MIXED": 0}
for name, c in PARENTS:
    pts, edges = fiber_hull_edges(c)
    for ei, e in enumerate(edges):
        sl = e[0]
        if sl.denominator != 1:
            continue  # integer-slope only in this sweep (coverings separate)
        if ei == 0:
            continue  # first edge: generic-simple (anchor edge), never charted
        E = edge_poly(pts, e)
        if E == 0:
            continue
        roots = sp.roots(sp.Poly(E, z))
        for rho, mult in roots.items():
            if mult >= 2 and rho != 0:
                child = chart_child(c, int(sl), rho)
                cls, ne, sp1, cedges = classify_child(child, mult)
                n1 = e[2][1] - e[1][1]
                flag = " [mu=N1-boundary]" if mult == n1 else ""
                key = cls.split("(")[0].split("-")[0] if cls.startswith("SINGLE") else cls
                found["SINGLE" if cls.startswith("SINGLE") else cls] = \
                    found.get("SINGLE" if cls.startswith("SINGLE") else cls, 0) + 1
                print(f"{name:>16} | e{ei+1}:{str(sl):>3} | {str(rho):>4} {mult}  | {cls:<24}{flag} | "
                      + str([(str(x2[0]), x2[1], x2[2]) for x2 in cedges]))
print("\nclass counts:", found)
print("\nNOTE: quadchotomy is exhaustive BY CONSTRUCTION (edge count x first-edge "
      "span x deeper-root simplicity); this sweep exhibits witnesses per class.")
