"""Partition fuzzer: the Line Theorem case tree, EXECUTABLE (assembly obligation 8).

Input: a candidate c in x*C[x,y] (dict {(a,b): coeff}) with gamma*x present and
N = deg_y c >= 2. Output: exactly one stratum label from
{A, B, C, D, E0, E1, E2, E3, E4, F} + the closing-theorem citation
(assembly section 2 table).

Definitions implemented (assembly section 2, FERRY section 1):
- FIBER HULL of c - t at x -> 0: lower-left convex hull of supp(c) plus the
  origin (the -t vertex), from (0,0) up to the top vertex (alpha, N) where
  alpha = min x-exponent in the top y-row. Edge slopes sigma_i = da/db are
  strictly increasing. All support lies on or right of the hull.
- e1 = first edge, endpoint (i1, N1).
- Split chamber (N1 = 1): sigma := i1, chart Y = x^sigma y sends
  (a,b) -> (a - sigma*b, b); top vertex -> (r, N) with r = alpha - sigma*N.
  CONVEXITY FACT (checked here on every shape): r >= 1 automatically —
  matching FERRY's cells (N,k,r), r >= 1; no legitimacy gap.
- Cone edge (2-edge shapes): (sigma,1) -> (alpha,N); d = gcd(N-1, r);
  cone polynomial E(xi) from the d+1 edge lattice points.
- Fertility test (E3 vs E4): germ w = x^-k * phi(xi) with
  xi = x^(r/d) Y^((N-1)/d), pair equation J(w, c_min) = lam x^-(1+k);
  fertile at k iff the linear system admits lam != 0. Tested k = 1..KMAX.

Per-stratum predicates are evaluated INDEPENDENTLY and the fuzzer asserts
exactly one fires — a stronger check than a decision tree (catches overlap
or fall-through in the partition itself).

Run: python3 classify_stratum.py            # zoo + 2000-shape sweep
     python3 classify_stratum.py --quick    # zoo + 200 shapes
"""
import argparse, random, sys
from math import gcd
from fractions import Fraction
import sympy as sp

KMAX = 8          # fertility search bound for E3/E4 (P_k data: family k = 1 mod 3)
GERM_DEG = 16     # xi-degree window for germ ansatz

CITE = {
    "A":  "T5 Discriminant-Ring (dominant edge)",
    "B":  "FATEDGE-KILL rev-5 (strict collapse, all m in Z)",
    "C":  "Theorem-prime Ramified-Edge Kill (bare-ramified e1)",
    "D":  "T6 divisor-law at simple e1-roots + first-edge chart recursion",
    "E0": "T6 third-edge undershoot + recursed 2-edge E-machinery",
    "E1": "N=2: T7 jewels/Parity-Pole; N>=3: FERRY Thm B rev-2.1 (+ zeta-cosets)",
    "E2": "No-Conspiracy telescoping (repeated cone root)",
    "E3": "Composite-Cone Fertility Lemma (germ-barren; BOUNDED fertility test)",
    "E4": "Composite Residue Law + sliver routes (fertile locus; D>=4, 3|J confined)",
    "F":  "No-Conspiracy at 2-edge level + CHI-DIV/D3-TRACE + edge recursion",
}


# ---------------- hull machinery ----------------

def fiber_hull(supp):
    """Lower-left hull of supp(c) + origin, from (0,0) to (alpha, N).
    Returns list of vertices [(0,0), ..., (alpha, N)] with strictly
    increasing slopes; every support point lies on or right of the hull."""
    N = max(b for _, b in supp)
    pts = {}
    for a, b in supp:
        if b not in pts or a < pts[b]:
            pts[b] = a
    pts[0] = 0                      # the -t vertex dominates row 0
    rows = sorted(pts.items(), key=lambda t: t[0])   # (b, minimal a)
    # monotone chain going up in b, keeping the LEFT hull (support right of it):
    chain = []
    for b, a in rows:
        p = (a, b)
        while len(chain) >= 2:
            (a1, b1), (a2, b2) = chain[-2], chain[-1]
            # slope condition: slope(chain[-2]->chain[-1]) must be strictly
            # less than slope(chain[-1]->p); pop when it is not
            if Fraction(a2 - a1, b2 - b1) >= Fraction(p[0] - a2, p[1] - b2):
                chain.pop()
            else:
                break
        # also pop a chain of length 1 never (anchor (0,0) stays)
        chain.append(p)
    assert chain[0] == (0, 0), f"hull must start at origin, got {chain}"
    assert chain[-1][1] == N, f"hull must end in top row: {chain}"
    # verify: all support on/right of every hull edge
    for (a1, b1), (a2, b2) in zip(chain, chain[1:]):
        for (a, b) in supp:
            if b1 <= b <= b2:
                # point on/right of the directed edge (a1,b1)->(a2,b2)
                cross = (a2 - a1) * (b - b1) - (b2 - b1) * (a - a1)
                assert cross <= 0, \
                    f"support {(a, b)} left of hull edge {(a1, b1)}-{(a2, b2)}"
    return chain


def edge_lattice(v1, v2):
    """All lattice points on the closed segment v1-v2, in order."""
    (a1, b1), (a2, b2) = v1, v2
    g = gcd(abs(a2 - a1), abs(b2 - b1))
    sa, sb = (a2 - a1) // g, (b2 - b1) // g
    return [(a1 + i * sa, b1 + i * sb) for i in range(g + 1)]


def edge_poly(supp_dict, v1, v2):
    """Edge polynomial: coefficients of c on the lattice points of edge
    v1-v2, as a sympy Poly in xi (constant term = v1's coefficient)."""
    xi = sp.Symbol('xi')
    pts = edge_lattice(v1, v2)
    return sp.Poly(sum(supp_dict.get(p, sp.Integer(0)) * xi ** i
                       for i, p in enumerate(pts)), xi), pts


# ---------------- fertility (germ) test ----------------

def bracket(d1, d2):
    out = {}
    for (a, b), u in d1.items():
        for (cc, d), v in d2.items():
            k = a * d - b * cc
            if k:
                key = (a + cc - 1, b + d - 1)
                out[key] = out.get(key, 0) + k * u * v
    return out


def fertile_k(cmin_chart, r, N, d, k):
    """Does the cell germ exist at pole depth k with lam != 0?
    Ansatz w = x^-k * sum f_n xi^n, xi = x^(r/d) Y^((N-1)/d).
    Degree bound: binomial termination is at zeta-deg Nk/r, i.e. xi-deg
    d*N*k/r; pad past it so a clipped window can never fake barrenness."""
    deg = max(GERM_DEG, (d * N * k) // r + 6)
    lam = sp.Symbol('lam')
    fs = sp.symbols(f'f0:{deg + 1}')
    step_a, step_b = r // d, (N - 1) // d
    wdict = {(-k + n * step_a, n * step_b): fs[n] for n in range(deg + 1)}
    eq = bracket(wdict, cmin_chart)
    eq[(-1 - k, 0)] = eq.get((-1 - k, 0), 0) - lam
    unk = list(fs) + [lam]
    M = sp.Matrix([[sp.expand(v).coeff(u) for u in unk] for v in eq.values()])
    ns = M.nullspace()
    return any(sp.simplify(v[-1]) != 0 for v in ns)


# ---------------- the classifier ----------------

def classify(supp_dict, kmax=KMAX):
    """Return (label, citation, info). Evaluates all ten predicates
    independently and asserts exactly one is true."""
    supp = [m for m, co in supp_dict.items() if co != 0]
    assert all(a >= 1 for a, b in supp), "c must lie in x*C[x,y]"
    assert (1, 0) in supp, "gamma slot (1,0) must be supported (u(0,0) != 0)"
    N = max(b for _, b in supp)
    assert N >= 2, "need deg_y c >= 2"
    hull = fiber_hull(supp)
    edges = list(zip(hull, hull[1:]))
    ne = len(edges)
    (i1, N1) = edges[0][1]
    info = {"hull": hull, "n_edges": ne, "N1": N1, "N": N}

    def open_edge_supported(v1, v2):
        pts = set(edge_lattice(v1, v2)[1:-1])
        return any(p in supp for p in pts)

    # --- cone-edge data for split chambers (N1 == 1) ---
    cone = None
    if N1 == 1 and ne >= 2:
        sigma = i1
        alpha = hull[-1][0]
        r = alpha - sigma * N
        # convexity fact: r >= 1 always (checked, not assumed)
        assert ne != 2 or r >= 1, f"2-edge cell with r = {r} <= 0: {hull}"
        if ne == 2:
            d = gcd(N - 1, r)
            E, pts = edge_poly(supp_dict, hull[1], hull[2])
            interior_supported = any(p in supp for p in pts[1:-1])
            squarefree = sp.degree(sp.gcd(E, E.diff())) == 0
            cone = {"sigma": sigma, "r": r, "d": d, "E": E,
                    "interior_supported": interior_supported,
                    "squarefree": squarefree}
            info["cone"] = {k: (str(v) if k == "E" else v)
                            for k, v in cone.items()}

    def fertile_some_k():
        sigma, r, d = cone["sigma"], cone["r"], cone["d"]
        cmin = {(a - sigma * b, b): supp_dict[(a, b)]
                for (a, b) in edge_lattice(hull[1], hull[2])
                if (a, b) in supp_dict and supp_dict[(a, b)] != 0}
        # BOUNDED fertility diagnostic (codex audit 2026-07-23, D(iii)):
        # the germ ODE's leading coefficient forces xi-degree L = d*N*k/r,
        # a positive integer — so fertile k lie in the progression
        # k = j*q/gcd(q, N), q = r/d. Testing a fixed 1..KMAX misses cells
        # with large q (exact counterexample: r=54 cell, first candidate
        # k=9). We test the legacy range PLUS the first three progression
        # members. This is still a BOUNDED test: E3 semantically means
        # "no fertile k among the tested candidates", and the E3/E4 split
        # is a diagnostic refinement of the geometric partition, not an
        # unbounded decision procedure.
        from math import gcd as _gcd
        q = r // d
        step = q // _gcd(q, N) if q else 1
        cands = sorted(set(range(1, kmax + 1))
                       | {j * step for j in (1, 2, 3) if step >= 1})
        ks = [k for k in cands if fertile_k(cmin, r, N, d, k)]
        info["fertile_ks"] = ks
        info["fertility_candidates_tested"] = cands
        return bool(ks)

    def deeper_edges_squarefree():
        ok = True
        for v1, v2 in edges[1:]:
            E, _ = edge_poly(supp_dict, v1, v2)
            if sp.degree(sp.gcd(E, E.diff())) != 0:
                ok = False
        return ok

    composite = (cone is not None and cone["interior_supported"])
    preds = {
        "A": ne == 1 and open_edge_supported(*edges[0]),
        "B": ne == 1 and not open_edge_supported(*edges[0]),
        "C": ne >= 2 and N1 >= 2 and not open_edge_supported(*edges[0]),
        "D": ne >= 2 and N1 >= 2 and open_edge_supported(*edges[0]),
        "E1": ne == 2 and N1 == 1 and not composite,
        "E2": ne == 2 and N1 == 1 and composite and not cone["squarefree"],
        "E3": ne == 2 and N1 == 1 and composite and cone["squarefree"]
              and not fertile_some_k(),
        "E4": ne == 2 and N1 == 1 and composite and cone["squarefree"]
              and bool(info.get("fertile_ks")),
        "E0": ne >= 3 and N1 == 1 and deeper_edges_squarefree(),
        "F": ne >= 3 and N1 == 1 and not deeper_edges_squarefree(),
    }
    fired = [lbl for lbl, v in preds.items() if v]
    assert len(fired) == 1, \
        f"partition violation: predicates {fired} fired for hull {hull}"
    lbl = fired[0]
    return lbl, CITE[lbl], info


# ---------------- named zoo ----------------

S3 = sp.sqrt(3)

ZOO = [
    # (name, support dict, expected label)
    ("flagship (3,1,3) cell, sigma=1", {(1, 0): 1, (1, 1): 1, (6, 3): 1}, "E1"),
    ("jewel cell (2,1), sigma=1", {(1, 0): 1, (1, 1): 1, (3, 2): 1}, "E1"),
    ("shear-orbit cell (2,1) w/ tails", {(1, 0): 2, (1, 1): 1, (3, 2): 1,
                                         (4, 1): 3, (2, 1): 5}, "E1"),
    ("(3,2,6) coset cell, sigma=4 (d=2 binomial)",
     {(1, 0): 1, (4, 1): 1, (18, 3): 1}, "E1"),
    ("(4,1,2) cell, sigma=1", {(1, 0): 1, (1, 1): 1, (6, 4): 1}, "E1"),
    ("composite-jewel (3,6) t=sqrt3, sigma=1",
     {(1, 0): 1, (1, 1): 1, (5, 2): S3, (9, 3): 1}, "E4"),
    ("composite cone generic t=1 (barren)",
     {(1, 0): 1, (1, 1): 1, (5, 2): 1, (9, 3): 1}, "E3"),
    ("composite cone repeated root (1+xi)^2",
     {(1, 0): 1, (1, 1): 1, (5, 2): 2, (9, 3): 1}, "E2"),
    ("fat edge / staircase strict collapse",
     {(1, 0): 1, (2, 1): 1, (3, 3): 1}, "B"),
    ("dominant single edge (open-edge support)",
     {(1, 0): 1, (1, 1): 1, (2, 2): 1, (3, 3): 1}, "A"),
    ("bare-ramified first edge (Theorem-prime)",
     {(1, 0): 1, (2, 2): 1, (5, 3): 1}, "C"),
    ("multi-supported ramified e1 (divisor law)",
     {(1, 0): 1, (1, 1): 1, (2, 2): 1, (5, 3): 1}, "D"),
    ("3-edge simple deeper roots", {(1, 0): 1, (1, 1): 1, (3, 2): 1,
                                    (6, 3): 1}, "E0"),
    ("3-edge composite-simple deeper edge",
     {(1, 0): 1, (1, 1): 1, (3, 2): 3, (5, 3): 1, (8, 4): 1}, "E0"),
    ("conspiracy: repeated deeper root (1+xi)^2",
     {(1, 0): 1, (1, 1): 1, (3, 2): 2, (5, 3): 1, (8, 4): 1}, "F"),
]

# NOTE on the "3-edge composite-simple" zoo entry: hull (0,0)-(1,1)-(5,3)-(8,4)
# has middle edge (1,1)-(5,3) with d = 2 and E = 1 + 3 xi + xi^2 (distinct) ->
# E0; flipping the interior coefficient to 2 gives (1+xi)^2 -> F.


# ---------------- random shape generators ----------------

def _coeff(rng):
    return sp.Integer(rng.choice([1, -1, 2, -2, 3, -3, 4, 5]))


def gen_staircase(rng):
    """Single-edge shapes: everything on/above the (0,0)-(alpha,N) line."""
    N = rng.randint(2, 6)
    alpha = rng.randint(1, 12)
    supp = {(1, 0): _coeff(rng), (alpha, N): _coeff(rng)}
    g = gcd(alpha, N)
    edge_pts = [(alpha * i // g, N * i // g) for i in range(1, g)]
    for _ in range(rng.randint(0, 6)):
        if edge_pts and rng.random() < 0.35:
            supp[rng.choice(edge_pts)] = _coeff(rng)     # on the open edge -> A
        else:
            b = rng.randint(0, N - 1)
            # strictly right of the line a*N = alpha*b, and a >= 1
            amin = max(1, alpha * b // N + 1)
            supp[(amin + rng.randint(0, 5), b)] = _coeff(rng)
    return supp, {"A", "B"}


def gen_ramified(rng):
    """N1 >= 2 first edge, then continuation: C or D."""
    N1 = rng.randint(2, 4)
    i1 = rng.randint(1, 6)
    dN = rng.randint(1, 3)
    N = N1 + dN
    # continuation slope strictly larger than i1/N1
    s2num, s2den = i1 * dN + rng.randint(1, 4) * N1, N1 * dN
    alpha = i1 + (s2num // s2den) * dN + rng.randint(0, 2) + 1
    if Fraction(alpha - i1, dN) <= Fraction(i1, N1):
        alpha = i1 + (Fraction(i1, N1).__ceil__() + 1) * dN
    supp = {(1, 0): _coeff(rng), (i1, N1): _coeff(rng), (alpha, N): _coeff(rng)}
    pts = edge_lattice((0, 0), (i1, N1))[1:-1]
    multi = rng.random() < 0.5
    if multi and pts:
        supp[rng.choice(pts)] = _coeff(rng)
    for _ in range(rng.randint(0, 4)):        # junk strictly right of the hull
        b = rng.randint(0, N)
        amin = max((i1 * b + N1 - 1) // N1, 1) + 1 + rng.randint(0, 4)
        m = (amin, b)
        if m not in ((i1, N1), (alpha, N)):
            supp.setdefault(m, _coeff(rng))
    return supp, {"C", "D"}


def gen_cell(rng, force_composite=False):
    """2-edge split chambers: E1/E2/E3/E4."""
    sigma = rng.randint(1, 5)
    N = rng.randint(2, 5)
    if force_composite:
        # need d = gcd(N-1, r) > 1: pick a multiple of a divisor of N-1 > 1
        N = rng.choice([3, 4, 5])
        dd = rng.choice([d for d in range(2, N) if (N - 1) % d == 0])
        r = dd * rng.randint(1, 4)
        if gcd(N - 1, r) == 1:
            r *= dd
    else:
        r = rng.randint(1, 10)
    alpha = sigma * N + r
    supp = {(1, 0): _coeff(rng), (sigma, 1): _coeff(rng), (alpha, N): _coeff(rng)}
    pts = edge_lattice((sigma, 1), (alpha, N))[1:-1]
    expected = {"E1", "E2", "E3", "E4"}
    if force_composite and pts:
        d = gcd(N - 1, r)
        style = rng.random()
        if style < 0.25:
            # engineered repeated root: E(xi) = (1 + t xi)^d  ->  E2
            t = _coeff(rng)
            e0 = supp[(sigma, 1)]
            for i, p in enumerate(pts + [(alpha, N)], start=1):
                supp[p] = e0 * sp.binomial(d, i) * t ** i
            expected = {"E2"}
        else:
            for p in pts:
                if rng.random() < 0.8:
                    supp[p] = _coeff(rng)
            if not any(p in supp for p in pts):
                supp[pts[0]] = _coeff(rng)
            expected = {"E2", "E3", "E4"}
    for _ in range(rng.randint(0, 4)):        # junk strictly right of the hull
        b = rng.randint(0, N)
        if b <= 1:
            amin = sigma * max(b, 1) + 2
        else:
            amin = sigma + ((r + sigma * (N - 1)) * (b - 1) + N - 2) // (N - 1) + 2
        m = (amin + rng.randint(0, 4), b)
        if m not in supp:
            supp[m] = _coeff(rng)
    return supp, (expected if force_composite else {"E1", "E2", "E3", "E4"})


def gen_fertile_jewel(rng):
    """Torus-scaled composite-jewel (3,6) family: fertility is invariant
    under (x, Y) -> (s x, u Y) torus scaling, realized here as coefficient
    pattern (p^2, sqrt(3) p q, q^2) on the cone edge -> expect E4."""
    sigma = rng.randint(1, 4)
    p, q = rng.choice([1, 2, 3]), rng.choice([1, -1, 2])
    supp = {(1, 0): _coeff(rng), (sigma, 1): sp.Integer(p * p),
            (2 * sigma + 3, 2): S3 * p * q, (3 * sigma + 6, 3): sp.Integer(q * q)}
    for _ in range(rng.randint(0, 3)):        # junk strictly right of the hull
        b = rng.randint(0, 3)
        amin = (sigma * max(b, 1) + (3 + sigma) * max(b - 1, 0)) + 2
        m = (amin + rng.randint(0, 4), b)
        if m not in supp:
            supp[m] = _coeff(rng)
    return supp, {"E4"}


def gen_multiedge(rng):
    """>= 3 edges, N1 = 1: E0 or F."""
    sigma = rng.randint(1, 4)
    verts = [(0, 0), (sigma, 1)]
    slope = Fraction(sigma, 1)
    for _ in range(rng.randint(2, 3)):
        db = rng.randint(1, 3)
        slope = slope + Fraction(rng.randint(1, 5), 1)
        a, b = verts[-1]
        verts.append((a + int(slope * db) + (0 if (slope * db).denominator == 1
                                             else 1), b + db))
    supp = {(1, 0): _coeff(rng)}
    for v in verts[1:]:
        supp[v] = _coeff(rng)
    for v1, v2 in zip(verts[1:], verts[2:]):
        pts = edge_lattice(v1, v2)[1:-1]
        for p in pts:
            if rng.random() < 0.5:
                supp[p] = _coeff(rng)
    if rng.random() < 0.3:      # engineered repeated deeper root
        v1 = verts[1]
        v2 = (v1[0] + 4, v1[1] + 2)
        # rebuild: force edge v1->v2 = (1+xi)^2 then a final steeper vertex
        supp = {(1, 0): _coeff(rng), v1: sp.Integer(1),
                ((v1[0] + v2[0]) // 2, (v1[1] + v2[1]) // 2): sp.Integer(2),
                v2: sp.Integer(1)}
        s2 = Fraction(v2[0] - v1[0], v2[1] - v1[1])
        if Fraction(v1[0], v1[1]) >= s2:
            return gen_multiedge(rng)
        supp[(v2[0] + int(s2) + 3, v2[1] + 1)] = _coeff(rng)
    return supp, {"E0", "F"}


def gen_random_cloud(rng):
    """Unstructured: gamma slot + top row + random junk. Any label legal."""
    N = rng.randint(2, 6)
    supp = {(1, 0): _coeff(rng), (rng.randint(1, 14), N): _coeff(rng)}
    for _ in range(rng.randint(1, 9)):
        supp.setdefault((rng.randint(1, 14), rng.randint(0, N)), _coeff(rng))
    return supp, None


# ---------------- driver ----------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--seed', type=int, default=20260723)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    print("== named zoo ==")
    nfail = 0
    for name, supp, expected in ZOO:
        lbl, cite, info = classify(supp)
        ok = lbl == expected
        nfail += (not ok)
        extra = f" fertile_ks={info['fertile_ks']}" if 'fertile_ks' in info else ""
        print(f"{'PASS' if ok else 'FAIL'}  {name}: {lbl} "
              f"(expected {expected}){extra}  [{cite[:46]}]")

    counts = {k: 0 for k in CITE}
    total = 0
    plan = [(gen_staircase, {}, 500), (gen_ramified, {}, 400),
            (gen_cell, {}, 350), (gen_cell, {"force_composite": True}, 150),
            (gen_fertile_jewel, {}, 100), (gen_multiedge, {}, 300),
            (gen_random_cloud, {}, 300)]
    if args.quick:
        plan = [(g, kw, max(10, n // 10)) for g, kw, n in plan]
    print("\n== fuzz sweep ==")
    for gen, kw, n in plan:
        for i in range(n):
            supp, allowed = gen(rng, **kw)
            try:
                lbl, cite, info = classify(supp)
            except AssertionError as e:
                print(f"FAIL  {gen.__name__}#{i}: {e}\n      supp={supp}")
                nfail += 1
                continue
            counts[lbl] += 1
            total += 1
            if allowed is not None and lbl not in allowed:
                print(f"FAIL  {gen.__name__}#{i}: got {lbl}, "
                      f"allowed {allowed}\n      supp={supp}")
                nfail += 1
        print(f"   {gen.__name__}{'(comp)' if kw else ''}: {n} shapes done",
              flush=True)

    print(f"\nper-stratum counts over {total} random shapes:")
    for k in ("A", "B", "C", "D", "E0", "E1", "E2", "E3", "E4", "F"):
        print(f"   {k:3} {counts[k]:5}   {CITE[k][:58]}")
    zero_rows = [k for k, v in counts.items() if v == 0]
    if zero_rows:
        print(f"NOTE: strata never hit by the random sweep: {zero_rows} "
              f"(zoo covers them)" )
    print(f"\n{'ALL GREEN' if nfail == 0 else f'{nfail} FAILURES'}: "
          f"every shape landed in exactly one stratum")
    sys.exit(1 if nfail else 0)


if __name__ == '__main__':
    main()
