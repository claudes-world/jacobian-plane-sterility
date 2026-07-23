"""Fat-edge kill (Staircase Sterility) — machine witnesses, 2026-07-23 jc lane.

Claim under test (FATEDGE-KILL): let c be a c-staircase — every non-top
monomial (i,j) of c satisfies omega := N*i - alpha*j > 0, where (alpha, N) is
the top vertex — i.e. the fiber hull of c - t is the single fat edge
(0,0)->(alpha,N), N >= 2. Then NO line pair (w = A/x^d, d >= 1, x !| A,
J(w,c) = lam*x^m, m >= 1, lam != 0) exists, for ANY above-hull polynomial
tail, any multiplicities, any (d, m).

Proof skeleton being witnessed:
  (K1) bracket grading: omega(J(M1, M2)) = omega(M1) + omega(M2) + (alpha - N),
       and omega(c-terms) >= 0 with equality ONLY at the top vertex;
  (K2) the leading bracket J(., x^alpha y^N) shifts y-degree by N-1 >= 1, so
       it can never source lam*x^m (y-degree 0) — cokernel functional = [y^0];
  (K3) kernel of J(., x^alpha y^N) = monomials on the ray N*p = alpha*j,
       all at omega = 0;
  (K4) bottom-up: the minimal omega-level of w is killed into the kernel,
       hence omega >= 0 on supp(w); poles (p < 0, j >= 0) force omega < 0.
       Contradiction with d >= 1. No pair.

Witness layers:
  A. exact identity checks for K1-K3 (both flagship staircases + chart form);
  B. windowed pair hunts with FULLY SYMBOLIC tails over QQ(knobs):
     lam forced to 0 (2-edge staircase; genuine T6 conspiracy with symbolic
     u, v + knobs); plus adversarial + random knob specializations;
  C. controls the kill must NOT touch: split cell (3,1,3) sheared pair is
     FOUND by the same hunter (lam free); T7 jewel's fiber polygon is
     genuinely two-edged (hypothesis fails there);
  D. mini fat-edge march in the chart (Y = x^3 y): watch the terminal
     obstruction fire at bucket m-3 with condition -lam, no shield reachable.
"""
import time
import sympy as sp

x, y, t = sp.symbols('x y t')
lam = sp.Symbol('lam')
u, v = sp.symbols('u v')

npass = nfail = 0
def check(label, ok, extra=""):
    global npass, nfail
    tag = "PASS" if ok else "FAIL"
    if ok: npass += 1
    else: nfail += 1
    print(f"{tag}  {label}" + (f"  [{extra}]" if extra else ""))

def J(f, g):
    return sp.expand(sp.diff(f, x) * sp.diff(g, y) - sp.diff(f, y) * sp.diff(g, x))

def monomials(expr):
    """dict {(i,j): coeff} of a polynomial in x, y (other symbols in coeffs)."""
    p = sp.Poly(sp.expand(expr), x, y)
    return {(mi, mj): co for (mi, mj), co in zip(p.monoms(), p.coeffs())}

# ---------------------------------------------------------------- staircases
# 2-edge simple staircase (slopes 3 > 2), top vertex (6,2), N=2, alpha=6
C2 = x + x**4*y + x**6*y**2
# genuine T6 conspiracy (N1=1, slopes 4>3>2, double+double), top (15,5), N=5
C5 = x + x**5*y + 2*u*x**8*y**2 + u**2*x**11*y**3 + 2*u*v*x**13*y**4 + v**2*x**15*y**5

def top_vertex(c):
    ms = monomials(c)
    N = max(j for (_, j) in ms)
    alpha = min(i for (i, j) in ms if j == N)
    return alpha, N

def omega(alpha, N, i, j):
    return N*i - alpha*j

# ---- A. exact identity checks -------------------------------------------
print("== A. graded-structure identities ==")
for cname, c in (("2-edge", C2), ("conspiracy", C5)):
    alpha, N = top_vertex(c)
    ms = monomials(c)
    bad = [(i, j) for (i, j) in ms if (i, j) != (alpha, N) and omega(alpha, N, i, j) <= 0]
    check(f"{cname}: fiber collapse (omega>0 off top vertex, alpha={alpha}, N={N})",
          not bad, f"bad={bad}")

# K1: bracket omega-shift, exact on random monomials (2-edge parameters)
alpha, N = 6, 2
ok = True
import random
random.seed(313)
for _ in range(40):
    p1, j1 = random.randint(-6, 8), random.randint(0, 5)
    p2, j2 = random.randint(0, 9), random.randint(0, 5)
    br = J(x**p1*y**j1, x**p2*y**j2)
    if br == 0:
        continue
    mons = sp.Poly(sp.expand(br*x**12*y**12), x, y).monoms()  # shift to clear negatives
    for (mi, mj) in mons:
        i, j = mi - 12, mj - 12
        if omega(alpha, N, i, j) != omega(alpha, N, p1, j1) + omega(alpha, N, p2, j2) + (alpha - N):
            ok = False
check("K1: omega(J(M1,M2)) = omega1 + omega2 + (alpha-N) on random monomials", ok)

# K2: leading bracket shifts y-degree by exactly N-1
ok = True
for cname, c in (("2-edge", C2), ("conspiracy", C5)):
    alpha, N = top_vertex(c)
    topmon = x**alpha*y**N
    for _ in range(20):
        p1, j1 = random.randint(-6, 8), random.randint(0, 5)
        br = J(x**p1*y**j1, topmon)
        if br == 0:
            continue
        pb = sp.Poly(sp.expand(br*x**20*y**2), x, y)
        ydegs = {mj - 2 for (_, mj) in pb.monoms()}
        if ydegs != {j1 + N - 1}:
            ok = False
check("K2: J(., top vertex) shifts y-degree by N-1 >= 1 (never reaches y^0)", ok)

# K3: kernel of J(., x^alpha y^N) = ray N*p = alpha*j (omega = 0)
ok = True
for (alpha, N) in ((6, 2), (15, 5)):
    for p1 in range(-8, 9):
        for j1 in range(0, 6):
            k = J(x**p1*y**j1, x**alpha*y**N)
            in_ker = (k == 0)
            on_ray = (N*p1 == alpha*j1)
            if in_ker != on_ray:
                ok = False
check("K3: ker J(., x^a y^N) == ray N*p = alpha*j (subset of omega=0)", ok)

# chart form of the 2-edge staircase: Y = x^3 y  =>  c = Y^2 + x + x*Y
Ych = sp.Symbol('Y')
c_chart = sp.expand(C2.subs(y, Ych/x**3))
check("chart: c(x, Y/x^3) == Y^2 + x + x*Y exactly",
      sp.expand(c_chart - (Ych**2 + x + x*Ych)) == 0, f"got {c_chart}")

# ---- B. windowed symbolic-tail pair hunts -------------------------------
print("== B. windowed pair hunts, symbolic tails: lam must be forced to 0 ==")

def hull_knob_slots(c, extra=6):
    """Above-hull tail slots (i,j), j <= N, modest menu near the hull."""
    alpha, N = top_vertex(c)
    ms = set(monomials(c))
    slots = []
    for j in range(0, N + 1):
        for i in range(0, alpha + extra + 1):
            if (i, j) in ms:
                continue
            om = omega(alpha, N, i, j)
            if om > 0 and above_hull(c, i, j):
                slots.append((i, j))
    return slots

def above_hull(c, i, j):
    """strictly above the c-polygon lower-left boundary (keeps Newton data):
    require (i,j) strictly above every supporting edge line of the c-hull,
    approximated by: strictly above each consecutive-vertex edge, and not
    left of the bottom vertex / below the top vertex column."""
    ms = monomials(c)
    verts = sorted({(min(i2 for (i2, j2) in ms if j2 == jj), jj)
                    for jj in {j2 for (_, j2) in ms}})
    # consecutive edges (low y -> high y)
    for (i1, j1), (i2, j2) in zip(verts, verts[1:]):
        # strictly above edge through (i1,j1),(i2,j2): cross product test
        cr = (i - i1) * (j2 - j1) - (j - j1) * (i2 - i1)
        if j1 <= j <= j2 and cr <= 0:
            return False
    if j == 0 and i <= verts[0][0]:
        return False
    if j == verts[-1][1] and i <= verts[-1][0]:
        return False
    return True

def pair_hunt(c, d, m, awin, bwin, knob_syms=(), label="", quiet=False):
    """Solve x*J(A,c) - d*A*c_y = lam*x^(m+d+1) for polynomial A in window.
    Returns (lam_forced_zero, nullspace_dim, witness) over QQ(knobs)."""
    A_syms = {}
    A = 0
    for i in range(awin + 1):
        for j in range(bwin + 1):
            s = sp.Symbol(f'A_{i}_{j}')
            A_syms[(i, j)] = s
            A += s * x**i * y**j
    eq = sp.expand(x*J(A, c) - d*A*sp.diff(c, y) - lam*x**(m + d + 1))
    p = sp.Poly(eq, x, y)
    unknowns = list(A_syms.values()) + [lam]
    M = sp.Matrix([[co.coeff(s) for s in unknowns] for co in p.coeffs()])
    ns = M.nullspace()
    lam_comps = [sp.simplify(vec[-1]) for vec in ns]
    forced = all(lc == 0 for lc in lam_comps)
    return forced, len(ns), lam_comps

# B1: 2-edge staircase, 2 SYMBOLIC knobs (tractable fraction field), all (d,m)
k1s, k2s = sp.symbols('kk1 kk2')
slots2 = hull_knob_slots(C2)
C2k2 = C2 + k1s*x**5*y + k2s*x**3          # one mid-edge knob + one low knob
allok = True
for (d, m) in ((1, 1), (1, 2), (2, 1), (1, 3)):
    t0 = time.time()
    forced, nsd, _ = pair_hunt(C2k2, d, m, 11, 4)
    allok = allok and forced
    print(f"   2-edge sym 2-knob (d={d},m={m}): lam forced 0 = {forced} "
          f"(null dim {nsd}, {round(time.time()-t0,1)}s)", flush=True)
check("2-edge + SYMBOLIC knobs (x^5 y, x^3): lam == 0 for all (d,m) tried", allok)

# B2: 6-knob family at specializations (random + adversarial), fast QQ solves
k_syms = sp.symbols('k1:7')
slots6 = []
for jj in (0, 1, 2):                       # spread across heights, 2 per row
    slots6 += [s for s in slots2 if s[1] == jj][:2]
slots6 = slots6[:6]
C2k6 = C2 + sum(ks * x**i * y**j for ks, (i, j) in zip(k_syms, slots6))
allok = True
specs = [{k: sp.Rational(random.randint(-9, 9), random.randint(1, 5)) for k in k_syms}
         for _ in range(3)]
specs += [{k: 1 for k in k_syms}, {k: (3 if n % 2 else -2) for n, k in enumerate(k_syms)}]
t0 = time.time()
for spec in specs:
    for (d, m) in ((1, 1), (2, 1)):
        forced, _, _ = pair_hunt(C2k6.subs(spec), d, m, 12, 4)
        allok = allok and forced
check(f"2-edge 6-knob family {slots6}, 5 specializations x (d,m) grid: lam == 0",
      allok, f"{round(time.time()-t0,1)}s")

# B3: genuine conspiracy — mod-p rank certificates over many parameter points
# lam forced 0 at a point  <=>  rank([M_A | m_lam]) = rank(M_A) + 1 there.
# One point with rank_p(M_A) at its observed max and augmented rank +1 certifies
# generic forcing; a sweep over primes x random points makes special-locus
# escape astronomically unlikely. (The PROOF is the graded argument; this layer
# is the falsifier.)
import numpy as np

def modp_rank(Mint, p):
    A = np.array(Mint, dtype=np.int64) % p
    r = 0
    rows, cols = A.shape
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i, c] % p), None)
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        for i in range(rows):
            if i != r and A[i, c]:
                A[i] = (A[i] - A[i, c] * A[r]) % p
        r += 1
        if r == rows:
            break
    return r

def build_system(c, d, m, awin, bwin, syms):
    A_syms, A = {}, 0
    for i in range(awin + 1):
        for j in range(bwin + 1):
            s = sp.Symbol(f'A_{i}_{j}')
            A_syms[(i, j)] = s
            A += s * x**i * y**j
    eq = sp.expand(x*J(A, c) - d*A*sp.diff(c, y) - lam*x**(m + d + 1))
    p = sp.Poly(eq, x, y)
    unknowns = list(A_syms.values()) + [lam]
    rows = [[co.coeff(s) for s in unknowns] for co in p.coeffs()]
    return sp.Matrix(rows), unknowns

kc1, kc2 = sp.symbols('kc1 kc2')
slots5 = [s for s in hull_knob_slots(C5) if s[1] in (1, 2)][:2]
C5k = C5 + kc1*x**slots5[0][0]*y**slots5[0][1] + kc2*x**slots5[1][0]*y**slots5[1][1]
PRIMES = (1000003, 2000003)
t0 = time.time()
allok, worst = True, []
for (d, m) in ((1, 1), (1, 2), (2, 1), (2, 2)):
    Msym, unk = build_system(C5k, d, m, 18, 7, (u, v, kc1, kc2))
    lam_col = Msym.shape[1] - 1
    for p in PRIMES:
        for trial in range(12):
            spec = {u: random.randint(1, p - 1), v: random.randint(1, p - 1),
                    kc1: random.randint(0, p - 1), kc2: random.randint(0, p - 1)}
            Mnum = Msym.subs(spec)
            Mint = [[int(e) for e in row] for row in Mnum.tolist()]
            rA = modp_rank([r[:lam_col] for r in Mint], p)
            rAug = modp_rank(Mint, p)
            if rAug != rA + 1:
                allok = False
                worst.append((d, m, p, spec))
check("conspiracy + sym knobs, mod-p certificates: rank jump +1 at ALL "
      "96 (d,m,p,point) samples (win 18,7)", allok,
      f"{round(time.time()-t0,1)}s bad={worst[:2]}")

# B4: conspiracy exact QQ solves at adversarial/special loci (small window)
t0 = time.time()
allok = True
for spec in ({u: 1, v: 1, kc1: 0, kc2: 0}, {u: 1, v: 1, kc1: 2, kc2: -3},
             {u: -1, v: 1, kc1: 1, kc2: 1}, {u: 2, v: 4, kc1: 0, kc2: 0}):
    forced, _, _ = pair_hunt(C5k.subs(spec), 1, 1, 15, 6)
    allok = allok and forced
check("conspiracy exact QQ at 4 special loci (incl. v=u^2 resonance): lam == 0",
      allok, f"{round(time.time()-t0,1)}s")

# ---- C. controls (the kill must NOT fire here) --------------------------
print("== C. controls ==")
# C1: split cell (3,1,3) sheared pair EXISTS and hunter finds it (gamma=1)
q = 1 + x**2*y
Csplit = sp.expand(x*q*(1 + x**5*q**2))
t0 = time.time()
forced, nsd, lam_comps = pair_hunt(Csplit, 1, 1, 14, 4)
check("split-cell control: hunter FINDS pair (lam NOT forced 0)",
      (not forced) and nsd >= 1, f"null dim {nsd}, {round(time.time()-t0,1)}s")
# C2: hypothesis fails on split cell: some non-top monomial has omega <= 0
alpha, N = top_vertex(Csplit)
ms = monomials(Csplit)
viol = [(i, j) for (i, j) in ms if (i, j) != (alpha, N) and omega(alpha, N, i, j) <= 0]
check("split-cell control: fat-edge hypothesis correctly FAILS (omega<=0 exists)",
      len(viol) > 0, f"violations={viol[:4]}")
# C3: T7 jewel u — fiber polygon genuinely 2-edged (not collapsed)
Ujewel = 1 + x**3 + (x**2 + 2*x**5)*y + x**7*y**2
alpha, N = top_vertex(Ujewel)
msj = monomials(Ujewel)
viol = [(i, j) for (i, j) in msj if (i, j) != (alpha, N) and j > 0 and omega(alpha, N, i, j) <= 0]
check("jewel control: hypothesis correctly FAILS on T7 jewel (collapse violated)",
      len(viol) > 0, f"violations={viol}")

# ---- D. mini fat-edge march in the chart --------------------------------
print("== D. chart march: watch the terminal obstruction ==")
# c_chart = Y^2 + x + x*Y (+ chart knobs, all with x-shift alpha_c >= 1).
# march levels: w = sum_{a} x^a * phi_a(Y). Bucket beta equation:
#   2(beta+1) * Y * phi_{beta+1} + [beta*phi_beta - (1+Y)*phi_beta'] + knobs = lam*delta(beta, m-3)
# bottom-up from below: all phi_a = 0 for a < 0 forced; bucket m-3 reads 0 = lam.
def chart_march(mm, bwin=6, floor=-8):
    """Build the linear system over the negative chart levels; the claim is
    that its ONLY solution is zero (all negative levels dead, lam dead)."""
    phis = {}          # a -> list of Y-coeff symbols (length bwin+1)
    for a in range(floor, 0):
        phis[a] = [sp.Symbol(f'p{a}_{b}') for b in range(bwin + 1)]
    unknowns = [s for a in sorted(phis) for s in phis[a]] + [lam]
    conds = []
    # buckets floor-1 .. -1: contributions L0[w_{beta+1}] + L1[w_beta] only
    # (level-0 w has L0 == 0 and feeds no negative bucket; chart knobs all
    # shift x-level UP so they feed no negative bucket from dead levels)
    for beta in range(floor - 1, 0):
        eq = {}        # Y-degree -> expr
        nxt = phis.get(beta + 1)
        if nxt is not None:
            for b in range(bwin + 1):       # L0: 2(beta+1) Y^{b+1} phi_{beta+1,b}
                eq[b + 1] = eq.get(b + 1, 0) + 2*(beta + 1)*nxt[b]
        cur = phis.get(beta)
        if cur is not None:                  # L1: (beta - b) phi_b - (b) phi_b -> Y^{b-1}
            for b in range(bwin + 1):
                eq[b] = eq.get(b, 0) + (beta - b)*cur[b]
                if b >= 1:
                    eq[b - 1] = eq.get(b - 1, 0) - b*cur[b]
        rhs = {0: lam} if beta == mm - 3 else {}
        for b in sorted(set(eq) | set(rhs)):
            conds.append(sp.expand(eq.get(b, 0) - rhs.get(b, 0)))
    M = sp.Matrix([[sp.diff(cnd, s) for s in unknowns] for cnd in conds])
    return M, unknowns

for mm in (1, 2):
    M, unknowns = chart_march(mm)
    ns = M.nullspace()
    check(f"chart march m={mm}: nullspace of negative-level system is ZERO "
          f"(all levels + lam killed)", len(ns) == 0,
          f"{M.shape[0]} eqs, {len(unknowns)} unknowns")

# ---- E. Theorem' (Ramified-Edge Kill): bare-ramified FIRST edge ----------
print("== E. ramified-edge generalization (Theorem') ==")
# convex two-edge shape: first edge (0,0)-(2,2) bare-ramified (N1=2),
# second edge to (4,3); NOT collapsed w.r.t. top vertex (omega(2,2) = -2).
CE = x + x**2*y**2 + x**4*y**3
alpha, N = top_vertex(CE)
ms = monomials(CE)
check("E-shape genuinely outside Theorem (omega<=0 off top exists)",
      any((i, j) != (alpha, N) and omega(alpha, N, i, j) <= 0 for (i, j) in ms),
      f"top={alpha, N}")
# omega1-grading: only the first-edge endpoint at level 0, all else > 0
i1, N1 = 2, 2
w1 = {m: N1*m[0] - i1*m[1] for m in ms}
check("E-shape: omega1 = 0 only at first-edge endpoint, > 0 elsewhere",
      w1[(2, 2)] == 0 and all(v > 0 for m, v in w1.items() if m != (2, 2)),
      f"{w1}")
CEt = CE + 3*x**3*y + 2*x**3*y**2 + 5*x**5*y**3
allok = True
for cname, c in (("bare", CE), ("tailed", CEt)):
    for (d, m) in ((1, 1), (1, 2), (2, 1)):
        forced, _, _ = pair_hunt(c, d, m, 12, 5)
        allok = allok and forced
check("E-shape bare + tailed, (d,m) grid: lam == 0 (Theorem' kill)", allok)

print(f"\n{npass + nfail} checks, {npass} passed, {nfail} failed")
if nfail:
    raise SystemExit(1)
