"""Composite-cone stratum probe (assembly session, 2026-07-23 ~03:20 ET).

Stratum: split-edge chamber (N1 = 1, two fiber edges), cone edge with
d := gcd(N-1, r) > 1 and SUPPORTED interior lattice points — the cone edge
polynomial is composite. With xi := x^(r/d) Y^((N-1)/d), the cone is
c_min = Y * E(xi), deg E = d. FERRY rev-2.1 treats only the binomial case
E = 1 + xi^d (no interior support).

Questions:
 Q1 (fertility): does the germ equation J(x^-k phi, c_min) = lam x^-(1+k)
    admit a polynomial solution phi when E carries a symbolic interior
    coefficient t1? Compare with the binomial fertility law r | Nk & r !| k.
 Q2 (pairs): windowed hunts on plane realizations, symbolic t1 + numeric
    specials, bare + tails: does anything live?
 Q3 (residues): if germs exist, compute phi at the roots of the cone-fiber
    components (xi-roots of E) — the multi-root residue data for the clash.

Test cell family: (N, r) = (3, 6), d = 2, xi = x^3 Y, interior point (3,2)
= xi*Y. c_min = Y + t1 x^3 Y^2 + x^6 Y^3 = Y*E(xi), E = 1 + t1 xi + xi^2.
k over the binomial-fertile range {1, 2, 4, 5} (r | 3k fails at k=1,2,4,5?
3k: r=6 | 3k iff k even; binomial-fertile k: r|Nk & r!|k: k=2: 6|6 yes,
6!|2 yes -> fertile; k=4: 6|12, 6!|4 -> fertile; k=1: 6|3? no -> barren;
k=5: 6|15? no -> barren). So compare k = 2, 4 (binomial-fertile) with the
trinomial cone.
"""
import time
import sympy as sp

x, Y = sp.symbols('x Y')
lam, t1 = sp.symbols('lam t1')

npass = nfail = 0
def check(label, ok, extra=""):
    global npass, nfail
    tag = "PASS" if ok else "FAIL"
    if ok: npass += 1
    else: nfail += 1
    print(f"{tag}  {label}" + (f"  [{extra}]" if extra else ""))

def Jch(f, g):
    return sp.expand(sp.diff(f, x)*sp.diff(g, Y) - sp.diff(f, Y)*sp.diff(g, x))

N, r = 3, 6
c_min_tri = Y + t1*x**3*Y**2 + x**6*Y**3

# ---- Q1: germ existence, trinomial cone, symbolic t1 --------------------
print("== Q1: trinomial-cone germ / fertility ==")
for k in (1, 2, 3, 4):
    # ansatz: w = x^-k * phi with phi a polynomial in xi = x^3 Y of deg <= D
    D = 8
    coeffs = sp.symbols(f'f0:{D+1}')
    xi = x**3 * Y
    phi = sum(coeffs[n] * xi**n for n in range(D + 1))
    w = x**(-k) * phi
    eq = sp.expand(Jch(w, c_min_tri) - lam * x**(-1-k))
    p = sp.Poly(sp.expand(eq * x**(k+1) * Y**1), x, Y)   # clear negatives
    unk = list(coeffs) + [lam]
    M = sp.Matrix([[co.coeff(s) for s in unk] for co in p.coeffs()])
    ns = M.nullspace()
    lam_free = any(sp.simplify(v[-1]) != 0 for v in ns)
    print(f"   k={k}: germ with lam != 0 exists (symbolic t1): {lam_free} "
          f"(null dim {len(ns)})")

# also: does the germ ansatz need HALF-steps? xi-degree vs zeta-degree:
# binomial germ used zeta = xi^2 only (even xi-powers). Q: with t1 != 0 do
# odd xi-powers activate? (checked implicitly above by full-xi ansatz)

# ---- Q2: windowed pair hunts, plane realization -------------------------
print("== Q2: plane pair hunts on the composite-cone shape ==")
# realization sigma: first edge (0,0)-(sigma,1); minimal sigma with the cone
# edge landing right: chart c = gamma*x + c_min_tri, plane via Y = x^sigma y,
# c_plane = x*u. Use sigma = 3 (mirrors (3,2,6) realization sigma=4? keep
# simple: try sigma = 4 like the (3,2,6) cell).
def plane_c(sigma, t1val, gam=1):
    cch = gam*x + c_min_tri.subs(t1, t1val)
    return sp.expand(cch.subs(Y, x**sigma * sp.Symbol('y')))

y = sp.Symbol('y')
def pair_hunt(c, d, m, awin, bwin):
    A_syms, A = {}, 0
    for i in range(awin + 1):
        for j in range(bwin + 1):
            s_ = sp.Symbol(f'A_{i}_{j}')
            A_syms[(i, j)] = s_
            A += s_ * x**i * y**j
    eq = sp.expand(x*(sp.diff(A, x)*sp.diff(c, y) - sp.diff(A, y)*sp.diff(c, x))
                   - d*A*sp.diff(c, y) - lam*x**(m + d + 1))
    p = sp.Poly(eq, x, y)
    unk = list(A_syms.values()) + [lam]
    M = sp.Matrix([[co.coeff(s_) for s_ in unk] for co in p.coeffs()])
    ns = M.nullspace()
    return all(sp.simplify(v[-1]) == 0 for v in ns), len(ns)

t0 = time.time()
allok = True
for t1val in (sp.Symbol('t1'), 1, -2, sp.Rational(1, 2)):
    c = plane_c(4, t1val)
    for (d, m) in ((1, 1), (2, 1), (2, 2)):
        forced, nsd = pair_hunt(c, d, m, 14, 5)
        tag = 'sym' if t1val == sp.Symbol('t1') else str(t1val)
        if not forced:
            allok = False
            print(f"   !! t1={tag} (d={d},m={m}): lam FREE (null {nsd})")
    print(f"   t1={tag}: all (d,m) forced lam=0 so far ({round(time.time()-t0,1)}s)", flush=True)
check("Q2: composite-cone plane hunts (sym t1 + 3 specials): lam == 0", allok)

# ---- Q3: distinct-root residue data (only meaningful if Q1 finds germs) --
print("== Q3: root structure of E(xi) = 1 + t1 xi + xi^2 ==")
xi_s = sp.Symbol('xi_s')
E = 1 + t1*xi_s + xi_s**2
disc = sp.discriminant(E, xi_s)
print(f"   disc(E) = {disc}  (distinct roots iff t1^2 != 4; t1 = ±2 -> double root, No-Conspiracy territory)")

print(f"\n{npass + nfail} checks, {npass} passed, {nfail} failed")
