"""Residue-clash verification for cell (3,1,3).

Claim chain:
  (1) The sheared-germ pair  w = -lam*(1/x + (3/2)x^2*(Y+gam*x)^2),
      c = gam*x + Y + x^3*(Y+gam*x)^3  (chart), i.e. in the plane
      c = x*q*(1 + x^5*q^2) with q := gam + x^2*y, solves J(w,c) = lam*x
      exactly, u(0,y) = gam != 0.
  (2) Each fiber {c = t} has exactly three x->0 punctures: the main branch
      (y ~ t x^-3) and a companion pair (y ~ +-i x^-9/2, the sheared
      companion curve {1 + x^5 q^2 = 0}).  Machine: solve the Puiseux
      series order by order and verify c - t vanishes to high order.
  (3) w's x^{-1}-coefficient is -lam at the main puncture and +lam/2 at
      the companions.  Master identity P_W(w,c) = mu/x with
      P_W = 2 p2(C) W + p1(C)  (D = 2 forced by sigma_1 = 3) then gives
      -2 lam p2(t) = mu  and  +lam p2(t) = mu  for generic t, so mu = 0:
      NO SCAFFOLD over any pair with these residues.
  (4) Window-free robustness for arbitrary polynomial tails:
      (a) every admissible knob monomial lies strictly above both Newton
          edges of the fiber hull (hull is cone-rigid: same 3 punctures);
      (b) in the full T11 march (all weight<=12 knobs symbolic), w's only
          pole monomial is (-1,0) with coefficient -lam and the v=-2 line
          is exactly the germ, knob- and E1-free; corrections have v > -2
          and are subleading at all three punctures (size x^{v/2} at the
          companions, x^a with a >= 0 at the main branch);
      (c) kernel additions g(c) are constant on fibers: residues untouched.
  (5) Independent windowed check: the linear system P_W(w,c) = mu/x over
      deg p2 <= 8, deg p1 <= 9 has only the zero solution (mu = 0).

Exit 0 iff every check passes.
"""
import sys
import sympy as sp

x, y, Y, s, t, lam, gam = sp.symbols('x y Y s t lam gam')
q = gam + x**2*y
c_plane = sp.expand(x*q*(1 + x**5*q**2))
w_plane = sp.expand(-lam*(1/x + sp.Rational(3, 2)*x**4*q**2))

results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))

# ---- (1) the pair, once more, from the plane-level closed form ----------
Jp = sp.simplify(sp.diff(w_plane, x)*sp.diff(c_plane, y)
                 - sp.diff(w_plane, y)*sp.diff(c_plane, x))
check("pair equation J(w,c) = lam*x (plane, exact)", sp.simplify(Jp - lam*x) == 0)
u_of = sp.expand(sp.cancel(c_plane/x))
check("legitimacy u(0,y) = gam", sp.simplify(u_of.subs(x, 0) - gam) == 0)

# ---- (2)+(3) Puiseux branches and residues ------------------------------
ORD = 14   # series depth

PAD = 24   # Laurent padding: multiply residuals by param**PAD before Poly

def laurent_dict(expr, param):
    """{exponent: coeff} for a Laurent polynomial in param."""
    P = sp.Poly(sp.expand(expr*param**PAD), param)
    return {e - PAD: co for (e,), co in P.as_dict().items()}

def branch_series(y_lead_pow, u0, param, x_of_s, n_terms=ORD):
    """Solve y(s) = s**y_lead_pow * (u0 + sum u_i s^i) from c(x(s), y) = t,
    order by order (linear in each new u_i). Returns y-series."""
    us = [u0] + [sp.Symbol(f'u{i}') for i in range(1, n_terms)]
    yser = sum(ui*param**i for i, ui in enumerate(us))*param**y_lead_pow
    F = sp.expand(c_plane.subs([(x, x_of_s), (y, yser)]) - t)
    sol = {}
    for _ in range(n_terms - 1):
        terms = sorted(laurent_dict(F.subs(sol), param).items())
        done = True
        for e, co in terms:
            co = sp.expand(co)
            if co == 0:
                continue
            free = [u for u in us[1:] if co.has(u) and u not in sol]
            if not free:
                # must be zero identically or we fail
                if sp.simplify(co) != 0:
                    return None, sol, (e, co)
                continue
            uu = free[0]
            sol[uu] = sp.expand(sp.solve(co, uu)[0])
            done = False
            break
        if done:
            break
    return sp.expand(yser.subs(sol)), sol, None

print("\n-- main branch: x = s, y = s^-3 (t + ...) --")
y_main, sol_m, err_m = branch_series(-3, t, s, s)
check("main branch series solves c = t", err_m is None,
      f"err={err_m}" if err_m else "")
res_main = None
if err_m is None:
    F = sp.expand(c_plane.subs([(x, s), (y, y_main)]) - t)
    lowr = laurent_dict(F, s)
    minexp = min((e for e, co in lowr.items() if sp.expand(co) != 0), default=99)
    check("main branch residual pushed high", minexp >= 8, f"first nonzero at s^{minexp}")
    w_m = sp.expand(w_plane.subs([(x, s), (y, y_main)]))
    wd = laurent_dict(w_m, s)
    minw = min(e for e, co in wd.items() if sp.expand(co) != 0)
    check("main branch w pole order 1", minw == -1, f"min exponent {minw}")
    res_main = sp.expand(wd.get(-1, 0))
    check("main-branch w residue = -lam", sp.simplify(res_main + lam) == 0,
          f"res={res_main}")

print("\n-- companion branches: x = s^2, y = s^-9 (u0 + ...), u0^2 = -1 --")
I = sp.I
res_comp = {}
for u0 in (I, -I):
    y_c, sol_c, err_c = branch_series(-9, u0, s, s**2)
    check(f"companion({u0}) series solves c = t", err_c is None,
          f"err={err_c}" if err_c else "")
    if err_c is None:
        F = sp.expand(c_plane.subs([(x, s**2), (y, y_c)]) - t)
        lowr = laurent_dict(F, s)
        minexp = min((e for e, co in lowr.items() if sp.expand(co) != 0), default=99)
        check(f"companion({u0}) residual pushed high", minexp >= 8,
              f"first nonzero at s^{minexp}")
        w_c = sp.expand(w_plane.subs([(x, s**2), (y, y_c)]))
        wd = laurent_dict(w_c, s)
        minw = min(e for e, co in wd.items() if sp.expand(co) != 0)
        check(f"companion({u0}) w pole order 2 in s (order 1 in x)", minw == -2,
              f"min exponent {minw}")
        rc = sp.expand(wd.get(-2, 0))     # coeff of s^-2 = coeff of x^-1
        res_comp[str(u0)] = rc
        check(f"companion({u0}) w residue = +lam/2",
              sp.simplify(rc - lam/2) == 0, f"res={rc}")

# residue clash => mu = 0
p2t, mu = sp.symbols('p2t mu')
eqs = [sp.Eq(2*p2t*res_main, mu)] + [sp.Eq(2*p2t*rc, mu) for rc in res_comp.values()]
solmu = sp.solve(eqs, [p2t, mu], dict=True)
clash = all(sd.get(mu, 0) == 0 and sd.get(p2t, 0) == 0 for sd in solmu)
check("master identity forces mu = 0 (residue clash)", clash, f"sol={solmu}")

# ---- (4a) Newton-hull rigidity over the FULL knob lattice ---------------
# fiber hull edges (plane exponents (i,j)): edge1 (0,0)-(3,1): j = i/3;
# edge2 (3,1)-(12,3): j - 1 = 2(i-3)/9. A knob x^a Y^b maps to
# (i,j) = (a+3b-1, b). Strictly above edge1 <=> a >= 2; strictly above
# edge2 <=> 2a - 3b >= 0 > -1. Verify for every admissible slot, any weight.
bad = []
for b in range(0, 4):
    for a in range(0, 40):
        wt = 2*a - 3*b
        if wt < 0 or (a, b) in [(0, 0), (0, 1), (3, 3), (1, 0)]:
            continue
        i, j = a + 3*b - 1, b
        above1 = sp.Rational(j, 1) < sp.Rational(i, 3)
        above2 = (j - 1) < sp.Rational(2*(i - 3), 9) if i >= 3 else True
        # for i < 3 edge2 does not bound; edge1 check suffices there
        if not (above1 and above2):
            bad.append((a, b))
check("every admissible knob strictly above both hull edges (hull cone-rigid)",
      not bad, f"violations={bad}")

# ---- (4b) full-march w: pole monomial and v=-2 line are knob/E1-free ----
sys.path.insert(0, '.')
from marcher2 import march, wt as WT, lam as LAM, gam as GAM, E1 as E1sym
records, sub, wdict = march(11, 12, 15, 25, minimal=False, verbose=False)
pole_mons = {m: v for m, v in wdict.items() if m[0] < 0}
check("full march: unique pole monomial (-1,0)", set(pole_mons) == {(-1, 0)},
      f"pole mons={sorted(pole_mons)}")
c_m1 = sp.expand(wdict.get((-1, 0), 0) - E1sym*0)
# E1*c_min injection sits at (0,1),(3,3); germ coefficient check:
check("full march: w[(-1,0)] = -lam (knob/E1-free)",
      sp.simplify(c_m1 + LAM) == 0, f"got {c_m1}")
c_22 = sp.expand(wdict.get((2, 2), 0))
check("full march: w[(2,2)] = -(3/2)lam (knob/E1-free)",
      sp.simplify(c_22 + sp.Rational(3, 2)*LAM) == 0, f"got {c_22}")
v2line_extra = {m: v for m, v in wdict.items()
                if WT(m) == -2 and m not in [(-1, 0), (2, 2)] and sp.expand(v) != 0}
check("full march: v=-2 line is exactly the germ", not v2line_extra,
      f"extras={list(v2line_extra)[:5]}")
below = {m: v for m, v in wdict.items() if WT(m) < -2 and sp.expand(v) != 0
         and m not in [(0, 1), (3, 3)]}   # E1*c_min injection sits at v=-3
check("full march: nothing below v=-2 except the E1 injection", not below,
      f"below={sorted(below)[:5]}")

# ---- (5) windowed P_W emptiness on the closed-form pair -----------------
print("\n-- windowed P_W solve: 2 p2(c) w + p1(c) = mu/x, deg p2<=8, deg p1<=9 --")
P2 = [sp.Symbol(f'P2_{i}') for i in range(9)]
P1 = [sp.Symbol(f'P1_{i}') for i in range(10)]
MU = sp.Symbol('MU')
p2c = sum(co*c_plane**i for i, co in enumerate(P2))
p1c = sum(co*c_plane**i for i, co in enumerate(P1))
G = sp.expand(2*p2c*sp.expand(x*w_plane) + x*p1c - MU)   # x * (identity)
gam0, lam0 = sp.Integer(3), sp.Integer(5)   # generic nonzero rationals
G0 = sp.expand(G.subs([(gam, gam0), (lam, lam0)]))
eqs = sp.Poly(G0, x, y).coeffs()
solP = sp.solve(eqs, P2 + P1 + [MU], dict=True)
empty = bool(solP) and all(sp.simplify(sd.get(MU, MU)) == 0 for sd in solP)
if solP:
    nonzero_p2 = [sp.simplify(sd.get(v, 0)) for sd in solP for v in P2]
    empty = empty and all(vv == 0 for vv in nonzero_p2)
check("windowed P_W solve forces mu = 0 and p2 = 0", empty, f"sol={solP}")

fails = [n for n, ok in results if not ok]
print(f"\n{len(results)} checks, {len(results)-len(fails)} passed, {len(fails)} failed")
sys.exit(1 if fails else 0)
