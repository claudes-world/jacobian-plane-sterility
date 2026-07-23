"""Explicit residue-clash witnesses for cell (2,3) (k = 2, sigma = 5, m = 2),
mirroring verify_residue.py at (3,1,3), plus the trace check and the hull
sweep extended to N = 5, 6.

Sheared germ pair at (2,3):
  chart:  c = c_min(x, Yt), w = x^-2 * phi(zeta_t),  Yt = Y + gam*x,
          phi(z) = -lam*(1/2 + 3z + (9/4)z^2),  zeta_t = x^3*Yt^2
  plane:  Y = x^5*y.
Predictions: J(w,c) = lam*x^-3 chart / lam*x^2 plane (m = sigma-1-k = 2);
main-puncture x^-2-coefficient phi(0) = -lam/2; companion (one cycle, ram 2)
x^-2-coefficient phi(-1) = +lam/4; master identity P_W = 2p2(C)W + p1(C) with
P_W(w,c) = mu*x^-2 forces mu = 0; Tr_c(w) is x-free.
"""
import sys
import sympy as sp

x, y, Y, s, t, lam, gam = sp.symbols('x y Y s t lam gam')
z = sp.Symbol('z')

phi = -lam*(sp.Rational(1, 2) + 3*z + sp.Rational(9, 4)*z**2)
Yt = Y + gam*x
c_chart = Yt + x**3*Yt**3
w_chart = sp.expand((x**-2*phi).subs(z, x**3*Yt**2))
c_plane = sp.expand(c_chart.subs(Y, x**5*y))
w_plane = sp.expand(w_chart.subs(Y, x**5*y))

results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))

Jc = sp.simplify(sp.diff(w_chart, x)*sp.diff(c_chart, Y) - sp.diff(w_chart, Y)*sp.diff(c_chart, x))
check("(2,3) chart pair eq = lam x^-3", sp.simplify(Jc - lam*x**-3) == 0)
Jp = sp.simplify(sp.diff(w_plane, x)*sp.diff(c_plane, y) - sp.diff(w_plane, y)*sp.diff(c_plane, x))
check("(2,3) plane pair eq = lam x^2 (m = 2)", sp.simplify(Jp - lam*x**2) == 0)
u_of = sp.expand(sp.cancel(c_plane/x))
check("legitimacy u(0,y) = gam", sp.simplify(u_of.subs(x, 0) - gam) == 0)

PAD = 40
def laurent_dict(expr, param):
    P = sp.Poly(sp.expand(expr*param**PAD), param)
    return {e - PAD: co for (e,), co in P.as_dict().items()}

def branch_series(y_lead_pow, u0, param, x_of_s, n_terms=16):
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
                if sp.simplify(co) != 0:
                    return None, sol, (e, co)
                continue
            sol[free[0]] = sp.expand(sp.solve(co, free[0])[0])
            done = False
            break
        if done:
            break
    return sp.expand(yser.subs(sol)), sol, None

print("\n-- main branch: x = s, y ~ t s^-5 --")
y_m, _, err = branch_series(-5, t, s, s)
check("main series solves c = t", err is None, f"err={err}" if err else "")
if err is None:
    wd = laurent_dict(sp.expand(w_plane.subs([(x, s), (y, y_m)])), s)
    minw = min(e for e, co in wd.items() if sp.expand(co) != 0)
    check("main pole order 2", minw == -2, f"min {minw}")
    r_main = sp.expand(wd.get(-2, 0))
    check("main residue = phi(0) = -lam/2", sp.simplify(r_main + lam/2) == 0, f"{r_main}")

print("\n-- companion: x = s^2, y ~ u0 s^-13, u0^2 = -1 --")
r_comps = []
for u0 in (sp.I, -sp.I):
    y_c, _, err = branch_series(-13, u0, s, s**2)
    check(f"companion({u0}) series solves c = t", err is None, f"err={err}" if err else "")
    if err is None:
        wd = laurent_dict(sp.expand(w_plane.subs([(x, s**2), (y, y_c)])), s)
        minw = min(e for e, co in wd.items() if sp.expand(co) != 0)
        check(f"companion({u0}) pole order 4 in s (2 in x)", minw == -4, f"min {minw}")
        rc = sp.expand(wd.get(-4, 0))
        r_comps.append(rc)
        check(f"companion({u0}) residue = phi(-1) = +lam/4",
              sp.simplify(rc - lam/4) == 0, f"{rc}")

print("\n-- windowed P_W solve: 2 p2(c) w + p1(c) = mu x^-2, deg <= 5/6 --")
MU = sp.Symbol('MU')
P2 = [sp.Symbol(f'P2_{i}') for i in range(6)]
P1 = [sp.Symbol(f'P1_{i}') for i in range(7)]
gam0, lam0 = sp.Integer(3), sp.Integer(5)
cp = sp.expand(c_plane.subs([(gam, gam0), (lam, lam0)]))
wp = sp.expand(w_plane.subs([(gam, gam0), (lam, lam0)]))
p2c = sum(co*cp**i for i, co in enumerate(P2))
p1c = sum(co*cp**i for i, co in enumerate(P1))
G = sp.expand(2*p2c*sp.expand(x**2*wp) + x**2*p1c - MU)
eqs = sp.Poly(G, x, y).coeffs()
solP = sp.solve(eqs, P2 + P1 + [MU], dict=True)
empty = bool(solP) and all(sp.simplify(v) == 0 for sd in solP for v in sd.values())
check("windowed P_W forces mu = 0, p2 = 0", empty, f"nsol={len(solP)}")

print("\n-- trace check: Tr_c(w) x-free (full Newton identities, deg_y w = 4) --")
p = sp.Poly(cp - t, y)
Np = p.degree()
ac = [p.coeff_monomial(y**i) for i in range(Np + 1)]
ee = [sp.Integer(1)] + [sp.cancel((-1)**k*ac[Np - k]/ac[Np]) for k in range(1, Np + 1)]
ps = [sp.Integer(Np)]
for k in range(1, 9):
    ssum = sum((-1)**(i - 1)*ee[i]*ps[k - i] for i in range(1, min(k - 1, Np) + 1))
    if k <= Np:
        ssum += (-1)**(k - 1)*k*ee[k]
    ps.append(sp.cancel(ssum))
wpol = sp.Poly(wp, y)
Tr = sp.simplify(sp.cancel(sum(wpol.coeff_monomial(y**i)*ps[i]
                               for i in range(wpol.degree() + 1))))
check("Tr(w) x-free", not Tr.has(x), f"Tr={Tr}")

print("\n-- hull-inequality equivalence, N = 5, 6 --")
bad = []
for (N, r) in [(5, 5), (5, 10), (6, 6)]:
    for sigma in range(2, 10):
        for a in range(0, 26):
            for b in range(0, N + 1):
                if (a, b) in [(0, 0), (1, 0), (0, 1), (r, N)]:
                    continue
                v = (N - 1)*a - r*b
                i, j = a + sigma*b - 1, b
                ab1 = j*(sigma - 1) < i
                di, dj = r + sigma*(N - 1), N - 1
                ab2 = (i - (sigma - 1))*dj - (j - 1)*di > 0
                if (ab1 and ab2) != ((a + b >= 2) and (v > -r)):
                    bad.append((N, r, sigma, a, b))
check("hull inequalities hold for N = 5, 6 cells", not bad, f"bad={bad[:4]}")

fails = [n for n, ok in results if not ok]
print(f"\n{len(results)} checks, {len(results)-len(fails)} passed, {len(fails)} failed")
sys.exit(1 if fails else 0)
