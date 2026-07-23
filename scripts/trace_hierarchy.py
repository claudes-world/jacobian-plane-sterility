"""Trace hierarchy for line pairs — the T6 attack toolkit.

Objects: line pair (w, c), J(w,c) = lam*x^m, w polynomial in y over C[x,1/x]
(poles only on x = 0), deg_y c = N >= 2. On the fiber {c = t} with sheets
y_i(x, t):

  Tr(F)     := sum_i F(x, y_i)                (Galois trace over the fibration)
  S_j       := Tr(w^j / c_y)

Identities to establish (each machine-verified here):

  (I1)  d/dx Tr(w^j) = j * Tr(w^{j-1} * J(w,c) / c_y)   [chain rule; any (w,c)]
        For a pair: d/dx Tr(w^j) = j*lam*x^m * S_{j-1}.
  (I2)  S_j = -res_{y=inf} ( w^j dy / (c - t) )
        [all finite y-poles of the integrand are the fiber roots, since w is
         y-polynomial; sum of residues = 0 on P^1_y]
        In particular S_0 = 0 (Trace Lemma engine), and S_j depends ONLY on
        the top part of the y-expansions of w and c — staircase TOP-edge data.
  (I3)  D = 3 forced profile: trace the master identity
        3 p3(c) w^2 + 2 p2(c) w + p1(c) = mu x^{-m}, then d/dx:
            3 p3(t) * d/dx Tr(w^2) = -m N mu x^{-m-1}
            => S_1 = -m N mu / (6 lam p3(t)) * x^{-2m-1}.
        So a D = 3 scaffold FORCES the infinity-residue stream S_1 to the
        exact monomial profile x^{-2m-1} with t-dependence -mNmu/(6 lam p3(t)),
        a rational function of t whose poles are the roots of p3.
        But by (I2), S_1 is (for each fixed x) a polynomial-in-t divided by
        nothing — expanding res_inf in powers of t shows S_1 in C[x,1/x][[t]]
        with t-degree growth controlled by top-edge slopes. Poles in t at
        roots of p3 are impossible unless p3 is constant. p3 constant is the
        first structural narrowing; the x-valuation of S_1 then gives the
        second equation. This file verifies (I1), (I2) exactly on the sheared
        (3,1,3) pair and on a random NON-pair (w, c) [identity (I1) in its
        general form], and computes S_1, S_2 symbolically for the sheared pair
        as a sanity base point.
"""
import sys
import sympy as sp

x, y, t, lam, gam = sp.symbols('x y t lam gam')

results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))

def power_sums(c_expr, kmax):
    """Power sums p_k(x, t) of the fiber roots y_i of c(x,y) = t, exact."""
    p = sp.Poly(sp.expand(c_expr) - t, y)
    Np = p.degree()
    ac = [p.coeff_monomial(y**i) for i in range(Np + 1)]
    ee = [sp.Integer(1)] + [sp.cancel((-1)**k*ac[Np - k]/ac[Np]) for k in range(1, Np + 1)]
    ps = [sp.Integer(Np)]
    for k in range(1, kmax + 1):
        ssum = sum((-1)**(i - 1)*ee[i]*ps[k - i] for i in range(1, min(k - 1, Np) + 1))
        if k <= Np:
            ssum += (-1)**(k - 1)*k*ee[k]
        ps.append(sp.cancel(ssum))
    return ps, Np

def trace_of(F_expr, c_expr, kmax_hint=12):
    """Tr(F) for F polynomial in y (coeffs rational in x, t after substitution):
    express via power sums of the fiber roots."""
    Fp = sp.Poly(sp.expand(F_expr), y)
    d = Fp.degree()
    ps, Np = power_sums(c_expr, max(d, 1))
    return sp.cancel(sum(Fp.coeff_monomial(y**i)*ps[i] for i in range(d + 1)))

def res_inf(F_expr, c_expr):
    """-res_{y=inf}(F dy / (c - t)) computed by series at y = infinity:
    substitute y = 1/v, expand F/(c - t) * (-1/v^2) in v, take v^{-1} coeff,
    then negate (we return S = -res_inf)."""
    v = sp.Symbol('v')
    integ = sp.together(F_expr/(sp.expand(c_expr) - t))
    gv = integ.subs(y, 1/v) * (-1)/v**2
    gv = sp.cancel(sp.together(gv))
    num, den = sp.fraction(gv)
    num = sp.Poly(sp.expand(num), v)
    den = sp.Poly(sp.expand(den), v)
    # Laurent expansion around v = 0: den = v^s * unit
    sden = min(e for (e,) in den.as_dict())
    unit = sp.Poly({(e - sden,): co for (e,), co in den.as_dict().items()}, v)
    # need num/den v^{-1} coefficient: series of num/unit up to order sden
    ser = sp.series(num.as_expr()/unit.as_expr(), v, 0, sden + 2).removeO()
    ser = sp.Poly(sp.expand(ser), v)
    coeff = ser.as_dict().get((sden - 1,), 0)
    return sp.cancel(-coeff)

# ---------------- the sheared (3,1,3) pair as base point ----------------
q = gam + x**2*y
c_pair = sp.expand(x*q*(1 + x**5*q**2))
w_pair = sp.expand(-lam*(1/x + sp.Rational(3, 2)*x**4*q**2))
gsub = [(gam, sp.Integer(3)), (lam, sp.Integer(5))]
cP = sp.expand(c_pair.subs(gsub)); wP = sp.expand(w_pair.subs(gsub))
cy = sp.diff(cP, y)

print("== (I1) chain-rule trace identity, j = 1, 2, 3, on the pair ==")
Jxy = sp.expand(sp.diff(wP, x)*cy - sp.diff(wP, y)*sp.diff(cP, x))   # = 5x
for j in (1, 2, 3):
    lhs = sp.cancel(sp.diff(trace_of(wP**j, cP), x))
    # Tr(w^{j-1} J / c_y): rational in y -> clear c_y via multiplication trick:
    # F = w^{j-1} * Jxy / cy is rational; use res_inf identity instead:
    # Tr(F) = sum_i F(y_i) with F = G/cy, G polynomial: Tr(G/cy) = sum res of
    # G dy/(c - t) at finite poles ... = -res_inf(G dy/(c-t))  (simple roots)
    rhs = sp.cancel(j * res_inf(wP**(j-1)*Jxy, cP))
    check(f"(I1) j={j}: d/dx Tr(w^j) = j*Tr(w^(j-1) J/c_y)",
          sp.simplify(lhs - rhs) == 0)

print("\n== (I1general) on a random non-pair (w,c) ==")
crnd = sp.expand(2*x + x**3*y + 7*x**4*y**2 + x**9*y**3)
wrnd = sp.expand(3/x + x**2*y - 2*x*y**2)
Jr = sp.expand(sp.diff(wrnd, x)*sp.diff(crnd, y) - sp.diff(wrnd, y)*sp.diff(crnd, x))
for j in (1, 2):
    lhs = sp.cancel(sp.diff(trace_of(wrnd**j, crnd), x))
    rhs = sp.cancel(j * res_inf(wrnd**(j-1)*Jr, crnd))
    check(f"(I1g) j={j} random (w,c)", sp.simplify(lhs - rhs) == 0)

print("\n== (I2) S_j = -res_inf(w^j dy/(c-t)) vs direct trace, pair ==")
for j in (0, 1, 2):
    S_res = res_inf(wP**j, cP)
    # direct: Tr(w^j/c_y) = -res_inf by the same residue theorem; instead
    # verify NUMERICALLY at sample (x, t): sum over numeric roots
    import random
    ok = True
    for (xv, tv) in [(sp.Rational(1, 3), 2), (sp.Rational(1, 5), 7)]:
        pc = sp.Poly(sp.expand(cP.subs(x, xv)) - tv, y)
        roots = pc.nroots(n=35)
        direct = sum((wP.subs([(x, xv), (y, r)]))**j /
                     (cy.subs([(x, xv), (y, r)])) for r in roots)
        pred = sp.N(S_res.subs([(x, xv), (t, tv)]), 30)
        ok &= abs(complex(sp.N(direct - pred, 25))) < 1e-18
    check(f"(I2) j={j} residue formula matches numeric trace", ok,
          f"S_{j} = {sp.simplify(S_res)}")

print("\n== S_1, S_2 profiles for the sheared pair (base point data) ==")
for j in (1, 2):
    Sj = sp.simplify(res_inf(wP**j, cP))
    print(f"  S_{j} =", Sj)

print("\n== (I3) sanity: D = 2 trace consequence reproduced ==")
# For the pair (m = 1): trace of 2 p2(c) w + p1(c) = mu/x would give
# 2 p2(t) tau(t) + 3 p1(t) = 3 mu/x with tau = Tr(w) = 0 -> mu = 0. Confirm
# Tr(w) = 0 via the machinery here:
check("Tr(w) = 0 on the pair", sp.simplify(trace_of(wP, cP)) == 0)

fails = [n for n, ok in results if not ok]
print(f"\n{len(results)} checks, {len(results)-len(fails)} passed, {len(fails)} failed")
sys.exit(1 if fails else 0)
