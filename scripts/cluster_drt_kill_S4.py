#!/usr/bin/env python3
"""Cluster-DRT S4 kill core, lean staged run (tdeg=3, window -2..1).
Same mathematics as cluster_drt_kill.py S4/S5; smaller corroboration
window so the run fits interactive time. The ANY-degree claim is the
lemma's (discriminant denominator); this is its machine corroboration.
"""
import sympy as sp

x, Y, Z, t, lam = sp.symbols("x Y Z t lam")
TRUNC = 12

def check(name, cond):
    if not cond:
        raise AssertionError(name)
    print(f"PASS  {name}", flush=True)

def xtrunc(e, n=TRUNC):
    e = sp.expand(e)
    return sum(term for term in sp.Add.make_args(e)
               if term.as_powers_dict().get(x, 0) < n)

F = Y**4 + 2*Y**3 + (1 + x**3)*Y**2 + x**3*Y + x**7 - t*x**6
g1, g0 = sp.Integer(0), sp.Integer(0)
h1, h0 = sp.Integer(2), sp.Integer(1)
for ell in range(1, TRUNC):
    R = xtrunc(F - sp.expand((Y**2 + g1*Y + g0)*(Y**2 + h1*Y + h0)))
    r = [sp.Integer(0)]*4
    for term in sp.Add.make_args(sp.expand(R)):
        pw = term.as_powers_dict()
        if pw.get(x, 0) == ell and pw.get(Y, 0) <= 3:
            r[pw.get(Y, 0)] += term / x**ell / Y**pw.get(Y, 0)
    r = [sp.expand(v) for v in r]
    dg0 = r[0]; dg1 = sp.expand(r[1] - 2*r[0])
    dh1 = sp.expand(r[3] - dg1); dh0 = sp.expand(r[2] - dg0 - 2*dg1)
    g0 += dg0*x**ell; g1 += dg1*x**ell; h0 += dh0*x**ell; h1 += dh1*x**ell

G = Y**2 + g1*Y + g0
ctil = t + sp.expand(F/x**6)
print("Hensel done", flush=True)

def reduce_mod(f, Gm):
    f = sp.expand(f)
    p = sp.Poly(f, Y)
    while p.degree() >= 2:
        d = p.degree(); lc = p.coeff_monomial(Y**d)
        f = sp.expand(f - lc*Y**(d-2)*Gm)
        p = sp.Poly(f, Y)
    return f

def jac(a, b):
    return sp.expand(sp.diff(a, x)*sp.diff(b, Y) - sp.diff(a, Y)*sp.diff(b, x))

def pair_system(mp, tdeg, lam_val, window):
    unk = []
    w = sp.Integer(0)
    for ell in range(window[0], window[1]+1):
        a = sum(sp.Symbol(f"a_{ell}_{j}")*t**j for j in range(tdeg+1))
        b = sum(sp.Symbol(f"b_{ell}_{j}")*t**j for j in range(tdeg+1))
        unk += [sp.Symbol(f"a_{ell}_{j}") for j in range(tdeg+1)]
        unk += [sp.Symbol(f"b_{ell}_{j}") for j in range(tdeg+1)]
        w += x**ell*(a + b*x**-3*Y)
    L = lam if lam_val is None else lam_val
    T = reduce_mod(sp.expand(jac(w, ctil) - L*x**mp), G)
    T = sp.expand(T.subs(Y, x**3*Z))
    eqs = []
    for lv in range(-9, TRUNC-8):
        c0 = sp.Integer(0); c1 = sp.Integer(0)
        for term in sp.Add.make_args(T):
            pw = term.as_powers_dict()
            if pw.get(x, 0) == lv:
                zd = pw.get(Z, 0)
                piece = term/x**lv/Z**zd
                if zd == 0: c0 += piece
                elif zd == 1: c1 += piece
        for c in (c0, c1):
            c = sp.expand(c)
            if c != 0:
                for co in sp.Poly(c, t).all_coeffs():
                    eqs.append(sp.expand(co))
    return eqs, unk

print("S4a: m'=-5, lam=1, t-poly tdeg=3 ...", flush=True)
eqs, unk = pair_system(-5, 3, 1, (-2, 1))
sol = sp.linsolve(eqs, unk)
check("lambda=1, t-polynomial: NO solution (THE KILL)", sol == sp.EmptySet)

print("S4b: m'=-5, lam free ...", flush=True)
eqs, unk = pair_system(-5, 3, None, (-2, 1))
sol = sp.linsolve(eqs, unk + [lam])
check("lam-free: solutions exist", sol != sp.EmptySet)
check("every t-polynomial solution has lambda = 0",
      all(sp.simplify(s[-1]) == 0 for s in sol))

print("S5: deeper-polar-only window, lam=0 ...", flush=True)
eqs0, unk0 = pair_system(-5, 3, 0, (-4, -2))
sol0 = sp.linsolve(eqs0, unk0)
check("polar window solvable", sol0 != sp.EmptySet)
check("all deeper-polar coefficients forced to ZERO",
      all(v == 0 for v in list(sol0)[0]))

print("S6: control m'=-3 (M-hat=+1) ...", flush=True)
eqs2, unk2 = pair_system(-3, 3, 1, (-2, 1))
sol2 = sp.linsolve(eqs2, unk2)
print(f"   m'=-3 pair-level system solvable: {sol2 != sp.EmptySet}", flush=True)

print("Q(t) denominator exhibit:", flush=True)
u0, u1, c_sym = sp.symbols("u0 u1 c_sym")
Etil = Z**2 + Z - t
expr = sp.expand(-(u0 + u1*Z)*(2*Z + 1) - c_sym)
remv = sp.rem(expr, Etil, Z)
solqt = sp.solve([sp.Poly(remv, Z).coeff_monomial(Z**d) for d in (0, 1)],
                 [u0, u1], dict=True)
d0 = sp.denom(sp.together(solqt[0][u0])); d1 = sp.denom(sp.together(solqt[0][u1]))
check("Q(t)-solution carries the discriminant denominator (1+4t)",
      sp.rem(sp.Poly(d0, t), sp.Poly(1+4*t, t)) == 0 or
      sp.rem(sp.Poly(d1, t), sp.Poly(1+4*t, t)) == 0)
crit = sp.expand(expr.subs(Z, sp.Rational(-1, 2)))
check("critical-point evaluation: LHS(Z0) = -c (any rho) => c forced 0",
      crit == -c_sym)

print("\nS4 LEAN CORE: ALL PASSED", flush=True)
