#!/usr/bin/env python3
# RETIRED round-11 (audit E6, blocker 5): stale standalone runner — its direct
# extraction evaluates to nan and fails Ghat0 == Z**2 - t. The maintained peel
# faces are cluster_peel_anchoredness.py, codexchk_round10b_peel_rig.py, and
# codexchk_round10c_normalization.py (all green). Kept verbatim below for the
# audit trail; NOT part of the green gate. Expected-red by design.
"""Standalone peel-block congruence rigs (mixed example) — 3rd attempt.
Stage 1: rebuild Hensel+peel, cache G_anch and V = G_deep*H.
Stage 2: cheap congruence rigs mod G_anch (all deg <= 2 operations).
Unpiped output; heavy flush. If externally killed again, the rig is
recorded as deferred to round-10b (its brief soft-spot (e))."""
import sympy as sp

x, y, Y, W, Z, t, lam = sp.symbols("x y Y W Z t lam")
TR = 10  # slimmer truncation than the fast pass; enough for levels <= 4

def xtrunc(e, n=TR):
    e = sp.expand(e)
    return sum(term for term in sp.Add.make_args(e)
               if term.as_powers_dict().get(x, 0) < n)

def hensel(F, var, A0, B0, levels):
    dA = sp.degree(A0, var)
    degF = sp.degree(sp.Poly(F, var))
    A, B = sp.expand(A0), sp.expand(B0)
    for ell in range(1, levels):
        R = xtrunc(F - sp.expand(A*B))
        lvl = [sp.Integer(0)]*(degF+1)
        for term in sp.Add.make_args(sp.expand(R)):
            pw = term.as_powers_dict()
            if pw.get(x, 0) == ell:
                yd = int(pw.get(var, 0))
                lvl[yd] += term/x**ell/var**yd
        if all(v == 0 for v in lvl):
            continue
        aco = [sp.Symbol(f"qa{ell}_{j}") for j in range(dA)]
        bco = [sp.Symbol(f"qb{ell}_{j}") for j in range(degF-dA+1)]
        dAe = sum(aco[j]*var**j for j in range(dA))
        dBe = sum(bco[j]*var**j for j in range(degF-dA+1))
        expr = sp.expand(dAe*B.subs(x,0) + dBe*A.subs(x,0)
                         - sum(sp.expand(v)*var**j for j, v in enumerate(lvl)))
        eqs = [sp.Poly(expr, var).coeff_monomial(var**j) for j in range(degF+1)]
        sol = list(sp.linsolve(eqs, aco+bco))[0]
        sub = dict(zip(aco+bco, sol))
        A = sp.expand(A + x**ell*dAe.subs(sub))
        B = sp.expand(B + x**ell*dBe.subs(sub))
    return A, B

print("stage 1: Hensel + peel ...", flush=True)
z_expr = x**7*y
c_mixed = sp.expand(x + x**-8*z_expr**2*(z_expr-1)**3
                    + x**-6*z_expr**2*(z_expr-1)**2
                    + 2*x**-4*z_expr*(z_expr-1)**2)
ctil = sp.expand(c_mixed.subs(y, x**-7*(1+Y)))
F = sp.expand(x**8*(ctil - t))
G, H = hensel(F, Y, Y**3, (1+Y)**2, TR)
Gstar = sp.expand(G.subs(Y, x**2*W)/x**6)
lin = sp.expand(sp.simplify(sp.expand(Gstar.subs(x, 0))/W**2))
Astar, Bstar = hensel(Gstar, W, W**2, lin, TR-2)
G_anch = sp.expand(x**4*Astar.subs(W, Y/x**2))
G_deep = sp.expand(x**2*Bstar.subs(W, Y/x**2))
V = sp.expand(G_deep*H)
assert xtrunc(sp.expand(F - G_anch*V), TR-4) == 0
Ghat0 = sp.expand((G_anch.subs(Y, x**3*Z)/x**6).subs(x, 0))
assert sp.expand(Ghat0 - (Z**2 - t)) == 0
print("stage 1 done: G_anch edge = Z^2 - t; F = G_anch*V ok", flush=True)

def red2(f):
    p = sp.Poly(f, Y)
    while p.degree() >= 2:
        d = p.degree(); lc = p.coeff_monomial(Y**d)
        f = sp.expand(f - lc*Y**(d-2)*G_anch)
        p = sp.Poly(f, Y)
    return f

Vbar = red2(sp.expand(V))
print("stage 2: rigs ...", flush=True)

def rig(mp, tdeg, lam_val, window):
    unk = []; wv = sp.Integer(0)
    for ell in range(window[0], window[1]+1):
        a = sum(sp.Symbol(f"ra{ell}_{j}")*t**j for j in range(tdeg+1))
        b = sum(sp.Symbol(f"rb{ell}_{j}")*t**j for j in range(tdeg+1))
        unk += [sp.Symbol(f"ra{ell}_{j}") for j in range(tdeg+1)]
        unk += [sp.Symbol(f"rb{ell}_{j}") for j in range(tdeg+1)]
        wv += x**ell*(a + b*x**-3*Y)
    L = lam if lam_val is None else lam_val
    Jg = sp.expand(sp.diff(wv, x)*sp.diff(G_anch, Y)
                   - sp.diff(wv, Y)*sp.diff(G_anch, x))
    T = red2(sp.expand(Vbar*red2(Jg)))
    T = sp.expand(T - L*x**(mp+8))
    T = sp.expand(red2(T).subs(Y, x**3*Z))
    eqs = []
    for lv in range(-4, 6):
        c0 = sp.Integer(0); c1 = sp.Integer(0)
        for term in sp.Add.make_args(T):
            pw = term.as_powers_dict()
            if pw.get(x, 0) == lv:
                zd = pw.get(Z, 0)
                piece = term/x**lv/Z**zd
                if zd == 0: c0 += piece
                else: c1 += piece
        for cc in (c0, c1):
            cc = sp.expand(cc)
            if cc != 0:
                for co in sp.Poly(cc, t).all_coeffs():
                    eqs.append(sp.expand(co))
    return eqs, unk

eqs, unk = rig(-6, 2, 1, (-3, 1))
r1 = sp.linsolve(eqs, unk) == sp.EmptySet
print(f"RIG m'=-6 (M-hat=-2): pair-kill infeasible: {r1}", flush=True)
eqs0, unk0 = rig(-4, 2, 0, (-3, -1))
s0 = sp.linsolve(eqs0, unk0)
r2 = s0 != sp.EmptySet and all(v == 0 for v in list(s0)[0])
print(f"RIG m'=-4 polar window: forced zero: {r2}", flush=True)
print("PEEL RIGS: " + ("ALL PASS" if (r1 and r2) else "CHECK OUTPUT"), flush=True)
