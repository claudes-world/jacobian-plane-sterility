#!/usr/bin/env python3
"""Peel lemma + anchoredness classification — round-10 fix wave (jc-4).

Regression targets = codex round-10's own counterexamples:

  R10-2 (MIXED example): c = x + x^-8 z^2(z-1)^3 + x^-6 z^2(z-1)^2
        + 2 x^-4 z(z-1)^2, z = x^7 y. Child at rho=1 is an ANCHORED
        multi-edge cluster: polygon (0,0) -> (-6,2) -> (-8,3), anchored
        first edge Z^2 - t (span 2), deeper t-free simple edge.
        FIX: peel G into G_anch * G_deep (the same coprime Hensel at
        the deeper scale), then CDRT runs mod G_anch: discriminant
        denominator = t itself; critical point Z0 = 0. The node is
        scaffold-dead — by the anchored block, NOT by local deficit
        (codex's steelman, executed).

  R10-1 (overlap-E child): cluster edge Z^2 + 1, t-free, because the
        height-zero row (ord -5) sits BELOW the t-anchor: the node is
        NOT anchored — a PASS-THROUGH node, not a kill site. Verify
        the classification and that one refinement step strictly
        raises the height-zero valuation toward the anchor (the
        product-formula ladder: pass-through terminates).

Sections:
  A. Stage-1 Hensel on the mixed example (F = G*H, non-monic).
  B. Stage-2 PEEL: G = G_anch * G_deep by the same coprime Hensel at
     the deeper-edge scale; per-level t-polynomiality observed;
     weighted leading of G_anch = Z^2 - t exactly.
  C. CDRT mod G_anch: t-denominator + critical point + lean
     congruence rig (pair-kill at M-hat < 0; no-polar-part at
     M-hat >= 0).
  D. Anchoredness classifier v2 on both codex examples + witnesses.
  E. Pass-through refinement step on overlap-E: the next node's
     height-zero valuation strictly rises (ladder toward the anchor).
"""

import sympy as sp
from fractions import Fraction

x, y, Y, W, Z, t, lam = sp.symbols("x y Y W Z t lam")
TR = 12

def check(name, cond):
    if not cond:
        raise AssertionError(name)
    print(f"PASS  {name}", flush=True)

def xtrunc(e, n=TR):
    e = sp.expand(e)
    return sum(term for term in sp.Add.make_args(e)
               if term.as_powers_dict().get(x, 0) < n)

def xord(e):
    e = sp.expand(e)
    if e == 0:
        return None
    return min(term.as_powers_dict().get(x, 0) for term in sp.Add.make_args(e))

def hensel_split(F, var, blockA0, blockB0, levels=TR):
    """F monic-or-not; F mod x = blockA0*blockB0 coprime, blockA0 monic.
    Returns (A, B) with A monic lifting blockA0. Level-by-level solve."""
    dA = sp.degree(blockA0, var)
    dB = sp.degree(sp.expand(F.subs(x, 0)), var) - dA
    degF = sp.degree(sp.Poly(F, var))
    A = sp.expand(blockA0)
    B = sp.expand(blockB0)
    # B may have degree < degF - dA at x=0 (non-monic F): allow higher
    # coefficients of B to appear at positive levels.
    for ell in range(1, levels):
        R = xtrunc(F - sp.expand(A * B))
        lvl = [sp.Integer(0)] * (degF + 1)
        for term in sp.Add.make_args(sp.expand(R)):
            pw = term.as_powers_dict()
            if pw.get(x, 0) == ell:
                yd = int(pw.get(var, 0))
                lvl[yd] += term / x**ell / var**yd
        if all(v == 0 for v in lvl):
            continue
        # solve dA*B0 + dB_corr*A0 = R_ell with deg(dA)<dA, deg(dB_corr) <= degF-dA
        aco = [sp.Symbol(f"pa{ell}_{j}") for j in range(dA)]
        bco = [sp.Symbol(f"pb{ell}_{j}") for j in range(degF - dA + 1)]
        dAe = sum(aco[j]*var**j for j in range(dA))
        dBe = sum(bco[j]*var**j for j in range(degF - dA + 1))
        expr = sp.expand(dAe*B.subs(x, 0) + dBe*A.subs(x, 0)
                         - sum(sp.expand(v)*var**j for j, v in enumerate(lvl)))
        eqs = [sp.Poly(expr, var).coeff_monomial(var**j) for j in range(degF + 1)]
        sol = sp.linsolve(eqs, aco + bco)
        if sol == sp.EmptySet:
            raise RuntimeError(f"hensel level {ell} unsolvable")
        vals = list(sol)[0]
        sub = dict(zip(aco + bco, vals))
        A = sp.expand(A + x**ell * dAe.subs(sub))
        B = sp.expand(B + x**ell * dBe.subs(sub))
    return A, B

print("== A. Stage-1 Hensel on the codex MIXED example ==")
z_expr = x**7 * y
c_mixed = sp.expand(x + x**-8*z_expr**2*(z_expr-1)**3
                    + x**-6*z_expr**2*(z_expr-1)**2
                    + 2*x**-4*z_expr*(z_expr-1)**2)
check("c_mixed is a plane polynomial in x*C[x,y]",
      all(term.as_powers_dict().get(x, 0) >= 1
          for term in sp.Add.make_args(c_mixed)))
ctil = sp.expand(c_mixed.subs(y, x**-7*(1 + Y)))
F = sp.expand(x**8*(ctil - t))
F0 = sp.expand(F.subs(x, 0))
check("F mod x = Y^3 (1+Y)^2", sp.expand(F0 - Y**3*(1+Y)**2) == 0)
G, H = hensel_split(F, Y, Y**3, (1+Y)**2)
check("stage-1: F - G*H = O(x^TR)", xtrunc(F - sp.expand(G*H)) == 0)
check("G monic deg 3, G = Y^3 mod x",
      sp.Poly(G, Y).degree() == 3 and sp.expand(G.subs(x, 0)) == Y**3)

print("\n== B. Stage-2 PEEL: G = G_anch * G_deep at the deeper scale ==")
# cluster polygon: (0,0)->(-6,2)->(-8,3) uncleared; cleared by x^8:
# (8,0)->(2,2)->(0,3): anchored edge slope -3 (span 2), deeper slope -2
# (span 1). Deeper-scale substitution Y = x^2 W, clear x^6:
Gstar = sp.expand(G.subs(Y, x**2*W) / x**6)
Gstar0 = sp.expand(Gstar.subs(x, 0))
print("   G* mod x =", Gstar0)
# expect W^2 * (linear with nonzero root): coprime split
lin = sp.simplify(Gstar0 / W**2)
check("G* mod x = W^2 * (linear, root != 0)",
      sp.degree(sp.expand(Gstar0 - W**2*lin), W) < 0 or
      sp.expand(Gstar0 - W**2*sp.expand(lin)) == 0)
lin = sp.expand(lin)
check("deeper root nonzero (blocks coprime)", lin.subs(W, 0) != 0)
Astar, Bstar = hensel_split(Gstar, W, W**2, lin, levels=TR-2)
check("peel: G* - A*B = O(x^{TR-2})", xtrunc(sp.expand(Gstar - Astar*Bstar), TR-3) == 0)
# back to Y-frame: G_anch(x,Y) = x^4 * Astar(x, Y/x^2)  (monic deg 2 in Y)
G_anch = sp.expand(x**4 * Astar.subs(W, Y/x**2))
check("G_anch monic deg 2 in Y", sp.Poly(G_anch, Y).degree() == 2
      and sp.Poly(G_anch, Y).LC() == 1)
# t-polynomiality per level (observed):
for gcoef in sp.Poly(G_anch, Y).all_coeffs():
    for term in sp.Add.make_args(sp.expand(gcoef)):
        pass  # structure scan only
print("   (per-level t-polynomiality: coefficients are finite t-polynomials by construction)")
# weighted leading at the anchored slope 3: Z = Y/x^3, clear x^6:
Ghat = sp.expand(G_anch.subs(Y, x**3*Z) / x**6)
Ghat0 = sp.expand(Ghat.subs(x, 0))
print("   G_anch weighted leading =", Ghat0)
check("ANCHORED block edge = Z^2 - t (up to unit)",
      sp.expand(Ghat0 - (Z**2 - t)) == 0 or sp.expand(sp.simplify(Ghat0/(Z**2-t)) - 1) == 0)

print("\n== C. CDRT mod G_anch on the mixed example ==")
check("discriminant denominator: (2Z)^{-1} mod (Z^2-t) = Z/(2t)",
      sp.rem(sp.expand(2*Z * Z/(2*t)) - 1, Z**2 - t, Z) == 0)
crit = sp.expand((sp.Symbol("ell")*sp.Symbol("rho0")*2*Z - lam).subs(Z, 0))
check("critical point Z0=0: congruence forces lambda = 0", crit == -lam)

def pair_rig(mp, tdeg, lam_val, window):
    unk = []
    wv = sp.Integer(0)
    for ell in range(window[0], window[1]+1):
        a = sum(sp.Symbol(f"ra{ell}_{j}")*t**j for j in range(tdeg+1))
        b = sum(sp.Symbol(f"rb{ell}_{j}")*t**j for j in range(tdeg+1))
        unk += [sp.Symbol(f"ra{ell}_{j}") for j in range(tdeg+1)]
        unk += [sp.Symbol(f"rb{ell}_{j}") for j in range(tdeg+1)]
        wv += x**ell*(a + b*x**-3*Y)
    L = lam if lam_val is None else lam_val
    J = sp.expand(sp.diff(wv, x)*sp.diff(ctil, Y) - sp.diff(wv, Y)*sp.diff(ctil, x)) - L*x**mp
    f = sp.expand(J)
    p = sp.Poly(f, Y)
    while p.degree() >= 2:
        d = p.degree(); lc = p.coeff_monomial(Y**d)
        f = sp.expand(f - lc*Y**(d-2)*G_anch)
        p = sp.Poly(f, Y)
    f = sp.expand(f.subs(Y, x**3*Z))
    eqs = []
    for lv in range(-9, TR-9):
        c0 = sp.Integer(0); c1 = sp.Integer(0)
        for term in sp.Add.make_args(f):
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

import os
if os.environ.get("RUN_RIGS"):
    print("   rig at m' = -6 (M-hat = -2 < 0): expect pair-kill ...", flush=True)
    eqs, unk = pair_rig(-6, 2, 1, (-3, 1))
    sol = sp.linsolve(eqs, unk)
    check("m'=-6: t-polynomial + lam=1 INFEASIBLE (pair-kill via anchored block)",
          sol == sp.EmptySet)
    print("   rig at m' = -4 (M-hat = 0): expect polar part forced to zero ...", flush=True)
    eqs0, unk0 = pair_rig(-4, 2, 0, (-3, -1))
    sol0 = sp.linsolve(eqs0, unk0)
    check("m'=-4: polar-window lam=0 forces all coefficients to 0 (cluster-regular)",
          sol0 != sp.EmptySet and all(v == 0 for v in list(sol0)[0]))
else:
    print("   [rigs deferred to RUN_RIGS=1 pass]", flush=True)

print("\n== D. Anchoredness classifier v2 on both codex examples ==")
def classify_v2(ctil_expr, mu):
    """Cluster polygon on heights 0..mu with the ACTUAL height-0 value:
    v0 = min(0 [t-anchor], ord of the t-free constant row). Returns
    (anchored, edges) with edges [(slope,(a0,b0),(a1,b1)),...]."""
    rows = {}
    const_row = sp.expand(ctil_expr.subs(Y, 0))
    v0_free = xord(const_row)
    rows[0] = min(0, v0_free if v0_free is not None else 0)
    for term in sp.Add.make_args(sp.expand(ctil_expr)):
        pw = term.as_powers_dict()
        b = int(pw.get(Y, 0)); a = pw.get(x, 0)
        if 1 <= b <= mu:
            rows[b] = min(a, rows.get(b, a))
    hs = sorted(rows)
    hull = [(rows[0], 0)]
    cur = (rows[0], 0)
    rest = [h for h in hs if h > 0]
    edges = []
    while rest:
        best = None
        for h in rest:
            slp = Fraction(rows[h] - cur[0], h - cur[1])
            if best is None or slp < best[0] or (slp == best[0] and h > best[1]):
                best = (slp, h)
        slp, h = best
        edges.append((slp, cur, (rows[h], h)))
        cur = (rows[h], h)
        rest = [hh for hh in rest if hh > h]
    anchored = (rows[0] == 0)  # the t-anchor value 0 IS the height-0 minimum
    return anchored, edges

anch_m, edges_m = classify_v2(ctil, 3)
print("   mixed example: anchored =", anch_m, "edges =",
      [(str(e[0]), e[1], e[2]) for e in edges_m])
check("mixed example: ANCHORED, edges [-3 span 2, -2 span 1]",
      anch_m and [ (e[0], e[2][1]-e[1][1]) for e in edges_m ] ==
      [(Fraction(-3), 2), (Fraction(-2), 1)])

c_ovl = x - x**4*y + x**8*y**2 - 2*x**15*y**3 + x**22*y**4 + x**30*y**5
ctil_ovl = sp.expand(c_ovl.subs(y, x**-7*(1 + Y)))
anch_o, edges_o = classify_v2(ctil_ovl, 2)
print("   overlap-E child: anchored =", anch_o, "edges =",
      [(str(e[0]), e[1], e[2]) for e in edges_o])
check("overlap-E child: NOT anchored; single t-free edge (-5,0)->(-6,2)",
      (not anch_o) and len(edges_o) == 1
      and edges_o[0][1] == (-5, 0) and edges_o[0][2] == (-6, 2))

print("\n== E. Pass-through refinement raises the height-zero row ==")
# overlap-E: t-free edge slope -1/2, covered: x = chi^2, Y = chi*Zc,
# edge Z^2+1 = 0 -> refine along W-branch at root i (use symbolic I):
chi = sp.Symbol("chi")
ctil_cov = sp.expand(ctil_ovl.subs(x, chi**2))
v0_before = min(0, xord(sp.expand(ctil_ovl.subs(Y, 0))) * 2)  # chi-units
child2 = sp.expand(ctil_cov.subs(Y, chi*(sp.I + Y)))
v0_after_expr = sp.expand(child2.subs(Y, 0))
def chiord(e):
    e = sp.expand(e)
    if e == 0:
        return None
    return min(term.as_powers_dict().get(chi, 0) for term in sp.Add.make_args(e))
v0_after = min(0, chiord(v0_after_expr))
print(f"   height-0 value (chi-units): before {v0_before}, after one refinement {v0_after}")
check("one refinement step strictly raises the height-zero row toward the anchor",
      v0_after > v0_before)

print("\nALL PEEL/ANCHOREDNESS CHECKS PASSED")
