#!/usr/bin/env python3
"""T4 LEDGER SUM machine face — EXACT.
Product over the N sheets of c_y  =  Res_y(c - t, c_y) / a_N(x)^{N-1},
so  Sum_p ord_p(c_y)  =  val_x Res_y(c - t, c_y)  -  (N-1) val_x a_N,
where val_x = lowest x-power with t-generic nonzero coefficient.
Claim under test (T4 ledger): the sum equals N * sigma1."""
import sympy as sp
x, y, t = sp.symbols('x y t')

def ledger_sum(c, N, aN_val):
    res = sp.expand(sp.resultant(sp.expand(c - t), sp.diff(c, y), y))
    poly = sp.Poly(res, x)
    val = None
    for (e,), coeff in sorted(poly.terms()):
        if sp.simplify(coeff) != 0:      # generic t: nonzero as poly in t
            val = e; break
    return val - (N - 1) * aN_val

tau = sp.Rational(1, 3)
c1 = (x - x**4*y + x**14*y**2 + tau*x**11*y**2 - 2*tau*x**21*y**3
      + tau*x**31*y**4 - x**8*y**2 + 3*x**18*y**3 - 3*x**28*y**4 + x**38*y**5)
s1 = ledger_sum(c1, 5, 38)
print(f"Part-C  (N=5, sigma1=4): ledger sum = {s1}, predict 20")
assert s1 == 20

c2 = x - x**4*y + x**11*y**2 + x**8*y**2 - 2*x**15*y**3 + x**22*y**4
s2 = ledger_sum(c2, 4, 22)
print(f"part-B  (N=4, sigma1=4): ledger sum = {s2}, predict 16")
assert s2 == 16

c3 = x + x**2*y**2
s3 = ledger_sum(c3, 2, 2)
print(f"control (N=2, sigma1=1): ledger sum = {s3}, predict 2")
assert s3 == 2

# staircase-ish 3-edge witness: c = x + x^3 y + x^7 y^2 + x^13 y^3 (slopes 2,4,6 from anchor: sigma1 = 2)
c4 = x + x**3*y + x**7*y**2 + x**13*y**3
s4 = ledger_sum(c4, 3, 13)
# TRUE hull: (0,0)->(3,1)->(7,2)->(13,3), slopes 3,4,6; sigma1 = 3.
# NAIVE claim N*sigma1 = 9 is FALSE here: generic (untuned) c sits at the
# per-branch FLOORS sigma_s + I_s = 3, 3, 1 -> sum 7. The gap 9 - 7 = 2 is
# EXACTLY Lemma U's telescoped deficit (N2-1)(sigma3-sigma2) = 1*2.
# => the ledger identity is Sum = N*sigma1 - deficit + cancellation-raises;
# Part-C matched 20 only BECAUSE its exact-hit tuning raises deep branches.
# This kills the naive Step-3/4 ledger input and replaces it with the
# UNCONDITIONAL local fact: first-edge branches sit at value sigma1 exactly
# (zero deficit, raises only for repeated edge-1 roots) >= 1 -> Keller
# differential polar there -> (A) holds with no global sum at all.
print(f"3-edge  (N=3, sigma1=3): ledger sum = {s4} = N*sigma1 - U_deficit = 9 - 2  [naive law FALSIFIED, U confirmed]")
assert s4 == 7
print("LEDGER CHECKS DONE: naive sum law falsified by c4; U-deficit accounting exact; first-edge-local route stands")
