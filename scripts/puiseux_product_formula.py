#!/usr/bin/env python3
"""Puiseux ladder product formula — machine face (jc-4, 2026-07-25).

Round-9 deep-fix A5. Replaces lem:puiseux's Taylor line with the exact
identity  c - t = a_N * prod_i (y - y_i), so that at any truncation T:

    v(T) := ord_x( c(x,T) - t ) = ord(a_N) + sum_i ord(y_i - T)   (PF)

— exact, no Taylor remainder anywhere (the Finding-9 O(dev^2) demand is
vacated, not repaired). Per-branch orders are computed as written;
case rules (gap-dominated vs tail-dominated) are proof bookkeeping,
not part of the formula.

Terminal-cell corollary faced here on synthetics:
  v_term - sig1 = ord(y_p - T)  whenever every other branch's gap to
  y_p is strictly deeper than y_p's own tail (true at terminal cells:
  gaps >= sig1 > tail scale), giving v_term = sig1 - sig_tail >= 1
  and deriving the observed +1 constant.

Higher terminal rows: [Z^b] row of the terminal cell has
  ord = ord( c^{(b)}(x,T)/b! ) - b*sig1,  c^{(b)}/b! = a_N * e_{N-b}
(elementary symmetric sums of (T - y_i) with b factors dropped):
  >= 0 for all b >= 1, EXACTLY 0 for b = 1 -> c~(0,Z) = aZ + cone rows.

Checks (all exact, self-contained):
  P1  (PF) itself on a 3-branch synthetic with a fine split ABOVE the
      dev scale and a shallow-tail branch BELOW it — the two dominance
      regimes — direct expansion vs formula.
  P2  Terminal bracket theorem: v - sig1_q = ord(y_p - T) when gaps
      all deeper than tail; +1 constant reproduced.
  P3  Higher rows at a synthetic terminal cell: b=1 row order exactly
      matches the demand-zero anchoring; b>=2 rows order >= 0 with
      equality exactly at cone rows (co-cluster at gap scale sig1).
  P4  Ramified/covering covariance: substitute x = chi^2 and rerun P2
      — every valuation doubles, conclusions unchanged (ladder ratios
      are covering-covariant).
"""

import sympy as sp

x, y, t = sp.symbols("x y t")

def check(name, cond):
    if not cond:
        raise AssertionError(name)
    print(f"PASS  {name}")

def xval(e):
    e = sp.expand(e)
    if e == 0:
        return None
    return min(term.as_powers_dict().get(x, 0) for term in sp.Add.make_args(e))

print("== P1. The product formula, both dominance regimes ==")
# tracked y1; sibling y2 splits at scale 3 = dev scale; y3 has shallow tail
y1 = x**-6 + x**-4 + x**-3          # tracked branch (tail after T: x^-3)
y2 = x**-6 + x**-4 - 2*x**-3        # splits from y1 at scale 3
y3 = x**-6 + x**-4 - 5*x**-2        # agrees with T fully; own tail at 2
aN = x**18
T  = x**-6 + x**-4                  # truncation; dev = y1 - T ~ x^-3
c_syn = sp.expand(aN*(y - y1)*(y - y2)*(y - y3)) + t

v_direct = xval(sp.expand(c_syn.subs(y, T) - t))
v_formula = 18 + xval(sp.expand(y1 - T)) + xval(sp.expand(y2 - T)) + xval(sp.expand(y3 - T))
check(f"(PF) exact: direct {v_direct} = sum-of-orders {v_formula}", v_direct == v_formula)
check("regime check: ord(y2-T) = -3 (gap-dominated), ord(y3-T) = -2 (tail-dominated)",
      xval(sp.expand(y2 - T)) == -3 and xval(sp.expand(y3 - T)) == -2)

print("\n== P2. Terminal bracket theorem and the +1 constant ==")
# terminal-style cluster: tracked branch SIMPLE, all gaps at scale >= sig1,
# tail shallower than sig1. sig1 = 3 here; gaps at 3 and 4; tail at 2.
p  = x**-7 + x**-3 + 9*x**-2        # tracked; tail scale 2 < sig1 = 3
q1 = x**-7 - x**-3                  # gap scale 3 (= sig1: the cone partner)
q2 = x**-7 + 4*x**-4                # gap scale 4 (deeper split)
aN2 = x**15
c_term = sp.expand(aN2*(y - p)*(y - q1)*(y - q2)) + t
T2 = x**-7 + x**-3                  # the chain's full polar truncation
# value law surrogate: ord(c_y(p)):
cy = sp.diff(c_term, y)
sig1_q = xval(sp.expand(cy.subs(y, p)))
v_t = xval(sp.expand(c_term.subs(y, T2) - t))
check(f"terminal bracket: v - ord(c_y|p) = ord(p - T) "
      f"({v_t} - {sig1_q} = {xval(sp.expand(p - T2))})",
      v_t - sig1_q == xval(sp.expand(p - T2)))
check("the +1-class constant: v_term = sig1 - sig_tail = "
      f"{sig1_q} - 2 = {v_t} >= 1", v_t == sig1_q - 2 and v_t >= 1)

print("\n== P3. Higher terminal rows on DEMAND-CONSISTENT (exact-hit) cells ==")
# Exact-hit = chart slope equals the demand sig1 := ord(c_y|p)/q — I1's
# terminal property, imposed here by choosing ord(a_N). First a
# NEGATIVE CONTROL: the P2 cell charted at a NON-demand slope has
# [Z^1] order != 0 — exact-hit is a genuine hypothesis, not automatic:
Zs = sp.Symbol("Z")
cell_bad = sp.expand(c_term.subs(y, T2 + x**-3*Zs) - t)
rows_bad = sp.Poly(cell_bad, Zs).all_coeffs()[::-1]
ords_bad = [xval(r) for r in rows_bad]
print("   NEGATIVE CONTROL (chart slope 3 != demand 8): rows", ords_bad)
check("mismatched frame: height-1 row NOT at 0 (exact-hit fails)", ords_bad[1] != 0)

# (a) two-branch exact-hit cell, no cone partner:
#     p tail scale 1 < sig1 = 2; gap 3 > sig1; ord(a_N)=5 => demand 2.
p3a = x**-5 + 3*x**-1
q3a = x**-5 - x**-3
c3a = sp.expand(x**5*(y - p3a)*(y - q3a)) + t
sig1a = xval(sp.expand(sp.diff(c3a, y).subs(y, p3a)))
check("demand = 2 by construction", sig1a == 2)
cell_a = sp.expand(c3a.subs(y, x**-5 + x**-2*Zs) - t)
rows_a = [xval(r) for r in sp.Poly(cell_a, Zs).all_coeffs()[::-1]]
print("   exact-hit cell (no cone partner): rows", rows_a)
check("[Z^0] = sig1 - sig_tail = 2 - 1 = 1 (the +1 constant, derived)",
      rows_a[0] == 1)
check("[Z^1] = 0 exactly (aZ anchoring)", rows_a[1] == 0)
check("[Z^2] > 0 (no cone partner: gap 3 > sig1)", rows_a[2] is None or rows_a[2] > 0)

# (b) three-branch exact-hit cell WITH a cone partner at gap = sig1:
p3b = x**-5 + 3*x**-1
q3b1 = x**-5 - x**-3
q3b2 = x**-5 + 2*x**-2
c3b = sp.expand(x**7*(y - p3b)*(y - q3b1)*(y - q3b2)) + t
sig1b = xval(sp.expand(sp.diff(c3b, y).subs(y, p3b)))
check("demand = 2 by construction (three branches)", sig1b == 2)
cell_b = sp.expand(c3b.subs(y, x**-5 + x**-2*Zs) - t)
rows_b = [xval(r) for r in sp.Poly(cell_b, Zs).all_coeffs()[::-1]]
print("   exact-hit cell (cone partner at gap = sig1): rows", rows_b)
check("[Z^0] = 1, [Z^1] = 0, [Z^2] = 0 (cone row), [Z^3] >= 0",
      rows_b[0] == 1 and rows_b[1] == 0 and rows_b[2] == 0
      and (rows_b[3] is None or rows_b[3] >= 0))

print("\n== P4. Covering covariance (on the exact-hit cell) ==")
chi = sp.Symbol("chi")
c_cov = sp.expand(c3a.subs(x, chi**2))
def chival(e):
    e = sp.expand(e)
    if e == 0:
        return None
    return min(term.as_powers_dict().get(chi, 0) for term in sp.Add.make_args(e))
cell_cov = sp.expand(c_cov.subs(y, chi**-10 + chi**-4*Zs) - t)
rows_cov = [chival(r) for r in sp.Poly(cell_cov, Zs).all_coeffs()[::-1]]
print("   covered cell rows:", rows_cov)
check("rows double: [2,0,>0] (ratios/conclusions covering-covariant)",
      rows_cov[0] == 2 and rows_cov[1] == 0
      and (rows_cov[2] is None or rows_cov[2] > 0))

print("\nALL PRODUCT-FORMULA CHECKS PASSED")
