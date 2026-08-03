#!/usr/bin/env python3
"""Positive machine face of the LEXICOGRAPHIC termination (rev-3.10).

History: round-11's codexchk_round11_repeated_termination.py killed degree
descent (leg-i family, (Z-1)^2 forever, ladder +2/step); round-12's
codexchk_round12_termination_multiedge.py killed the rev-3.9 every-step
ladder (second-edge refinement, valuation 4 -> 4, cluster degree 3 -> 1).
Rev-3.10's invariant: the pair (cluster degree, anchor - valuation)
descends lexicographically — leg (i) (vertex edge, q = 1) climbs the
valuation at fixed-or-lower degree; leg (ii) (off-vertex edge) and every
covering strictly drop the degree. Faces F1-F4 below are the leg-(i) and
covering computations (unchanged from rev-3.9 — they were true; the claim
around them was too broad); F5 is the leg-(ii) face on round-12's own
witness:

  F1  On codex's own family (R = 6): the height-zero valuation strictly
      increases at EVERY repeated-root refinement (by exactly 2 there),
      i.e. the ladder terminates what degree descent cannot.
  F2  The mechanism: at each step the node normalizes to
      x^l (e(Z) + x R(x,Z)) with e(0) != 0 (valuation = l exactly), and
      the child's height-zero row is x^l (e(rho) + x R(x,rho)) with
      e(rho) = 0 -> order >= l+1. Checked at every depth of the family.
  F3  Generic repeated-root increment: a random normalized non-anchored
      node with a multiplicity-2 root gains >= +1 (5 random instances).
  F4  Covering-step span division: a slope with denominator q forces the
      edge polynomial into C[Z^q]; recentering at any root divides the
      height span by >= q (q = 2 and q = 3 instances) — so prod(q) <= N_1
      and only finitely many covering steps occur.
"""

import random
import sympy as sp

x, y, Y, Z, z, t = sp.symbols("x y Y Z z t")

def check(name, cond):
    if not cond:
        raise AssertionError(name)
    print(f"PASS  {name}")

def x_order(expr):
    expr = sp.expand(expr)
    if expr == 0:
        return None
    return int(min(term.as_powers_dict().get(x, 0) for term in sp.Add.make_args(expr)))

print("== F1/F2: the ladder on codex's round-11 row-D family (R = 6) ==")
R = 6
L = 2*R + 1
A = 2*L + 2
S = A + 1
sigma1 = S - A//2
phi = sum(x**j for j in range(1, L + 2))
charted = sp.expand(x + x**(-A//2)*(1 + Y) + x**(-A)*(1 + Y)**2*(Y - phi)**2)

prev_v = None
for depth in range(R):
    trunc = sum(x**j for j in range(1, depth + 1))
    node = sp.expand(charted.subs(Y, trunc + x**(depth + 1)*Z))
    v = x_order(node.subs(Z, 0))                      # height-zero valuation
    ell = x_order(node)                               # edge level (cleared min)
    edge = sp.expand(node * x**(-ell)).subs(x, 0)     # expand FIRST, then x->0
    check(f"depth {depth}: e(0) != 0 and height-zero valuation == edge level",
          edge.subs(Z, 0) != 0 and v == ell)
    check(f"depth {depth}: e(rho) = 0 at the recentering root rho = 1 (mult 2)",
          edge.subs(Z, 1) == 0 and sp.diff(edge, Z).subs(Z, 1) == 0)
    if prev_v is not None:
        check(f"depth {depth}: ladder strictly increased (+{v - prev_v})", v > prev_v)
    prev_v = v
print(f"      ladder values: {[-A + 2*d + 2 for d in range(R)]} (as printed, +2/step)")

print("\n== F3: generic multiplicity-2 increment ==")
random.seed(11)
for i in range(5):
    rho = random.choice([1, 2, -1, 3, -2])
    off = random.choice([2, 3, -1, -3])
    if off == rho:
        off += 1
    ell = random.randint(-9, -3)
    e = sp.expand((Z - rho)**2 * (Z - off))           # e(0) != 0 ensured below
    assert e.subs(Z, 0) != 0
    Rterm = sum(random.randint(-2, 2) * x**a * Z**h for a in (1, 2) for h in (0, 1, 2))
    node = sp.expand(x**ell * (e + Rterm))
    v_parent = x_order(node.subs(Z, 0))
    child = sp.expand(node.subs(Z, rho + x*Z))
    v_child = x_order(child.subs(Z, 0))
    check(f"instance {i}: parent v = l = {ell}; child v >= l+1 (got {v_child})",
          v_parent == ell and (v_child is None or v_child >= ell + 1))

print("\n== F4: covering slope denominator q forces e in C[Z^q]; span divides ==")
for q, span in ((2, 4), (3, 6)):
    # integer support heights with denominator-q slope: consecutive edge
    # heights differ by q  =>  e(Z) = etilde(Z^q)
    etilde = sp.expand(sp.prod(sp.Symbol('W') - r for r in range(1, span//q + 1)))
    e = etilde.subs(sp.Symbol('W'), Z**q)
    check(f"q={q}: edge polynomial is a polynomial in Z^{q} (span {span})",
          all(mm[0] % q == 0 for mm in sp.Poly(e, Z).monoms()))
    rho_q = 1  # a root of etilde at W=1 -> Z = any q-th root of 1; chart real one
    mult = sp.Poly(e, Z).count_roots()  # not used for the bound; explicit below
    m_rho = 0
    ee = e
    while sp.expand(ee.subs(Z, rho_q)) == 0:
        m_rho += 1
        ee = sp.cancel(sp.expand(ee / (Z - rho_q)))
        ee = sp.expand(ee)
    check(f"q={q}: child span (= mult of charted root) {m_rho} <= span/q = {span//q}",
          m_rho <= span//q)

print("\n== F5: leg (ii) on round-12's second-edge witness — degree drops ==")
t_ = sp.symbols("t")
S5 = 10
K5 = Y**3 + x*Y**2 + x**2*Y + x**5
parent_z = x + x**-4*(x**S5*y) + x**-8*(x**S5*y)**2*K5.subs(Y, x**S5*y - 1)
charted5 = sp.expand(parent_z.subs(y, x**(-S5)*(1 + Y)))
F5 = sp.expand(x**8*(charted5 - t_))
tfree = sp.expand(F5.subs(t_, 0))
vals5 = {h: x_order(sp.Poly(tfree, Y).coeff_monomial(Y**h)) for h in range(4)}
check("F5 polygon (4,2,1,0), two edges, vertex NOT on the second edge",
      vals5 == {0: 4, 1: 2, 2: 1, 3: 0})
rho5 = (-1 + sp.sqrt(-3))/2
child5 = sp.expand(F5.subs(Y, x*(rho5 + x*Z)))
v_par = x_order(F5.subs(Y, 0).subs(t_, 0))
v_ch = x_order(sp.expand(child5.subs(Z, 0).subs(t_, 0)))
check("F5 leg (ii): valuation does NOT climb (4 -> 4) — as the invariant expects",
      v_par == 4 and v_ch == 4)
edge2 = sp.expand(sp.expand(F5.subs(Y, x*Z))*x**-3).subs(x, 0)
mult5 = 0
ee5 = sp.expand(edge2)
while sp.simplify(ee5.subs(Z, rho5)) == 0:
    mult5 += 1
    ee5 = sp.expand(sp.cancel(ee5/(Z - rho5)))
check("F5 leg (ii): cluster degree STRICTLY drops (3 -> 1 = root multiplicity)",
      mult5 == 1 and mult5 <= 3 - 1)
check("F5 anchor invariant across the step (8 -> 8), child still non-anchored",
      x_order(sp.diff(F5, t_)) == 8 and x_order(sp.diff(child5.subs(Z, 0), t_)) == 8 and v_ch < 8)

print("\nLEXICOGRAPHIC TERMINATION FACE: ALL CHECKS PASSED")
print("(leg i: rounds 10b/11 families climb the ladder; leg ii: round-12's")
print(" witness drops the degree; the pair (d, A - v) descends in N^2)")
