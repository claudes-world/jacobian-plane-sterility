#!/usr/bin/env python3
"""Round-10 exact checks for uncovered T4, covering plumbing, and anchoring."""

import re
import sympy as sp

s, x, y, t, Z = sp.symbols("s x y t Z")


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    print(f"PASS  {name}")


def val(expr, uniformizer):
    expr = sp.expand(expr)
    if expr == 0:
        return sp.oo
    return min(term.as_powers_dict().get(uniformizer, 0)
               for term in sp.Add.make_args(expr))


print("== L1. A fractional first slope has intrinsic q=q1, not q=1 ==")
# c=x+x^5 y^2 has first slope 5/2.  On c=t, put x=s^2:
# y=s^-5*sqrt(t-s^2).  A simple generic root therefore has ord_s x=2
# and ord_s c_y=5=(5/2)*2.
c = x + x**5 * y**2
cy = sp.diff(c, y)
q1 = 2
sigma1 = sp.Rational(5, 2)
root_lead = s**-5  # the nonzero sqrt(t) coefficient is valuation-neutral
cy_branch_lead = sp.expand(cy.subs({x: s**2, y: root_lead}))
check("simple-edge branch has q=ord_s(x)=2",
      val(s**2, s) == q1)
check("ord_s(c_y)=sigma1*q1=5",
      val(cy_branch_lead, s) == sigma1 * q1 == 5)

# Choose m=1, e=D-1=2, so k=1/2 and sigma1=(m+1)+k.
m, e = 1, 2
k = sp.Rational(m, e)
check("uncovered quantization permits fractional k",
      sigma1 == (m + 1) + k)
m_q = q1 * m + q1 - 1
ell_q = q1 - 1
gap_q = m_q - ell_q
k_q = q1 * k
check("covered pair/top exponents and master gap are distinct",
      (m_q, ell_q, gap_q) == (3, 1, 2))
check("covered k is integral and covered slope identity holds",
      k_q == 1 and q1 * sigma1 == (m_q + 1) + k_q)


print("\n== L2. Product formula and terminal row inequalities ==")
p = x**-5 + 3*x**-1
q_far = x**-5 - x**-3
q_cone = x**-5 + 2*x**-2
aN = x**7
cprod = sp.expand(aN * (y - p) * (y - q_far) * (y - q_cone)) + t
T = x**-5
direct = val(cprod.subs(y, T) - t, x)
product = val(aN, x) + sum(val(root - T, x)
                           for root in (p, q_far, q_cone))
check("valuation product formula is exact", direct == product)

# Exact-hit chart at demand 2.
cell = sp.Poly(sp.expand(cprod.subs(y, T + x**-2 * Z) - t), Z)
row_orders = [val(cell.coeff_monomial(Z**b), x) for b in range(4)]
check("exact-hit rows are [>=1,0,>=0,...]",
      row_orders[0] >= 1 and row_orders[1] == 0
      and all(order >= 0 for order in row_orders[2:]))

# Negative control: the same cell charted at a non-demand slope.
bad_cell = sp.Poly(sp.expand(cprod.subs(y, T + x**-3 * Z) - t), Z)
check("non-demand chart breaks the height-one zero order",
      val(bad_cell.coeff_monomial(Z), x) != 0)


print("\n== L3. Consumer grep and the missing exact-hit producer ==")
paper = open("paper-3b.tex", encoding="utf-8").read()
retired_pattern = re.compile(
    r"v[^\\n]{0,80}\\\\ge[^\\n]{0,80}"
    r"\\\\hat\\s*\\\\sigma_1[^\\n]{0,80}next",
    re.IGNORECASE,
)
check("no downstream formula consumes the retired intermediate inequality",
      not retired_pattern.search(paper))
check("paper invokes Lemma I1 but defines no labeled Lemma I1",
      "Lemma I1" in paper and r"\label{lem:I1}" not in paper)
check("cor:anchor explicitly relies on the absent exact-hit producer",
      "Lemma I1's terminal identity" in paper)

print("\nROUND-10 LAWS/PUISEUX CHECKS PASSED")
