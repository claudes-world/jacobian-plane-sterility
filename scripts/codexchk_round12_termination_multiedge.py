#!/usr/bin/env python3
"""Round-12 counterexample to rev-3.9's every-step valuation ladder.

The repaired proof in thm:leafkill(1) only works when the refined edge is
the edge incident to the height-zero vertex.  A non-anchored cluster may
have several edges, and the printed classification routes roots on all of
them.  This checker builds an exact row-D plane parent and follows a simple
root on the *second* cluster edge.  The actual height-zero valuation does
not increase.

This does not construct an infinite route or a scaffold.  It refutes the
printed implication "pass-through => refined edge contains height zero",
and therefore refutes the claimed proof of termination.
"""

from __future__ import annotations

import sympy as sp


x, y, z, Y, W, Z, t = sp.symbols("x y z Y W Z t")


def check(name: str, condition: bool) -> None:
    if not condition:
        raise AssertionError(name)
    print(f"PASS  {name}")


def x_order(expr: sp.Expr) -> int:
    expr = sp.expand(expr)
    if expr == 0:
        raise AssertionError("zero has no finite x-order")
    return int(
        min(term.as_powers_dict().get(x, 0) for term in sp.Add.make_args(expr))
    )


print("== A. Exact row-D plane parent ==")

# z = x^S y.  The first edge has span two and supported height one.
# Its deeper edge has the repeated root z=1 of multiplicity three.
S = 10
K = Y**3 + x * Y**2 + x**2 * Y + x**5
parent_z = x + x**-4 * z + x**-8 * z**2 * K.subs(Y, z - 1)
parent = sp.expand(parent_z.subs(z, x**S * y))

check(
    "parent lies in x*C[x,y]",
    all(
        term.as_powers_dict().get(x, 0) >= 1
        and term.as_powers_dict().get(y, 0) >= 0
        for term in sp.Add.make_args(parent)
    ),
)

py = sp.Poly(parent, y)
rows = {
    h: x_order(py.coeff_monomial(y**h))
    for h in range(py.degree() + 1)
    if py.coeff_monomial(y**h) != 0
}
check(
    "row D: fiber first edge (0,0)-(2,12), height one supported",
    rows[0] == 1 and rows[1] == 6 and rows[2] == 12,
)
check(
    "the next edge has slope 10",
    rows[3] - rows[2] == rows[4] - rows[3] == rows[5] - rows[4] == 10,
)

deeper = sp.expand(x**8 * (parent_z - x - x**-4 * z)).subs(x, 0)
check("deeper edge polynomial is z^2(z-1)^3", sp.expand(deeper - z**2 * (z - 1) ** 3) == 0)
check(
    "z=1 is a canonical root of multiplicity exactly three",
    deeper.subs(z, 1) == 0
    and sp.diff(deeper, z, 1).subs(z, 1) == 0
    and sp.diff(deeper, z, 2).subs(z, 1) == 0
    and sp.diff(deeper, z, 3).subs(z, 1) != 0,
)


print("\n== B. The charted non-anchored degree-three cluster ==")

charted = sp.expand(parent_z.subs(z, 1 + Y))
F = sp.expand(x**8 * (charted - t))
expected_F = sp.expand(x**9 + x**4 * (1 + Y) + (1 + Y) ** 2 * K - x**8 * t)
check("charted fiber polynomial is exact", sp.expand(F - expected_F) == 0)
check(
    "special fiber has the governed zero cluster Y^3 times a unit",
    sp.expand(F.subs(x, 0) - Y**3 * (1 + Y) ** 2) == 0,
)

t_free = sp.expand(F.subs(t, 0))
vals = {
    h: x_order(sp.Poly(t_free, Y).coeff_monomial(Y**h))
    for h in range(4)
}
check(
    "cluster polygon values at heights 0..3 are (4,2,1,0)",
    vals == {0: 4, 1: 2, 2: 1, 3: 0},
)
check(
    "node is strictly non-anchored: free height-zero row 4 < t-anchor 8",
    vals[0] == 4 and x_order(sp.diff(F, t)) == 8,
)
check(
    "cluster has two edges, and the second edge does not contain height zero",
    (vals[1] - vals[0]) == -2
    and (vals[2] - vals[1]) == (vals[3] - vals[2]) == -1,
)


print("\n== C. Legal simple refinement on the second edge ==")

scaled = sp.expand(F.subs(Y, x * W))
second_edge = sp.expand(scaled / x**3).subs(x, 0)
check(
    "second edge polynomial is W(1+W+W^2), so e(0)=0",
    sp.expand(second_edge - W * (1 + W + W**2)) == 0
    and second_edge.subs(W, 0) == 0,
)

rho = (-1 + sp.sqrt(-3)) / 2
check(
    "rho is a nonzero simple root on that edge",
    rho != 0
    and sp.simplify(second_edge.subs(W, rho)) == 0
    and sp.simplify(sp.diff(second_edge, W).subs(W, rho)) != 0,
)

parent_height_zero = x_order(F.subs({Y: 0, t: 0}))
child = sp.expand(F.subs(Y, x * (rho + x * Z)))
child_height_zero = x_order(child.subs({Z: 0, t: 0}))
check(
    "actual height-zero valuation does NOT increase (4 -> 4)",
    parent_height_zero == child_height_zero == 4,
)
check(
    "the t-anchor remains strictly above the unchanged free row",
    x_order(sp.diff(child.subs(Z, 0), t)) == 8 and child_height_zero < 8,
)


print("\n== D. Covering claim: the missing monomial prefactor ==")

# On an edge beginning at height j0, denominator q only gives
# e(Z)=Z^j0 E(Z^q), not generally e in C[Z^q].
q_edge = Z**3 * (Z**2 - 1)
check(
    "q=2 edge support is one congruence class but e is not in C[Z^2]",
    {m[0] % 2 for m in sp.Poly(q_edge, Z).monoms()} == {1},
)
span = 2
zero_mult = 3
check(
    "at the printed 'ANY root' rho=0, multiplicity exceeds span/q",
    q_edge.subs(Z, 0) == 0 and zero_mult > span // 2,
)

print("\nROUND-12 R1 COUNTEREXAMPLE CONFIRMED")
print("The repair samples only first-edge refinements; multi-edge routing breaks")
print("the asserted +1-at-every-step ladder.  Termination may still admit a")
print("lexicographic repair (off-first-edge steps lower cluster degree), but it")
print("is not the proof frozen in rev-3.9.")
