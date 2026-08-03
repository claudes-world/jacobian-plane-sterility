#!/usr/bin/env python3
"""Round-11 counterexample to the printed repeated-root degree descent.

For arbitrary R the constructed *row-D* plane parent has R consecutive
canonical non-anchored, t-free repeated-root refinements.  Every edge
polynomial is (Z-1)^2, so the cluster degree stays two rather than strictly
decreasing.  The height-zero valuation rises at each step, suggesting a
possible repaired termination argument, but that is not the degree-descent
proof printed in thm:T / thm:leafkill(1).
"""

from __future__ import annotations

import sympy as sp


x, y, Y, Z, z = sp.symbols("x y Y Z z")


def check(name: str, condition: bool) -> None:
    if not condition:
        raise AssertionError(name)
    print(f"PASS  {name}")


def x_order(expr: sp.Expr) -> int:
    expanded = sp.expand(expr)
    if expanded == 0:
        raise AssertionError("zero has no finite x-order in this checker")
    return int(
        min(term.as_powers_dict().get(x, 0) for term in sp.Add.make_args(expanded))
    )


def is_plane_x_multiple(expr: sp.Expr) -> bool:
    return all(
        term.as_powers_dict().get(x, 0) >= 1
        and term.as_powers_dict().get(y, 0) >= 0
        for term in sp.Add.make_args(sp.expand(expr))
    )


R = 12  # twice the prior six-step simple-root window
L = 2 * R + 1
A = 2 * L + 2
S = A + 1
sigma1 = S - A // 2
phi = sum(x**j for j in range(1, L + 2))

# In the z=x^S y frame:
#   x + x^(sigma1-S) z + x^-A z^2 ((z-1)-phi)^2.
# The first edge has height 2, slope sigma1, and an open-edge height-1
# support.  The second edge has slope S and polynomial z^2(z-1)^2.
parent_z = (
    x
    + x ** (sigma1 - S) * z
    + x ** (-A) * z**2 * ((z - 1) - phi) ** 2
)
parent = sp.expand(parent_z.subs(z, x**S * y))
check("the length-12 parent is a plane polynomial divisible by x", is_plane_x_multiple(parent))

poly_y = sp.Poly(parent, y)
row_orders = {
    height: x_order(poly_y.coeff_monomial(y**height))
    for height in range(1, 5)
}
check(
    "the parent is row D: first edge span 2 with open-edge support",
    row_orders[1] == sigma1
    and row_orders[2] == 2 * sigma1
    and row_orders[3] - row_orders[2] == S
    and row_orders[4] - row_orders[3] == S,
)

deeper_edge = sp.expand(x**A * (parent_z - x - x ** (sigma1 - S) * z)).subs(x, 0)
check(
    "the row-D deeper edge is z^2(z-1)^2",
    sp.expand(deeper_edge - z**2 * (z - 1) ** 2) == 0,
)
check(
    "the canonical deeper root z=1 has multiplicity exactly two",
    deeper_edge.subs(z, 1) == 0
    and sp.diff(deeper_edge, z).subs(z, 1) == 0
    and sp.diff(deeper_edge, z, 2).subs(z, 1) != 0,
)

charted = sp.expand(parent.subs(y, x ** (-S) * (1 + Y)))
advertised = (
    x
    + x ** (-A // 2) * (1 + Y)
    + x ** (-A) * (1 + Y) ** 2 * (Y - phi) ** 2
)
check("charted row-D parent has the advertised exact form", sp.expand(charted - advertised) == 0)

height_zero = []
for depth in range(0, R):
    truncation = sum(x**j for j in range(1, depth + 1))
    scale = x ** (depth + 1)
    node = sp.expand(charted.subs(Y, truncation + scale * Z))
    row = x_order(node.subs(Z, 0))
    height_zero.append(row)

    cleared = sp.expand(x ** (A - 2 * depth - 2) * node)
    edge = sp.expand(cleared.subs(x, 0))
    check(
        f"depth {depth}: non-anchored edge is the repeated polynomial (Z-1)^2",
        sp.expand(edge - (Z - 1) ** 2) == 0,
    )
    check(
        f"depth {depth}: repeated-root multiplicity equals full cluster degree 2",
        sp.Poly(edge, Z).degree() == 2
        and edge.subs(Z, 1) == 0
        and sp.diff(edge, Z).subs(Z, 1) == 0
        and sp.diff(edge, Z, 2).subs(Z, 1) != 0,
    )

expected_rows = [-A + 2 * depth + 2 for depth in range(R)]
check("height-zero valuations strictly rise by two", height_zero == expected_rows)
check("all twelve valuations remain below both the separator and t-anchor", height_zero[-1] < -A // 2 < 0)
check("twelve consecutive repeated refinements retain cluster degree two", len(height_zero) == R)

print("\nEXACT ROW-D COUNTEREXAMPLE: repeated-root degree descent is false as printed.")
print("The same construction works for every R; no bound from cluster degree exists.")
