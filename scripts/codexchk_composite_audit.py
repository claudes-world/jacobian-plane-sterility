"""Exact checks for codex-audit-composite-2026-07-23.md.

This script is intentionally independent of the shipped JSON conclusions.
It uses the shipped marcher only for the section-E window hunt.
"""

import json
from pathlib import Path

import sympy as sp

from classify_stratum import classify


x, Y, y, z = sp.symbols("x Y y z")
t = sp.sqrt(3)


def jac(f, g, second=Y):
    return sp.expand(
        sp.diff(f, x) * sp.diff(g, second)
        - sp.diff(f, second) * sp.diff(g, x)
    )


def monomial_dict(expr):
    """Convert a Laurent polynomial in x and polynomial in Y to a dict."""
    out = {}
    for term in sp.Add.make_args(sp.expand(expr)):
        powers = term.as_powers_dict()
        a = int(powers.get(x, 0))
        b = int(powers.get(Y, 0))
        coeff = sp.cancel(term / (x**a * Y**b))
        out[(a, b)] = sp.cancel(out.get((a, b), 0) + coeff)
    return {m: c for m, c in out.items() if c != 0}


def section_a():
    print("== A: minimal-chart kill ==")
    L, G = sp.symbols("L G", nonzero=True)
    E = 1 + t * z + z**2
    cmin = Y + t * x**3 * Y**2 + x**6 * Y**3
    wg1 = -L * (x**-1 + t * x**2 * Y)
    assert sp.expand(jac(wg1, cmin) - L * x**-2) == 0

    # At bucket 2 put the unknown weight-3 layer in the form x^3 Phi(z).
    # Direct differentiation gives
    #   J(x^3 Phi(z), cmin) = 3*x^2*(z*E*Phi)'.
    # The gamma source is +t*G*L*x^2.  Integration would require
    # E*Phi = -t*G*L/3, impossible for polynomial Phi when G*L != 0.
    deg = 12
    coeffs = sp.symbols(f"a0:{deg + 1}")
    Phi = sum(coeffs[i] * z**i for i in range(deg + 1))
    bucket2 = sp.Poly(
        sp.expand(3 * sp.diff(z * E * Phi, z) + t * G * L), z
    )
    rows = [
        [bucket2.coeff_monomial(z**i).coeff(u) for u in coeffs]
        + [-bucket2.coeff_monomial(z**i).subs({u: 0 for u in coeffs})]
        for i in range(deg + 3)
    ]
    aug = sp.Matrix(rows)
    assert sp.linsolve(
        (aug[:, :-1], aug[:, -1]), coeffs
    ) is sp.EmptySet
    print("k=1 bucket 2: 3*(z*E*Phi)' = -sqrt(3)*G*L has no polynomial solution")
    print("canonical obstruction: -sqrt(3)*G*L")

    # k=4: at the first gamma rung the unknown weight-0 layer is Phi(z),
    # whose cmin bracket is 3*x^-1*z*E*Phi'.  It has zero constant term.
    phi4 = L * (
        -sp.Rational(1, 4)
        + 2 * t * z
        + sp.Rational(27, 2) * z**2
        + 9 * t * z**3
        + sp.Rational(27, 4) * z**4
    )
    obstruction4 = sp.expand(G * sp.diff(phi4, z).subs(z, 0))
    assert obstruction4 == 2 * t * G * L
    print("k=4 bucket -1 constant obstruction:", obstruction4)

    expected = {
        "out/mc_k1_minimal.json": {
            -5: [],
            -4: [[(0, 1), (3, 2), (6, 3)]],
            -3: [],
            -2: [],
        },
        "out/mc_k1_min_e1.json": {
            -5: [],
            -4: [[(0, 1), (3, 2), (6, 3)]],
            -3: [],
            -2: [],
        },
        "out/mc_k4_minimal.json": {
            -8: [],
            -7: [[(0, 2), (3, 3), (6, 4), (9, 5), (12, 6)]],
            -6: [],
            -5: [],
        },
        "out/mc_k1_full_deep.json": {
            -5: [],
            -4: [[(0, 1), (3, 2), (6, 3)]],
            -3: [],
            -2: [],
        },
    }
    for filename, want in expected.items():
        data = json.loads(Path(filename).read_text())
        got = {}
        for rec in data["buckets"]:
            if "kernel_support" in rec:
                got[rec["bucket"]] = [
                    [tuple(m) for m in support]
                    for support in rec["kernel_support"]
                ]
        assert got == want, (filename, got, want)
        print(filename, "kernel supports match c_min powers")


def section_b_c():
    print("\n== B/C: shear law, polynomiality, and corrected P-layer ==")
    W, C = sp.symbols("W C")
    cmin = Y + t * x**3 * Y**2 + x**6 * Y**3
    wg = -(x**-1 + t * x**2 * Y)  # lambda = 1
    cS = sp.expand(cmin.subs(Y, Y + x))  # gamma = 1
    wS = sp.expand(wg.subs(Y, Y + x))
    assert sp.expand(jac(wS, cS) - x**-2) == 0

    # Relevant counterexample to the claimed polynomiality transfer:
    # P(W,C)=W+W^2 C at sigma=2.
    P = W + W**2 * C
    f0 = sp.expand(P.subs({W: wg, C: cmin}, simultaneous=True))
    fS = sp.expand(P.subs({W: wS, C: cS}, simultaneous=True))
    f0_plane = sp.expand(f0.subs(Y, x**2 * y))
    fS_plane = sp.expand(fS.subs(Y, x**2 * y))
    assert f0_plane.coeff(x, -1) == -1
    assert all(
        int(term.as_powers_dict().get(x, 0)) >= 0
        for term in sp.Add.make_args(fS_plane)
    )
    print("counterexample at sigma=2: P=W+W^2*C")
    print("bare f has term -x^-1; sheared f is polynomial")

    # Plane exponent bookkeeping.
    k, sigma = sp.symbols("k sigma", integer=True)
    print("m_plane = sigma + (-1-k) =", sp.expand(sigma - 1 - k))
    assert (sigma - 1 - k).subs(k, 1) == sigma - 2

    # Corrected exact sweep.  Put q=1+sqrt(3)z and
    # h=z(1+sqrt(3)z+z^2).  For a term C^s W^J in P,
    # its P_W contribution has weight -(J-1+3s).
    q = 1 + t * z
    h = z * (1 + t * z + z**2)
    bare_survivors = []
    shear_survivors = []
    for sig in range(1, 9):
        n = sig - 2
        for D in range(1, 7):
            active = [
                (J, s)
                for J in range(1, D + 1)
                for s in range(0, 15)
                if J - 1 + 3 * s == n
            ]
            aa = sp.symbols(f"a_{sig}_{D}_0:{len(active)}")
            MU = sp.Symbol(f"mu_{sig}_{D}")
            master = -MU
            fpoly = 0
            for acoef, (J, s) in zip(aa, active):
                master += (
                    J * acoef * (-1) ** (J - 1) * q ** (J - 1) * h**s
                )
                fpoly += acoef * (-1) ** J * q**J * h**s
            master_poly = sp.Poly(sp.expand(master), z)
            rows_master = [
                master_poly.coeff_monomial(z**ell)
                for ell in range(master_poly.degree() + 1)
            ] if master_poly.as_expr() != 0 else []
            f_poly = sp.Poly(sp.expand(fpoly), z)
            max_fdeg = max(0, f_poly.degree()) if f_poly.as_expr() != 0 else 0
            rows_bare = list(rows_master)
            rows_shear = list(rows_master)
            for ell in range(max_fdeg + 1):
                coeff = f_poly.coeff_monomial(z**ell)
                if (sig + 3) * ell < sig - 1:
                    rows_bare.append(coeff)
                if 4 * ell < sig - 1:
                    rows_shear.append(coeff)
            unknowns = list(aa) + [MU]

            def has_mu(rows):
                matrix = sp.Matrix(
                    [[sp.expand(row).coeff(u) for u in unknowns] for row in rows]
                )
                return any(v[-1] != 0 for v in matrix.nullspace())

            if has_mu(rows_bare):
                bare_survivors.append((sig, D))
            if has_mu(rows_shear):
                shear_survivors.append((sig, D))

    assert bare_survivors == []
    assert shear_survivors == []
    print("corrected exact sweep: no mu survivor for sigma<=8, D<=6")
    print("bare-system sweep:     no mu survivor for sigma<=8, D<=6")

    # The target-weight equation itself bounds every relevant C-degree.
    # For sigma<=8, J-1+3s=sigma-2 implies s<=2.
    max_s = max(
        [
            s
            for sig in range(2, 9)
            for J in range(1, 7)
            for s in range(15)
            if J - 1 + 3 * s == sig - 2
        ],
        default=0,
    )
    assert max_s == 2
    print("target-weight bound in swept range: degree_C(P_W) <=", max_s)


def section_d():
    print("\n== D: classifier ==")
    # First possible pole depth is k=9: N=3, d=2, r=54, q=r/d=27.
    # At t=sqrt(3), phi=1+sqrt(3)xi and lambda=-9 solve the germ ODE.
    support = {
        (1, 0): sp.Integer(1),
        (1, 1): sp.Integer(1),
        (29, 2): t,
        (57, 3): sp.Integer(1),
    }
    lbl8, _, info8 = classify(support, kmax=8)
    lbl9, _, info9 = classify(support, kmax=9)
    assert lbl8 == "E3" and lbl9 == "E4"
    assert info8["fertile_ks"] == [] and info9["fertile_ks"] == [9]

    E = 1 + t * z + z**2
    phi = 1 + t * z
    germ_ode = sp.expand(
        27 * z * E * sp.diff(phi, z)
        - 9 * (E + z * sp.diff(E, z)) * phi
    )
    assert germ_ode == -9
    print("same valid support: kmax=8 ->", lbl8, "; kmax=9 ->", lbl9)
    print("exact k=9 germ ODE output:", germ_ode)


def section_e():
    print("\n== E: exact non-shear b=2,3 window survivor ==")
    max_bucket = 5
    # This candidate was found by an independent triangular coefficient
    # solve with lambda=gamma=1 and the free b=3 slot k_(7,3)=1.  Every
    # unlisted admissible chart-tail slot through weight 7 is set to zero.
    cdict = {
        (0, 1): sp.Integer(1),
        (3, 2): t,
        (6, 3): sp.Integer(1),
        (4, 2): t / 2,
        (7, 3): sp.Integer(1),
        (2, 1): sp.Rational(2, 9),
        (3, 1): sp.Rational(7, 27),
        (9, 3): -sp.Rational(14, 27),
        (1, 0): sp.Integer(1),
        (4, 1): sp.Rational(47, 162) + 2 * t,
        (7, 2): 3 + 5 * t / 18,
        (2, 0): -1 + 47 * t / 729,
        (5, 1): sp.Rational(10, 81) - t,
        (3, 0): sp.Rational(11, 9) - 2425 * t / 78732,
        (6, 1): sp.Rational(49, 2916) + t,
        (4, 0): -sp.Rational(4, 9) + 9713 * t / 157464,
        (7, 1): t / 27 + sp.Rational(1153, 8748),
    }
    wspec = {
        (-1, 0): -1,
        (2, 1): -t,
        (3, 1): -t / 3,
        (1, 0): -sp.Rational(2, 9),
        (4, 1): t / 9,
        (2, 0): -sp.Rational(7, 54),
        (5, 1): t / 9,
        (3, 0): -t - sp.Rational(13, 162),
        (6, 1): -2 * t / 27,
        (4, 0): -sp.Rational(8, 81) + 2 * t / 3,
        (7, 1): 16 * t / 243,
        (5, 0): sp.Rational(277, 4374) - 5 * t / 9,
        (8, 1): -70 * t / 2187,
        (6, 0): -sp.Rational(529, 6561) + 4 * t / 27,
        (9, 1): 62 * t / 6561,
    }
    cexpr = sum(value * x**a * Y**b for (a, b), value in cdict.items())
    wexpr = sum(value * x**a * Y**b for (a, b), value in wspec.items())
    residual = monomial_dict(sp.expand(jac(wexpr, cexpr) - x**-2))
    low_residual = {
        m: value
        for m, value in residual.items()
        if m[0] - 3 * m[1] <= max_bucket
    }
    assert low_residual == {}
    assert cdict[(4, 2)] == t / 2
    assert cdict[(7, 3)] == 1
    print("lambda=gamma=1; free k_(7,3)=1; all other free knobs=0")
    print("nonzero tails (all unlisted slots through weight 7 are zero):")
    for m, value in sorted(
        (
            (m, v)
            for m, v in cdict.items()
            if m not in {(0, 1), (3, 2), (6, 3)}
        ),
        key=lambda item: (
            item[0][0] - 3 * item[0][1],
            item[0][1],
            item[0][0],
        ),
    ):
        print(f"  {m}: {value}")
    print("residual is exactly zero through bucket", max_bucket)
    print("first residual weight:", min(m[0] - 3 * m[1] for m in residual))


if __name__ == "__main__":
    section_a()
    section_b_c()
    section_d()
    section_e()
    print("\nALL EXACT CHECKS PASSED")
