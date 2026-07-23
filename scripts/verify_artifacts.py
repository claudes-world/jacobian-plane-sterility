"""A1-A6 artifact re-verification (Line Theorem handoff, Section 4 + T10).

Every check is exact (sympy, rational arithmetic). Exit nonzero on any failure.
Run: python3 verify_artifacts.py
"""
import sys
import sympy as sp

x, y, Y, lam, gam = sp.symbols('x y Y lam gam')

def J(f, g, v1, v2):
    return sp.simplify(sp.diff(f, v1) * sp.diff(g, v2) - sp.diff(f, v2) * sp.diff(g, v1))

results = []

def check(name, ok, detail=""):
    results.append((name, bool(ok), detail))
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))

# ---- A1: swap automorphism scaffold: w = y/x, c = x, P = WC, f = y
w = y / x
c = x
f = sp.expand(w * c)
check("A1 f = P(w,c) = y", sp.simplify(f - y) == 0)
check("A1 J(f,c) = -1 (kappa)", sp.simplify(J(f, c, x, y) + 1) == 0)
check("A1 pair eq J(w,c) = -x^-1", sp.simplify(J(w, c, x, y) + 1/x) == 0)

# ---- A2: nonconstant-u pair dying at P-layer: w = 1/x, c = x + x^3 y
w = 1 / x
c = x + x**3 * y
check("A2 J(w,c) = -x", sp.simplify(J(w, c, x, y) + x) == 0)

# ---- A3: jewel pair, rung nu*=0
u3 = 1 + x**3 + (x**2 + 2*x**5)*y + x**7*y**2
c = sp.expand(x * u3)
w = -(1 + 2*x**3 + 2*x**5*y) / x
check("A3 J(w,c) = x", sp.simplify(J(w, c, x, y) - x) == 0)
check("A3 jewel w^2 - 4c = x^-2", sp.simplify(w**2 - 4*c - x**-2) == 0)

# ---- A4: jewel pair, rung nu*=1
u4 = sp.Rational(-4, 3) - sp.Rational(2, 3)*x**3 + (x**6 + x**9)*y - sp.Rational(3, 8)*x**15*y**2
c = sp.expand(x * u4)
v4 = 1 + x**3 - sp.Rational(3, 4)*x**9*y
w = (sp.Rational(1, 3)*x**-3 + c/x) * v4
xi = x**-2
check("A4 J(w,c) = -x^3", sp.simplify(J(w, c, x, y) + x**3) == 0)
check("A4 jewel w^2 = (1/9)xi^3 + (1/2)c xi^2 - (3/2)c^3",
      sp.simplify(w**2 - (sp.Rational(1, 9)*xi**3 + sp.Rational(1, 2)*c*xi**2 - sp.Rational(3, 2)*c**3)) == 0)

# ---- A5: jewel pair, rung nu*=2 (blind-built)
u5 = 1 + x**3/2 + x**10*(1 + x**3)*y + sp.Rational(1, 2)*x**23*y**2
c = sp.expand(x * u5)
v5 = 1 + x**3 + x**13*y
B = sp.Rational(3, 32)*x**-5 - sp.Rational(1, 4)*c*x**-3 + c**2*x**-1
w = B * v5
Phi = sp.Rational(3, 32)*xi**2 - sp.Rational(1, 4)*xi*c + c**2
check("A5 J(w,c) = -(15/32)x^5", sp.simplify(J(w, c, x, y) + sp.Rational(15, 32)*x**5) == 0)
check("A5 jewel w^2 = Phi^2 (xi + 2c)", sp.simplify(w**2 - Phi**2*(xi + 2*c)) == 0)

# ---- A6: the (3,1,3) germ (chart variables x, Y)
c_min = Y + x**3*Y**3
w_germ = -lam*(x**-1 + sp.Rational(3, 2)*x**2*Y**2)
check("A6 J(w_germ,c_min) = lam x^-2", sp.simplify(J(w_germ, c_min, x, Y) - lam*x**-2) == 0)

# ---- A6b (T10): bare-tail level-1 resummation
# gamma^1 equation: J(w1, c_min) = -J(w_germ, gam*x) = -3 gam lam x^2 Y
w1 = -lam*gam*x**3*Y/(1 + x**3*Y**2)
lhs = J(w1, c_min, x, Y)
check("T10 w1 solves J(w1,c_min) = -3 gam lam x^2 Y",
      sp.simplify(lhs + 3*gam*lam*x**2*Y) == 0)
check("T10 J(w_germ, gam*x) = 3 gam lam x^2 Y",
      sp.simplify(J(w_germ, gam*x, x, Y) - 3*gam*lam*x**2*Y) == 0)

fails = [r for r in results if not r[1]]
print(f"\n{len(results)} checks, {len(results)-len(fails)} passed, {len(fails)} failed")
sys.exit(1 if fails else 0)
