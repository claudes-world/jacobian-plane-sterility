"""Closed-form fertility criterion + germ + residue ratio, swept against the
shipped cone solver.

Theory (zeta-operator, base x^{-k}, b = 0, c_min = Y + x^r Y^N, zeta = x^r Y^{N-1}):
  J(x^a phi(zeta), c_min) = x^{a-1}[(a + N*a*zeta)phi + r*zeta*(1+zeta)*phi']
  Germ equation L[phi] = lam with a = -k. On monomials:
  L[zeta^n] = (a + r n) zeta^n + (N a + r n) zeta^{n+1}
  => phi_0 = lam/a;  (a + r n) phi_n = -(N a + r(n-1)) phi_{n-1}.
  * terminates (phi polynomial, degree n*-1) iff N a + r(n*-1) = 0 has an
    integer solution n* >= 1:  n* = 1 + N k / r  =>  r | N k.
  * fatal diagonal zero iff a + r n = 0 for some n >= 1 with nonzero source:
    n = k/r  =>  r | k  (source coefficient -(N-1)k - r != 0 always).
  FERTILITY (b=0 germ exists with lam != 0)  <=>  r | N k  and  r NOT| k.
  Residue ratio: L[phi]=lam evaluated at zeta = 0 and -1 (phi'-coefficient
  vanishes at both):  phi(0) = -lam/k,  phi(-1) = lam/((N-1)k):
  ratio = -1/(N-1), independent of k and r.

This script:
  (1) builds phi for each criterion-fertile (N,k,r) in a grid, machine-checks
      J(x^{-k} phi, c_min) = lam x^{-1-k} exactly and the ratio -1/(N-1);
  (2) sweeps the shipped cone_solver.cone() over the grid and diffs verdicts
      against the criterion (windows widened with k and r; Vcap grown to
      reach the phi-degree so 'none' vs 'SOLVABLE' is honest);
  (3) N=2 control: criterion must reproduce the jewel ladder law
      (r even, k an odd multiple of r/2).
"""
import sys
import sympy as sp

sys.path.insert(0, '/home/claude/claudes-world/jc')
from cone_solver import cone

x, Y, lam = sp.symbols('x Y lam')

def criterion(N, k, r):
    return (N*k) % r == 0 and k % r != 0

def germ_phi(N, k, r):
    """Polynomial phi(zeta) solving the germ equation, or None."""
    if not criterion(N, k, r):
        return None
    a = -k
    nstar = 1 + N*k // r
    phis = [sp.Rational(1, a)*lam]
    for n in range(1, nstar):
        phis.append(sp.cancel(-(N*a + r*(n-1))*phis[n-1]/(a + r*n)))
    z = sp.Symbol('z')
    return sum(co*z**n for n, co in enumerate(phis))

results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))

z = sp.Symbol('z')
print("== (1) germ construction + exact pair equation + ratio, criterion-fertile cells ==")
GRID3 = [(3, k, r) for r in (3, 6, 9) for k in range(1, 8)]
GRID4 = [(4, k, r) for r in (4, 8) for k in range(1, 6)]
for (N, k, r) in GRID3 + GRID4:
    phi = germ_phi(N, k, r)
    if phi is None:
        continue
    c_min = Y + x**r*Y**N
    w = sp.expand((x**-k*phi).subs(z, x**r*Y**(N-1)))
    J = sp.diff(w, x)*sp.diff(c_min, Y) - sp.diff(w, Y)*sp.diff(c_min, x)
    ok_eq = sp.simplify(J - lam*x**(-1-k)) == 0
    ratio = sp.cancel(phi.subs(z, -1)/phi.subs(z, 0))
    ok_ratio = sp.simplify(ratio + sp.Rational(1, N-1)) == 0
    check(f"({N},{k},{r}) germ exact + ratio -1/{N-1}", ok_eq and ok_ratio,
          f"deg phi = {sp.degree(phi, z)}, ratio = {ratio}")

print("\n== (2) criterion vs shipped cone_solver (windows scaled to phi-degree) ==")
mism = []
for (N, k, r) in [(3, k, r) for r in (3, 6) for k in range(1, 6)] + \
                 [(4, k, 4) for k in range(1, 5)] + [(2, k, 2) for k in range(1, 6)]:
    want = 'SOLVABLE' if criterion(N, k, r) else 'none'
    deg = (1 + N*k//r) if (N*k) % r == 0 else 2
    Vcap = max(3, 2*deg)               # reach b = (N-1)*(deg-1) via c^nu * Y^{<=N-1}
    got = cone(N, k, r, Vcap=Vcap, lo_pad=3 + r, hi_pad=1 + r)
    tag = "OK" if got == want else "MISMATCH"
    if got != want:
        mism.append((N, k, r, want, got))
    print(f"  ({N},{k},{r}): criterion={want:8s} solver={got:8s} {tag}")
check("criterion matches cone_solver on full grid", not mism, f"mismatches={mism}")

print("\n== (3) N=2 jewel-ladder control ==")
bad2 = []
for k in range(1, 13):
    for r in range(1, 9):
        law = (r % 2 == 0) and (k % (r//2) == 0 if r >= 2 else False) \
              and ((k // (r//2)) % 2 == 1) if r % 2 == 0 else False
        crit = criterion(2, k, r)
        if bool(law) != bool(crit):
            bad2.append((k, r, law, crit))
check("criterion == jewel law (r even, k odd multiple of r/2) on k<=12, r<=8",
      not bad2, f"bad={bad2}")

fails = [n for n, ok in results if not ok]
print(f"\n{len(results)} checks, {len(results)-len(fails)} passed, {len(fails)} failed")
sys.exit(1 if fails else 0)
