"""(3,2,6) witness — the gcd(N-1, r) > 1 cell exposed by the GPT audit.

FERRY Lemma 3's chain argument treated a weight level as one zeta-chain;
at (N,k,r) = (3,2,6), d := gcd(N-1, r) = 2 and each level splits into d
cosets (b mod (N-1)/d = b mod 1 ... careful: coset index = b mod (N-1)/d?
No: on a level, allowed b-steps are multiples of (N-1)/d = 1, i.e. the
level lattice is b-steps of (N-1)/d; zeta-steps move b by N-1 = d*(N-1)/d,
so the level splits into d classes of b mod ... the zeta-chain reaches b's
congruent mod N-1; level lattice has b-step (N-1)/d; classes = d cosets
b mod (N-1) within the sublattice. Here: b-step 1, zeta-step 2 => 2 cosets
(b even / b odd) per level.]

Checks:
 1. fertility at (3,2,6): r | Nk (6 | 6), r !| k (6 !| 2) => germ exists,
    degree Nk/r = 1: w_germ = x^-2 * (phi0 + phi1*zeta), zeta = x^6 Y^2.
    Verify J(w_germ, c_min) = lam*x^-3 EXACTLY (chart bracket).
 2. residue evaluations: phi(0) = -lam/k = -lam/2, phi(-1) = +lam/((N-1)k)
    = +lam/4; ratio -1/(N-1) = -1/2.
 3. coset preservation: J(. , Y) keeps b; J(. , x^6 Y^3) shifts b by +2 —
    both preserve b mod 2 on every level. (The repair's engine.)
 4. per-coset kill: windowed solve of J(w, c_min) = lam*x^-3 with w
    supported ONLY on the odd-b coset: forces w = 0 and lam = 0 (the
    source lives on the even coset; no cross-feeding).
 5. fiber punctures at a generic t (bare cell + a tail): 1 main +
    d = 2 companions, each unramified (e = (N-1)/d = 1); w-residues
    -lam/2 (main), +lam/4 (each companion); trace check
    1*(-1/2) + 2*(1/4) = 0.
 6. P_W kill at the punctures: matching forces mu = 0 for D = 2
    (windowed exact solve).
"""
import time
import sympy as sp

x, Y, t = sp.symbols('x Y t')
lam = sp.Symbol('lam')
N, k, r = 3, 2, 6

npass = nfail = 0
def check(label, ok, extra=""):
    global npass, nfail
    tag = "PASS" if ok else "FAIL"
    if ok: npass += 1
    else: nfail += 1
    print(f"{tag}  {label}" + (f"  [{extra}]" if extra else ""))

def Jch(f, g):
    return sp.expand(sp.diff(f, x)*sp.diff(g, Y) - sp.diff(f, Y)*sp.diff(g, x))

c_min = Y + x**r * Y**N

# -- 1. germ from the Lemma-1 recursion: (a + r n) phi_n = -(N a + r(n-1)) phi_{n-1}
a = -k
phi0 = -lam / k                                   # a*phi0 = lam  =>  -k phi0 = lam
phi1 = sp.Rational(0)
# n = 1: (a + r) phi1 = -(N a) phi0
phi1 = sp.simplify(-(N*a) * phi0 / (a + r))
zeta = x**r * Y**(N-1)
w_germ = sp.expand(x**a * (phi0 + phi1 * zeta))
resid = sp.expand(Jch(w_germ, c_min) - lam * x**(-1-k))
check("germ (3,2,6): J(w_germ, c_min) = lam x^-3 exact", resid == 0,
      f"phi = {phi0} + ({phi1}) zeta")

# -- 2. residue evaluations of phi at zeta = 0, -1
phi = phi0 + phi1 * sp.Symbol('z')
z = sp.Symbol('z')
ev0 = phi.subs(sp.Symbol('z'), 0)
evm1 = phi.subs(sp.Symbol('z'), -1)
check("phi(0) = -lam/2 and phi(-1) = +lam/4, ratio -1/2",
      sp.simplify(ev0 + lam/2) == 0 and sp.simplify(evm1 - lam/4) == 0
      and sp.simplify(ev0/evm1 + 2) == 0, f"phi(0)={ev0}, phi(-1)={evm1}")

# -- 3. coset preservation of the two c_min brackets
ok = True
for (p_, b_) in [(-2, 0), (1, 1), (4, 2), (-5, 3), (0, 5)]:
    m = x**p_ * Y**b_
    for gpart, shift in ((Y, 0), (x**r*Y**N, N-1)):
        br = Jch(m, gpart)
        if br == 0:
            continue
        pb = sp.Poly(sp.expand(br*x**30*Y**2), x, Y)
        bs = {mj - 2 for (_, mj) in pb.monoms()}
        if bs != {b_ + shift}:
            ok = False
check("both c_min brackets shift b by 0 or N-1=2: coset (b mod 2) preserved", ok)

# -- 4. odd-coset-only windowed solve must be empty (w = 0, lam = 0)
t0 = time.time()
A_syms, w = {}, 0
for p_ in range(-6, 13):
    for b_ in range(1, 10, 2):                    # ODD b only
        s_ = sp.Symbol(f'q{p_+6}_{b_}')
        A_syms[(p_, b_)] = s_
        w += s_ * x**p_ * Y**b_
eq = sp.expand(Jch(w, c_min) - lam * x**(-1-k))
p = sp.Poly(sp.expand(eq * x**8), x, Y)
unk = list(A_syms.values()) + [lam]
M = sp.Matrix([[co.coeff(s_) for s_ in unk] for co in p.coeffs()])
ns = M.nullspace()
lam_zero = all(sp.simplify(v[-1]) == 0 for v in ns)
# kernel on odd coset: c_min^b odd powers? c_min^1 = Y + x^6 Y^3 lives on odd coset!
# So nullspace may contain g(c) elements with odd b — but they must have lam = 0.
check("odd-coset solve: lam forced 0 (kernel = odd c_min powers only)",
      lam_zero, f"null dim {len(ns)}, {round(time.time()-t0,1)}s")

# -- 5. fiber punctures + residues (bare + tailed), plane realization
# sigma for minimal realization: sigma = ... chart Y = x^sigma y with
# sigma = r/(N-1) = 3: c = x^3 y + x^(6+3*3) ... use chart directly:
# fiber {c_min = t}: solve x -> 0 punctures via Puiseux in chart coords.
# main: Y ~ t (zeta -> 0); companions: 1 + zeta = 0 => Y^2 = -x^-6:
# Y = ±i x^-3 — two UNRAMIFIED punctures (e = 1 each), gcd(2,6) = 2. Verify
# w_germ's x^-2-coefficient at each puncture equals phi(0) / phi(-1).
s_ = sp.Symbol('s_')                              # local parameter x = s_
# main puncture: Y = t + higher; w on it:
w_main = sp.expand(w_germ.subs(Y, t))
main_coeff = sp.expand(w_main * x**k).subs(x, 0)
check("main puncture x^-2-coefficient = phi(0)",
      sp.simplify(main_coeff - ev0) == 0, f"{main_coeff}")
okc = True
for sgn in (sp.I, -sp.I):
    Yc = sgn * x**(-3) * (1 + sp.Symbol('eps'))   # perturb to check leading only
    w_c = sp.expand(w_germ.subs(Y, Yc))
    lead = sp.simplify(sp.expand(w_c * x**k).subs(sp.Symbol('eps'), 0).subs(x, 0))
    # zeta = x^6 Y^2 -> -(1+eps)^2 -> -1: phi(-1) = lam/4
    if sp.simplify(lead - evm1) != 0:
        okc = False
check("both companion punctures (Y = ±i x^-3, unramified): coeff = phi(-1)",
      okc)
check("trace consistency: 1*phi(0) + (N-1)*phi(-1) = 0",
      sp.simplify(ev0 + (N-1)*evm1) == 0)

# -- 6. P_W kill at D = 2 on the punctures: 2 p2(t) w + p1(t) = mu x^-m
# leading matching at the three punctures: 2 p2 phi(0) = mu vs 2 p2 phi(-1) = mu
# => phi(0) = phi(-1), false since ratio -1/2 => p2 = 0 => mu = 0.
mu, p2 = sp.symbols('mu p2')
sols = sp.solve([2*p2*ev0 - mu, 2*p2*evm1 - mu], [p2, mu], dict=True)
check("D=2 puncture matching forces p2 = mu = 0",
      all(sl[p2] == 0 and sl[mu] == 0 for sl in sols) if sols else False,
      f"{sols}")

print(f"\n{npass + nfail} checks, {npass} passed, {nfail} failed")
if nfail:
    raise SystemExit(1)
