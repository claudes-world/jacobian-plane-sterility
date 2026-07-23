"""Composite-jewel cell: shear pair + P-layer kill on the germ pair.

A. SHEAR PAIR exact: S: (x,Y) -> (x, Y + gamma x), Jacobian 1, fixes x, so
   (w_germ o S, c_min o S) solves the pair equation with the same lam, and
   c_min o S carries the gamma-x slot -> legitimate chart c. Its tail slots
   {(4,1),(5,0),(7,2),(8,1),(9,0)} are exactly the slots the full-knob
   marcher imposes at the gamma-rungs (buckets 2, 6, 10): the forced
   cascade is the shear orbit's shadow.

B. S-CONJUGATION REDUCTION + P-LAYER SWEEP: a scaffold on the shear pair
   pulls back to P on the bare germ pair:
       (master)  P_W(w_germ, c_min) = mu x^(2 - sigma)   [m_plane = sigma-2, k=1]
       (poles)   f = P(w_germ, c_min) polynomial in the sigma-realization
                 (every chart monomial (a,b) must have a + sigma*b >= 0)
   Both constraints are linear in P's coefficients p_{j,s} (P = sum p_js C^s W^j)
   and mu. For each (sigma, D) we build the FULL system over the exact
   monomial support (target row always included — an unreachable target
   reads 0 = mu) and test whether mu != 0 is consistent. Window guard: any
   mu-bearing null vector whose p-support touches the degree cap is
   reported, never silently trusted.

Sweep: sigma = 1..8, D = 1..6. Expected from theory: sigma = 1 dead by
Y-degree counting; sigma = 2 (m = 0) dead by pole-cancellation; sigma >= 3
is the genuine sliver question restricted to the shear locus.
"""
import time
import sympy as sp

x, Y, y = sp.symbols('x Y y')
lam, gam, mu = sp.symbols('lam gam mu')
S3 = sp.sqrt(3)

npass = nfail = 0
def check(label, ok, extra=""):
    global npass, nfail
    tag = "PASS" if ok else "FAIL"
    if ok: npass += 1
    else: nfail += 1
    print(f"{tag}  {label}" + (f"  [{extra}]" if extra else ""), flush=True)

def J(f, g, v=Y):
    return sp.expand(sp.diff(f, x) * sp.diff(g, v) - sp.diff(f, v) * sp.diff(g, x))

c_min = Y + S3 * x**3 * Y**2 + x**6 * Y**3
w_germ = -lam * (1 / x + S3 * x**2 * Y)

# ---- A. germ identity + shear pair --------------------------------------
check("A1: germ identity J(w_germ, c_min) = lam x^-2",
      sp.expand(J(w_germ, c_min) * x**2 - lam) == 0)

wS = w_germ.subs(Y, Y + gam * x)
cS = sp.expand(c_min.subs(Y, Y + gam * x))
check("A2: shear pair J(w o S, c o S) = lam x^-2 (gamma symbolic)",
      sp.expand(J(wS, cS) * x**2 - lam) == 0)

tails = sp.expand(cS - c_min - gam * x)
expected_tails = (2 * S3 * gam * x**4 * Y + S3 * gam**2 * x**5
                  + 3 * gam * x**7 * Y**2 + 3 * gam**2 * x**8 * Y
                  + gam**3 * x**9)
check("A3: shear tails = {(4,1),(5,0),(7,2),(8,1),(9,0)} slots exactly",
      sp.expand(tails - expected_tails) == 0)

w_p = wS.subs(Y, x * y)
c_p = sp.expand(cS.subs(Y, x * y))
check("A4: sigma=1 plane realization J = lam x^-1 and c legit (u(0,0)=gamma)",
      sp.expand(J(w_p, c_p, v=y) * x - lam) == 0
      and sp.expand(sp.expand(c_p / x).subs({x: 0, y: 0}) - gam) == 0)

# ---- B. P-layer sweep on the bare germ pair -----------------------------
print("== B: P-layer on (w_germ, c_min): sigma = 1..8, D = 1..6 ==", flush=True)
t0 = time.time()
CDEG = 14

def dict_mul(d1, d2):
    out = {}
    for (a, b), u in d1.items():
        for (c, dd), v in d2.items():
            key = (a + c, b + dd)
            out[key] = sp.expand(out.get(key, 0) + u * v)
    return {k: v for k, v in out.items() if v != 0}

wd = {(-1, 0): -lam, (2, 1): -S3 * lam}
cd = {(0, 1): sp.Integer(1), (3, 2): S3, (6, 3): sp.Integer(1)}
wpow = [{(0, 0): sp.Integer(1)}]
for i in range(1, 7):
    wpow.append(dict_mul(wpow[-1], wd))
cpow = [{(0, 0): sp.Integer(1)}]
for s in range(1, CDEG + 1):
    cpow.append(dict_mul(cpow[-1], cd))

# precompute per (j, s): the dicts c^s * w^j and c^s * w^(j-1)
prod = {}
for jj in range(0, 7):
    for s in range(CDEG + 1):
        prod[(jj, s)] = dict_mul(cpow[s], wpow[jj])

survivors = []
for D in range(1, 7):
    syms = {(j, s): sp.Symbol(f'p{j}_{s}')
            for j in range(0, D + 1) for s in range(CDEG + 1)}
    unk = list(syms.values()) + [mu]
    # linear-form dicts: monomial -> linear expr in p_{j,s}
    pw_tot, f_tot = {}, {}
    for (j, s), sym in syms.items():
        if j >= 1:
            for m_, v in prod[(j - 1, s)].items():
                pw_tot[m_] = pw_tot.get(m_, 0) + j * sym * v
        for m_, v in prod[(j, s)].items():
            f_tot[m_] = f_tot.get(m_, 0) + sym * v
    for sigma in range(1, 9):
        target = (2 - sigma, 0)
        rows = []
        support = set(pw_tot) | {target}
        for mon in support:
            expr = sp.expand(pw_tot.get(mon, 0) - (mu if mon == target else 0))
            if expr != 0:
                rows.append([expr.coeff(u) for u in unk])
        for (a, b), v in f_tot.items():
            if a + sigma * b < 0:               # plane pole: must vanish
                rows.append([sp.expand(v).coeff(u) for u in unk])
        M = sp.Matrix(rows)
        ns = M.nullspace()
        mu_vecs = [v for v in ns if sp.simplify(v[-1]) != 0]
        if mu_vecs:
            # window guard: does any mu-solution need p-support at the cap?
            capped = any(any(v[i] != 0 and list(syms)[i][1] >= CDEG - 1
                             for i in range(len(unk) - 1)) for v in mu_vecs)
            survivors.append((sigma, D, len(mu_vecs), capped))
            print(f"   !! sigma={sigma} (m={sigma-2}), D={D}: mu-consistent "
                  f"null dim {len(mu_vecs)}{' [AT CAP — enlarge]' if capped else ''}",
                  flush=True)
    print(f"   D={D} swept ({round(time.time() - t0, 1)}s)", flush=True)

check("B: mu forced 0 for all sigma <= 8, D <= 6 on the germ pair "
      "(shear locus scaffold-dead)", not survivors,
      f"survivors: {survivors}" if survivors else "clean")

print(f"\n{npass + nfail} checks, {npass} passed, {nfail} failed "
      f"({round(time.time() - t0, 1)}s)", flush=True)
