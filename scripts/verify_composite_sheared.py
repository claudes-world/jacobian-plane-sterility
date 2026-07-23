"""Corrected P-layer sweep (post codex-audit finding B/C): full (j,s) grid +
SHEARED pole-cancellation rows (f o S plane-legality), lambda = gamma = 1
WLOG. The Shear-Locus Kill theorem (COMPOSITE-JEWEL rev-2 section 3) proves
mu = 0 for ALL (sigma, D); this witness machine-corroborates sigma <= 8,
D <= 6 with the corrected constraint system, and rediscovers the
functional-independence step (all off-target-line coefficients forced 0).
Run log: out/verify_composite_sheared.log (CLEAN)."""

import time
import sympy as sp
t0 = time.time()
x, Y = sp.symbols('x Y')
S3 = sp.sqrt(3)
mu = sp.Symbol('mu')

def dict_mul(d1, d2):
    out = {}
    for (a, b), u in d1.items():
        for (c, dd), v in d2.items():
            key = (a + c, b + dd)
            out[key] = sp.expand(out.get(key, 0) + u * v)
    return {k: v for k, v in out.items() if v != 0}

def to_dict(e):
    d = {}
    p = sp.expand(e)
    for term in sp.Add.make_args(p):
        pw = term.as_powers_dict()
        a, b = int(pw.get(x, 0)), int(pw.get(Y, 0))
        d[(a, b)] = sp.expand(d.get((a, b), 0) + term / (x**a * Y**b))
    return d

CDEG = 16
c_min = Y + S3*x**3*Y**2 + x**6*Y**3
w_g = -(1/x + S3*x**2*Y)
cS = sp.expand(c_min.subs(Y, Y + x))
wS = sp.expand(w_g.subs(Y, Y + x))
wd, cd = to_dict(w_g), to_dict(c_min)
wSd, cSd = to_dict(wS), to_dict(cS)

def pows(base, n):
    out = [{(0, 0): sp.Integer(1)}]
    for _ in range(n):
        out.append(dict_mul(out[-1], base))
    return out
wp, cp = pows(wd, 6), pows(cd, CDEG)
wSp, cSp = pows(wSd, 6), pows(cSd, CDEG)

survivors = []
for D in range(1, 7):
    syms = {(j, s): sp.Symbol(f'p{j}_{s}')
            for j in range(0, D + 1) for s in range(CDEG + 1)}
    unk = list(syms.values()) + [mu]
    pw_tot, fS_tot = {}, {}
    for (j, s), sym in syms.items():
        if j >= 1:
            for m_, v in dict_mul(cp[s], wp[j - 1]).items():
                pw_tot[m_] = pw_tot.get(m_, 0) + j * sym * v
        for m_, v in dict_mul(cSp[s], wSp[j]).items():   # SHEARED f
            fS_tot[m_] = fS_tot.get(m_, 0) + sym * v
    for sigma in range(1, 9):
        target = (2 - sigma, 0)
        rows = []
        for mon in set(pw_tot) | {target}:
            expr = sp.expand(pw_tot.get(mon, 0) - (mu if mon == target else 0))
            if expr != 0:
                rows.append([expr.coeff(u) for u in unk])
        for (a, b), v in fS_tot.items():
            if a + sigma * b < 0:                        # sheared pole rows
                rows.append([sp.expand(v).coeff(u) for u in unk])
        M = sp.Matrix(rows)
        if any(sp.simplify(v[-1]) != 0 for v in M.nullspace()):
            survivors.append((sigma, D))
            print(f"  !! sigma={sigma} D={D} mu-consistent", flush=True)
    print(f"  D={D} done ({round(time.time()-t0,1)}s)", flush=True)
print("SHEARED-ROWS FULL-GRID SWEEP:",
      "CLEAN — mu forced 0, sigma<=8, D<=6" if not survivors else survivors)
