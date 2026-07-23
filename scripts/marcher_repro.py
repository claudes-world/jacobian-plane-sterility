"""C3 bucket-marcher, transcribed VERBATIM from the handoff doc (Section 7),
with non-invasive instrumentation: per-bucket JSON records, no-silent-skip
accounting, canonical cokernel normalization report, and a T11 diff.

Algorithm untouched: impose-and-continue, knob weights 0..12, buckets -4..11,
tower windows b<=15, a in [-3,25]. Run: python3 marcher_repro.py
"""
import json, time, sys
import sympy as sp

lam, gam, E1 = sp.symbols('lam gam E1')
knobs = {}
cd = {(0,1): sp.Integer(1), (3,3): sp.Integer(1), (1,0): gam}   # c_min + gamma*x
for b in range(0,4):
    for a in range(0,12):
        wt = 2*a-3*b
        if 0 <= wt <= 12 and (a,b) not in [(0,1),(3,3),(1,0),(0,0)]:
            s = sp.Symbol(f'k{a}{b}'); knobs[(a,b)] = s; cd[(a,b)] = s

def Jb(d1, d2):
    out = {}
    for (a,b),u in d1.items():
        for (c,d),v in d2.items():
            co = a*d - b*c
            if co:
                m = (a+c-1, b+d-1); out[m] = out.get(m, 0) + co*u*v
    return out

cmin = {(0,1): sp.Integer(1), (3,3): sp.Integer(1)}
tails = {m:v for m,v in cd.items() if m not in cmin}
w = {(0,1): E1, (3,3): E1}          # kernel injection E1*c_min
sub = {}
records = []
for beta in range(-4, 12):
    t0 = time.time()
    Jt = Jb(w, tails)
    src = {m: sp.expand(v.subs(sub)) for m,v in Jt.items() if 2*m[0]-3*m[1] == beta}
    if beta == -4: src[(-2,0)] = src.get((-2,0), 0) - lam
    mons = [(a,b) for b in range(0,16) for a in range(-3,26) if 2*a-3*b == beta+2 and (a,b) != (0,0)]
    us = {m: sp.Symbol(f'u{beta}_{m[0]}_{m[1]}') for m in mons}
    Jc = Jb({m: us[m] for m in mons}, cmin)
    eqm = sorted(set(list(Jc.keys()) + list(src.keys())))
    M = sp.Matrix([[sp.expand(Jc.get(m,0)).coeff(us[mm]) for mm in mons] for m in eqm])
    s = sp.Matrix([-src.get(m, 0) for m in eqm])
    rec = {"bucket": beta, "n_mons": len(mons), "n_eqm": len(eqm),
           "cokernel_dim": 0, "nonzero_conds": 0, "shield": None,
           "core": None, "core_canonical": None, "cond": None, "terminal": False}
    for vv in M.T.nullspace():                      # cokernel = obstruction functionals
        rec["cokernel_dim"] += 1
        cond = sp.expand(sp.expand((vv.T*s)[0]).subs(sub))
        if cond == 0: continue
        rec["nonzero_conds"] += 1
        fresh = [knobs[m] for m in knobs if 2*m[0]-3*m[1] == beta+1 and knobs[m] not in sub]
        core = sp.factor(cond.subs({f: 0 for f in fresh}))
        target = next((f for f in fresh if cond.has(f) and sp.degree(sp.Poly(cond, f)) == 1), None)
        print(f'bucket {beta}: shield={target}   core(fresh->0) = {core}')
        # canonical normalization: scale vv so its entry on the lowest
        # equation monomial with nonzero entry is 1, then re-extract the core
        piv = next(i for i in range(len(eqm)) if vv[i] != 0)
        vcan = vv / vv[piv]
        cond_can = sp.expand(sp.expand((vcan.T*s)[0]).subs(sub))
        core_can = sp.factor(cond_can.subs({f: 0 for f in fresh}))
        rec.update(shield=str(target), core=str(core), cond=str(cond),
                   core_canonical=str(core_can),
                   pivot_monomial=list(eqm[piv]), pivot_entry=str(vv[piv]),
                   fresh=[str(f) for f in fresh])
        if target is not None:
            sub[target] = sp.expand(sp.solve(cond, target)[0].subs(sub))
        else:
            print(f'  TERMINAL at bucket {beta}: {sp.factor(cond)}')
            rec["terminal"] = True
    s2 = sp.Matrix([sp.expand(xx.subs(sub)) for xx in s])
    solp, params = M.gauss_jordan_solve(s2)
    solp = solp.subs({p: 0 for p in params})        # particular solution, kernel params -> 0
    for mm, val in zip(mons, list(solp)):
        val = sp.expand(val.subs(sub))
        if val != 0: w[mm] = sp.expand(w.get(mm, 0) + val)
    rec["seconds"] = round(time.time() - t0, 2)
    records.append(rec)
    print(f'  [bucket {beta}: {rec["n_mons"]} unknowns, {rec["n_eqm"]} eqs, '
          f'coker {rec["cokernel_dim"]}, {rec["seconds"]}s]')

# ---- no-silent-skip accounting + T11 diff -------------------------------
EXPECTED = {
    0:  ("k53", "0"),
    1:  ("k42", "-27*gam*lam"),
    2:  ("k31", "-(177147/21505)*k32**2*lam"),
    3:  ("k20", "(531441/13585)*k21*k32*lam"),
    4:  ("k41", "45*lam*(k21**2 - 4*gam*k32)"),
    5:  ("k30", "-(39366/4301)*lam*(27*gam*k21 + k32**3 - 3*k32*k63)"),
    6:  ("k51", "-(177147/2470)*lam*(9*gam**2 + 2*k32*k52)"),
    7:  ("k40", "(165/2)*lam*(-2*gam*k32**2 + 2*gam*k63 + k21**2*k32 - 2*k21*k52 + k32*k73)"),
    8:  ("k61", "(6561/391)*lam*(36*gam*k52 + 9*k21**3 + 3*k21*k73 - 6*k32**2*k63 + 12*k32*k62 - 2*k63**2)"),
    9:  ("k50", "(19683/1235)*lam*(-54*gam**2*k32 + 27*gam*k21**2 + 27*gam*k73 + 18*k21*k32*k63 - 18*k21*k62 - 6*k32**2*k52 + 10*k32*k83 + 6*k52*k63)"),
    10: ("k71", "11*lam*(-90*gam*k32*k63 + 90*gam*k62 + 15*k21**2*k63 + 10*k21*k83 - 15*k32**2*k73 + 30*k32*k72 + 15*k52**2 - 9*k63*k73)"),
    11: ("k60", "(6561/782)*lam*(108*gam**3 - 108*gam*k21*k63 + 72*gam*k32*k52 - 60*gam*k83 + 18*k21**3*k32 - 18*k21**2*k52 - 30*k21*k32*k73 + 36*k21*k72 - 8*k32**3*k63 + 12*k32**2*k62 + 8*k32*k63**2 - 21*k32*k93 - 18*k52*k73 - 12*k62*k63)"),
}
ns = {str(v): v for v in list(knobs.values()) + [lam, gam, E1]}
print("\n==== T11 diff ====")
ok_all = True
by_bucket = {r["bucket"]: r for r in records}
for bkt in range(-4, 0):
    r = by_bucket[bkt]
    quiet = r["nonzero_conds"] == 0 and not r["terminal"]
    ok_all &= quiet
    print(f"bucket {bkt}: {'quiet (consistent), as expected' if quiet else 'UNEXPECTED OBSTRUCTION'}")
for bkt in range(0, 12):
    r = by_bucket[bkt]
    exp_shield, exp_core = EXPECTED[bkt]
    got_shield = r["shield"]
    shield_ok = (got_shield == exp_shield)
    core_run = sp.sympify(r["core"], locals=ns) if r["core"] is not None else None
    core_exp = sp.sympify(exp_core, locals=ns)
    if core_run is None:
        core_ok, ratio = False, None
    elif sp.simplify(core_run - core_exp) == 0:
        core_ok, ratio = True, sp.Integer(1)
    elif core_exp != 0 and core_run != 0 and sp.simplify(core_run/core_exp).is_Rational:
        core_ok, ratio = True, sp.simplify(core_run/core_exp)   # scale-equivalent
    else:
        core_ok, ratio = False, None
    one_cond = r["nonzero_conds"] == 1
    ok_all &= shield_ok and core_ok and one_cond and not r["terminal"]
    print(f"bucket {bkt}: shield {got_shield} {'==' if shield_ok else '!='} {exp_shield}; "
          f"core {'MATCH' if core_ok else 'MISMATCH'}"
          + (f" (scale {ratio})" if core_ok and ratio != 1 else "")
          + ("" if one_cond else f"; nonzero_conds={r['nonzero_conds']} (expected 1)"))

with open("out/repro_t11.json", "w") as fh:
    json.dump(records, fh, indent=1)
print(f"\nRESULT: {'REPRODUCED' if ok_all else 'DIVERGED'}; records -> out/repro_t11.json")
sys.exit(0 if ok_all else 1)
