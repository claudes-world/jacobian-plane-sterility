"""Value-level identification: the deep march's forced cascade == the shear.

Reads out/mc_k1_full_deep2.json (full-knob march, buckets -> 11, shields_all
recording), fixpoint-resolves the 25 forced imposition values (recorded
strings are imposition-time snapshots that may reference later-imposed
knobs), sets the 23 genuinely-free knobs to 0, and checks the result equals
the shear expansion c_min o (Y -> Y + gamma x) exactly:
   k4_1 = 2 sqrt3 gamma, k5_0 = sqrt3 gamma^2, k7_2 = 3 gamma,
   k8_1 = 3 gamma^2, k9_0 = gamma^3, all other imposed slots 0.
Exit 0 iff the match is exact. (See COMPOSITE-JEWEL-2026-07-23.md section 2.)
"""
import json, sys
import sympy as sp

d = json.load(open('out/mc_k1_full_deep2.json'))
gam = sp.Symbol('gam'); S3 = sp.sqrt(3)
shear = {"k4_1": 2*S3*gam, "k5_0": S3*gam**2, "k7_2": 3*gam,
         "k8_1": 3*gam**2, "k9_0": gam**3}
imposed = {}
for b in d["buckets"]:
    for slot, val in b.get("shields_all", []):
        imposed[slot] = sp.sympify(val)
symmap = {sp.Symbol(k): v for k, v in imposed.items()}
for _ in range(30):
    changed = False
    for k in list(imposed):
        v2 = sp.expand(sp.radsimp(imposed[k].subs(symmap)))
        if v2 != imposed[k]:
            imposed[k] = v2; symmap[sp.Symbol(k)] = v2; changed = True
    if not changed:
        break
assert not changed
allsyms = set().union(*({str(x) for x in v.free_symbols} for v in imposed.values()))
assert not [s for s in allsyms if s in imposed]
zero = {sp.Symbol(s): 0 for s in allsyms if s.startswith('k')}
bad = [(s, v) for s, v in imposed.items()
       if sp.simplify(sp.expand(v.subs(zero)) - shear.get(s, sp.Integer(0))) != 0]
print(f"{len(imposed)} impositions resolved; free knobs: "
      f"{len([s for s in allsyms if s.startswith('k')])}")
print("PASS: cascade at free->0 == shear exactly" if not bad else f"FAIL: {bad}")
sys.exit(1 if bad else 0)
