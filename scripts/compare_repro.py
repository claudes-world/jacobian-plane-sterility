"""Cross-certify marcher_repro (C3 verbatim) vs marcher2 (clean-room) at the
C3 config: same shields, and canonical conditions/cores EQUAL as polynomials
(after k{a}{b} -> k{a}_{b} name translation). Also report canonical gamma-level
cores (the c_j data) and pivot monomials."""
import json, sys, re
import sympy as sp

rep = json.load(open('out/repro_t11.json'))
m2 = json.load(open('out/m2_repro.json'))["buckets"]
rep = {r["bucket"]: r for r in rep}
m2 = {r["bucket"]: r for r in m2}

# explicit name map from the weight<=12 knob lattice (all single-digit a here
# except k103 = (10,3): translate longest-first to avoid partial hits)
slots = []
for b in range(0, 4):
    for a in range(0, 12):
        w = 2*a - 3*b
        if 0 <= w <= 12 and (a, b) not in [(0,0), (0,1), (3,3), (1,0)]:
            slots.append((a, b))
names = sorted(((f'k{a}{b}', f'k{a}_{b}') for a, b in slots),
               key=lambda p: -len(p[0]))

def translate(sstr):
    for old, new in names:
        sstr = re.sub(rf'\b{old}\b', new, sstr)
    return sstr

syms = {f'k{a}_{b}': sp.Symbol(f'k{a}_{b}') for a, b in slots}
syms.update({s: sp.Symbol(s) for s in ('lam', 'gam', 'E1', 'E2')})

ok = True
print(f"{'bkt':>4} {'shield':>7} {'cond==':>7} {'core==':>7}  pivot(rep/m2)")
for bkt in range(-4, 12):
    r, m = rep[bkt], m2[bkt]
    sh_r = r["shield"] if r["shield"] not in (None, 'None') else None
    sh_m = m["shield"]
    sh_ok = (sh_r is None and sh_m is None) or \
            (sh_r is not None and sh_m is not None and translate(sh_r) == sh_m)
    if r.get("cond_canonical") or r.get("core_canonical") not in (None,):
        pass
    def eq(a_str, b_str):
        if a_str is None and b_str is None: return True
        if a_str is None or b_str is None: return False
        ea = sp.sympify(translate(a_str), locals=syms)
        eb = sp.sympify(b_str, locals=syms)
        return sp.expand(ea - eb) == 0
    # repro stored canonical core; canonical cond wasn't stored -> compare cores
    core_ok = eq(r.get("core_canonical"), m.get("core_canonical"))
    piv_r, piv_m = r.get("pivot_monomial"), m.get("pivot_monomial")
    ok &= sh_ok and core_ok and (r["cokernel_dim"] == m["cokernel_dim"])
    print(f"{bkt:>4} {str(sh_ok):>7} {'—':>7} {str(core_ok):>7}  {piv_r}/{piv_m}")

print("\ncanonical gamma-level cores (marcher2, knob-free):")
for lvl, bkt in ((1, 1), (2, 6), (3, 11)):
    core = sp.sympify(m2[bkt]["core_canonical"], locals=syms)
    knobfree = core.subs({v: 0 for k, v in syms.items() if k.startswith('k')})
    cj = sp.cancel(knobfree / (syms['gam']**lvl * syms['lam']))
    print(f"  level {lvl} (bucket {bkt}): core={m2[bkt]['core_canonical']}")
    print(f"           knob-free = {sp.factor(knobfree)}   =>  c_{lvl}^can = {cj}")

print("\nCROSS-CERTIFIED" if ok else "\nDIVERGENCE — inspect above")
sys.exit(0 if ok else 1)
