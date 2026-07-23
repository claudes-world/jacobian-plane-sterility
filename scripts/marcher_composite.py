"""Bucket marcher for the COMPOSITE-JEWEL cell (assembly section 5, route 1).

Cell: (N, r) = (3, 6), d = gcd(N-1, r) = 2, xi = x^3 Y, cone polynomial
E(xi) = 1 + t1*xi + xi^2 at the FERTILE locus t1 = sqrt(3) (P_1(t) ~ t^2-3;
the t = -sqrt(3) twin is the image under x -> -x). Chart data:

    c = gamma*x + c_min + tails,   c_min = Y + sqrt(3) x^3 Y^2 + x^6 Y^3
    J(w, c) = lam * x^-(1+k)       (pole depth k; default k = 1)

Grading: wt(x^a Y^b) = a - 3b. All three c_min monomials sit at wt -3;
J(., c_min) is graded of degree -1 (bracket adds wts then +2), so the tower
at bucket beta solves unknowns of weight beta+1 against equations of weight
beta. gamma*x sits at wt +1. Admissible tail slots: wt > -3 (strictly above
the cone edge), 0 <= b <= 3 = N. The lam source lands at bucket -(1+k).

Design rules inherited verbatim from marcher2.py (the (3,1,3) marcher):
impose-or-halt; window-termination assertion; pinned (b,a)-lex shield rule;
canonical cokernel normalization; E1-free assertion on conditions. Ground
field Q(sqrt(3)) via sympy's canonical sqrt(3).

New here vs marcher2: PRE-BUCKETS below the lam bucket are marched with zero
source and their tower KERNEL dimension recorded — kernel 0 there certifies
no deeper-pole gauge freedom (w has pole exactly x^-k in-window). The germ
bucket's solved layer is printed for comparison with the closed form
w_germ = -lam (x^-1 + sqrt(3) x^2 Y) from probe_composite_cone.py (k=1).

Run examples:
  python3 marcher_composite.py --minimal --no-e1 --max-bucket 24 --bmax 20 --amax 90 \
      --json out/mc_min.json          # primary clean minimal ladder
  python3 marcher_composite.py --minimal --max-bucket 18 --bmax 20 --amax 90 \
      --json out/mc_min_e1.json       # E1*c_min kernel-seed robustness run
  python3 marcher_composite.py --minimal --no-e1 --k 4 --max-bucket 12 ...  # k=4 germ family
"""
import argparse, json, time
import sympy as sp

lam, gam, E1 = sp.symbols('lam gam E1')
SQ3 = sp.sqrt(3)
CMIN = {(0, 1): sp.Integer(1), (3, 2): SQ3, (6, 3): sp.Integer(1)}
EXCLUDED = {(0, 0), (0, 1), (3, 2), (6, 3), (1, 0)}   # const, c_min, gamma slot


def bracket(d1, d2):
    """Poisson/Jacobian bracket of two monomial dicts {(a,b): coeff}."""
    out = {}
    for (a, b), u in d1.items():
        for (c, d), v in d2.items():
            k = a * d - b * c
            if k:
                key = (a + c - 1, b + d - 1)
                out[key] = out.get(key, 0) + k * u * v
    return out


def selftest():
    """Certify bracket() against sympy differentiation, trinomial c_min included."""
    x, Y = sp.symbols('x Y')
    f = 2 * x ** 3 * Y ** 2 - 5 * x ** -1 + 7 * x * Y ** 5
    g = Y + SQ3 * x ** 3 * Y ** 2 + x ** 6 * Y ** 3 - 4 * x ** 2
    direct = sp.expand(sp.diff(f, x) * sp.diff(g, Y) - sp.diff(f, Y) * sp.diff(g, x))
    fd = {(3, 2): sp.Integer(2), (-1, 0): sp.Integer(-5), (1, 5): sp.Integer(7)}
    gd = {(0, 1): sp.Integer(1), (3, 2): SQ3, (6, 3): sp.Integer(1),
          (2, 0): sp.Integer(-4)}
    viab = sum(co * x ** a * Y ** b for (a, b), co in bracket(fd, gd).items())
    assert sp.expand(direct - viab) == 0, "bracket formula selftest FAILED"


def wt(m):
    return m[0] - 3 * m[1]


def knob_slots(max_wt, min_wt=-2):
    """Chart tail slots (a,b), 0 <= b <= 3, wt in [min_wt, max_wt], excluding
    c_min/gamma/constant, in (b, a)-lex order. min_wt=-2 is the full wt > -3
    admissible lattice (strictly above the cone edge)."""
    slots = []
    for b in range(0, 4):
        for a in range(0, max_wt + 3 * b + 1):
            m = (a, b)
            if m not in EXCLUDED and min_wt <= wt(m) <= max_wt and a + b >= 2:
                slots.append(m)
    return slots


def sname(m):
    return f'k{m[0]}_{m[1]}'


class WindowError(RuntimeError):
    pass


def march(k, max_bucket, max_knob_wt, bmax, amax, minimal=False,
          inject_e1=True, verbose=True):
    selftest()
    lam_bucket = -(1 + k)
    knobs = {}
    cd = dict(CMIN); cd[(1, 0)] = gam
    if not minimal:
        for m in knob_slots(max_knob_wt):
            knobs[m] = sp.Symbol(sname(m))
            cd[m] = knobs[m]
    tails = {m: v for m, v in cd.items() if m not in CMIN}
    w = dict()
    if inject_e1:
        for m, co in CMIN.items():          # E1 * c_min kernel seed
            w[m] = co * E1
    sub = {}
    gauge_syms = [E1] if inject_e1 else []
    records, shields_seen = [], []

    for beta in range(lam_bucket - 3, max_bucket + 1):
        t0 = time.time()
        src = {}
        for m, v in bracket(w, tails).items():
            if wt(m) == beta:
                sv = sp.expand(v.subs(sub))
                if sv != 0:
                    src[m] = sv
        if beta == lam_bucket:
            src[(-1 - k, 0)] = src.get((-1 - k, 0), 0) - lam
        if beta < lam_bucket and src:
            raise RuntimeError(f"pre-bucket {beta}: unexpected source {src}")
        mons = [(a, b) for b in range(0, bmax + 1) for a in range(-2 - k, amax + 1)
                if wt((a, b)) == beta + 1 and (a, b) != (0, 0)]
        us = {m: sp.Symbol(f'u{beta}_{m[0]}_{m[1]}') for m in mons}
        Jc = bracket({m: us[m] for m in mons}, CMIN)
        eqm = sorted(set(Jc) | set(src))
        M = sp.Matrix([[sp.expand(Jc.get(m, 0)).coeff(us[mm]) for mm in mons]
                       for m in eqm])
        rec = {"bucket": beta, "n_mons": len(mons), "n_eqm": len(eqm),
               "cokernel_dim": 0, "kernel_dim": None, "nonzero_conds": 0,
               "shield": None, "shield_value": None, "core_canonical": None,
               "cond_canonical": None, "gauge_leak": None, "terminal": False}
        if beta <= lam_bucket and len(mons):
            # record kernel SUPPORT at pre/germ buckets: certifies exactly
            # which gauge directions exist (expect only multiples of c_min)
            rec["kernel_support"] = [
                sorted([list(m) for m, uv in zip(mons, list(kv)) if uv != 0])
                for kv in M.nullspace()]

        # ---- one exact rref of [M | I] gives cokernel functionals AND the
        # solution template; M is knob-free (Q(sqrt3)), so no symbolic
        # pivoting ever happens (the EX-domain zero-test hazard is gone) ----
        neq, nun = len(eqm), len(mons)
        if neq:
            Raug, piv_cols = sp.Matrix.hstack(M, sp.eye(neq)).rref()
            sol_rows = [(i, pc) for i, pc in enumerate(piv_cols) if pc < nun]
            cok_rows = [i for i, pc in enumerate(piv_cols) if pc >= nun]
        else:
            Raug, sol_rows, cok_rows = None, [], []
        rec["kernel_dim"] = nun - len(sol_rows)

        def canon(e):
            return sp.cancel(sp.together(sp.expand(e)))

        def sub_add(s, v):
            """Keep sub RESOLVED (idempotent): new values are closed under
            sub, and every stored value learns the new pair. Without this,
            chained impositions leave stale knob references in sources."""
            v = canon(v.subs(sub))
            for kk in list(sub):
                sub[kk] = sub[kk].subs(s, v)
            sub[s] = v
            return v

        for irow in cok_rows:
            vv = Raug[irow, nun:].T                 # cokernel functional (rref-normalized)
            rec["cokernel_dim"] += 1
            rhs_now = sp.Matrix([canon(-sp.expand(src.get(m, 0)).subs(sub))
                                 for m in eqm])
            cond_c = canon((vv.T * rhs_now)[0])     # canonical rational form
            cond_den = sp.denom(cond_c)
            cond = sp.expand(sp.numer(cond_c))      # impose on the numerator
            if cond == 0:
                continue
            rec["nonzero_conds"] += 1
            leak = [str(gsym) for gsym in gauge_syms if cond.has(gsym)]
            rec["gauge_leak"] = leak or None
            piv = next(j for j in range(neq) if vv[j] != 0)
            rec["pivot_monomial"] = list(eqm[piv])
            if cond_den != 1:
                rec.setdefault("cond_denominators", []).append(str(cond_den))
            # ---- shield search (pinned rule) ----
            # minimal mode adds the raw shield delta to the RATIONAL form
            # cond_c (numerator-solving happens at the end); full mode works
            # with the numerator polynomial directly
            target_slot, target_sym, contrib_cond, target_dsrc = None, None, cond_c, None
            if minimal:
                for slot in knob_slots(beta + 1):
                    if wt(slot) != beta + 1 or slot in cd:
                        continue
                    sigma = sp.Symbol(sname(slot))
                    dsrc = bracket(w, {slot: sigma})
                    parts = []
                    for m, v in dsrc.items():
                        if wt(m) != beta:
                            continue
                        vs = sp.expand(v.subs(sub))
                        if vs == 0:
                            continue
                        if m not in eqm:      # P3: no silent truncation, ever
                            raise RuntimeError(
                                f"bucket {beta}: shield contribution at {m} "
                                f"outside equation support — enlarge windows")
                        parts.append(-vv[eqm.index(m)] * vs)
                    delta = sp.expand(sum(parts))
                    if delta.has(sigma) and sp.degree(sp.Poly(delta, sigma)) == 1:
                        target_slot, target_sym, target_dsrc = slot, sigma, dsrc
                        contrib_cond = canon(cond_c + delta)
                        break
            else:
                for slot in [m for m in knobs
                             if wt(m) == beta + 1 and knobs[m] not in sub]:
                    ksym = knobs[slot]
                    if cond.has(ksym) and sp.degree(sp.Poly(cond, ksym)) == 1:
                        target_slot, target_sym = slot, ksym
                        break
                if target_sym is None:
                    # No FRESH shield: the condition constrains OLDER chart
                    # data. Forced-imposition rule: pick the (b,a)-lex first
                    # unfixed knob appearing linearly whose linear coefficient
                    # is a manifestly nonzero monomial in lam, gam (never a
                    # branch split, since lam != 0 and gam != 0 are standing
                    # hypotheses). Impose-or-halt otherwise.
                    for slot in sorted(
                            [m for m in knobs if knobs[m] not in sub
                             and cond.has(knobs[m])],
                            key=lambda m: (m[1], m[0])):
                        ksym = knobs[slot]
                        pc = sp.Poly(cond, ksym)
                        if sp.degree(pc) != 1:
                            continue
                        coeff = sp.expand(sp.radsimp(pc.coeff_monomial(ksym)))
                        if coeff == 0 or (coeff.free_symbols - {lam, gam}) \
                                or len(sp.Add.make_args(coeff)) != 1:
                            continue
                        target_slot, target_sym = slot, ksym
                        if verbose:
                            print(f'   [bucket {beta}: forced imposition on '
                                  f'older knob {sname(slot)} '
                                  f'(coeff {coeff})]', flush=True)
                        break
            fresh_syms = ([] if minimal else
                          [knobs[m] for m in knobs
                           if wt(m) == beta + 1 and knobs[m] not in sub])
            core = sp.factor(sp.radsimp(cond.subs({f: 0 for f in fresh_syms})))
            rec["core_canonical"] = str(core)
            rec["cond_canonical"] = str(cond)
            if verbose:
                print(f'bucket {beta}: shield={target_slot} '
                      f'core(fresh->0) = {core}'
                      + (f'  GAUGE LEAK {leak}' if leak else ''), flush=True)
            if target_sym is None:
                rec["terminal"] = True
                rec["cokernel_functional"] = {str(m): str(vv[i])
                                              for i, m in enumerate(eqm)
                                              if vv[i] != 0}
                rec["shields_tried"] = [list(s) for s in knob_slots(beta + 1)
                                        if wt(s) == beta + 1 and s not in cd] \
                    if minimal else sorted(sname(m) for m in knobs
                                           if wt(m) == beta + 1)
                rec["cond_knob_deps"] = sorted(
                    str(s) for s in cond.free_symbols
                    if str(s).startswith('k'))
                rec["seconds"] = round(time.time() - t0, 2)
                records.append(rec)
                raise TerminalObstruction(beta, sp.factor(sp.radsimp(cond)), records)
            val = sp.solve(sp.expand(sp.numer(canon(contrib_cond))),
                           target_sym)[0]
            shields_seen.append((beta, target_slot))
            if minimal:
                val = canon(val.subs(sub))
                cd[target_slot] = val
                tails[target_slot] = val
                for m, v in target_dsrc.items():
                    if wt(m) == beta:
                        src[m] = sp.expand(src.get(m, 0)
                                           + v.subs(target_sym, val))
            else:
                val = sub_add(target_sym, val)
            rec["shield"], rec["shield_value"] = sname(target_slot), str(val)
            rec.setdefault("shields_all", []).append(
                [sname(target_slot), str(val)])
        # ---- solve tower via the template: x_pivot = (T b)_row, params 0;
        # consistency is exactly the cokernel conditions imposed above, so a
        # leftover inconsistency is impossible by construction (each cond
        # numerator was driven to 0 or the march halted) ----
        rhs2 = sp.Matrix([canon(-sp.expand(src.get(m, 0)).subs(sub))
                          for m in eqm])
        residual_check = [i for i in cok_rows
                          if sp.expand(sp.numer(canon(
                              (Raug[i, nun:] * rhs2)[0]))) != 0]
        if residual_check:
            raise RuntimeError(f"bucket {beta}: cokernel residual survived "
                               f"imposition at rows {residual_check} "
                               f"(silent-drop guard)")
        edge_hits = []
        layer = {}
        for irow, pc in sol_rows:
            m = mons[pc]
            val = canon((Raug[irow, nun:] * rhs2)[0])
            if val != 0:
                if m[1] >= bmax - 1 or m[0] >= amax - 1:
                    edge_hits.append(m)
                w[m] = canon(w.get(m, 0) + val)
                layer[m] = val
        if beta < lam_bucket and (layer or rec["kernel_dim"]):
            print(f'  PRE-BUCKET {beta}: kernel_dim={rec["kernel_dim"]} '
                  f'layer={layer}  (deep-pole check)', flush=True)
        if beta == lam_bucket and verbose:
            print(f'  GERM layer (bucket {beta}): '
                  f'{ {m: str(v) for m, v in layer.items()} }', flush=True)
        if edge_hits:
            raise WindowError(f"bucket {beta}: nonzero tower coefficients on "
                              f"window edge {edge_hits}; enlarge --bmax/--amax")
        rec["seconds"] = round(time.time() - t0, 2)
        records.append(rec)
        if verbose:
            print(f'  [bucket {beta}: {len(mons)} unk, {len(eqm)} eq, '
                  f'coker {rec["cokernel_dim"]}, ker {rec["kernel_dim"]}, '
                  f'{rec["seconds"]}s]', flush=True)
    return records, sub, w


class TerminalObstruction(RuntimeError):
    def __init__(self, bucket, core, records):
        super().__init__(f"TERMINAL obstruction at bucket {bucket}: {core}")
        self.bucket, self.core, self.records = bucket, core, records


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--k', type=int, default=1, help='pole depth (fertile: k=1 mod 3)')
    ap.add_argument('--max-bucket', type=int, default=12)
    ap.add_argument('--max-knob-wt', type=int, default=8)
    ap.add_argument('--bmax', type=int, default=16)
    ap.add_argument('--amax', type=int, default=60)
    ap.add_argument('--minimal', action='store_true')
    ap.add_argument('--no-e1', action='store_true',
                    help='skip the E1*c_min kernel injection (clean minimal ladder)')
    ap.add_argument('--json', type=str, default=None)
    args = ap.parse_args()
    t0 = time.time()
    outcome = {"terminal": None}
    try:
        records, sub, w = march(args.k, args.max_bucket, args.max_knob_wt,
                                args.bmax, args.amax, minimal=args.minimal,
                                inject_e1=not args.no_e1)
    except TerminalObstruction as e:
        records = e.records
        outcome = {"terminal": {"bucket": e.bucket, "core": str(e.core)}}
        print(f"\n*** {e} ***")
    meta = {"cell": "(N,r)=(3,6) t1=sqrt(3) composite-jewel", "k": args.k,
            "mode": "minimal" if args.minimal else "full",
            "e1": not args.no_e1, "max_bucket": args.max_bucket,
            "bmax": args.bmax, "amax": args.amax,
            "total_seconds": round(time.time() - t0, 1), **outcome}
    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"meta": meta, "buckets": records}, fh, indent=1)
        print(f"wrote {args.json}")
    print(f"done in {meta['total_seconds']}s")


if __name__ == '__main__':
    main()
