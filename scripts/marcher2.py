"""Clean-room bucket marcher for the (3,1,3) cell gamma-ladder.

Independent implementation (not a transcription of C3): re-derived from the
bracket formula J(x^a Y^b, x^c Y^d) = (ad - bc) x^{a+c-1} Y^{b+d-1}, which is
machine-checked in selftest(). Design rules:

  * impose-or-halt: an unresolved obstruction STOPS the march (P2);
  * window-termination assertion: solved tower coefficients on the window
    edge must vanish, else WindowError (P4);
  * pinned shield rule: among fresh knobs at weight beta+1 appearing linearly
    in the condition, take the (b, a)-lexicographically smallest slot;
  * canonical cokernel normalization: scale each obstruction functional so its
    entry on the lex-lowest equation monomial in its support is +1 — printed
    cores are then presentation-independent;
  * conditions are asserted E1-free (and E2-free with --e2); violations recorded.

Chart data: c = c_min + gamma*x + tails, c_min = Y + x^3 Y^3, weight
v(x^a Y^b) = 2a - 3b, pair equation J(w, c) = lam * x^-2 (k = 1).

Modes:
  full:    all knob slots of weight 0..max_knob_wt symbolic (default 12)
  minimal: no optional knobs; shield slots introduced lazily at obstructed
           buckets and solved to rational multiples of gamma powers.

Run examples:
  python3 marcher2.py --max-bucket 11 --json out/m2_repro.json          # cross-check C3
  python3 marcher2.py --max-bucket 16 --max-knob-wt 17 --bmax 20 --amax 35 --json out/m2_L4.json
  python3 marcher2.py --minimal --max-bucket 46 --bmax 40 --amax 80 --json out/m2_min.json
"""
import argparse, json, time, sys
import sympy as sp

lam, gam, E1, E2 = sp.symbols('lam gam E1 E2')
CMIN = {(0, 1): sp.Integer(1), (3, 3): sp.Integer(1)}
EXCLUDED = {(0, 0), (0, 1), (3, 3), (1, 0)}   # const, c_min, gamma slot


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
    """Certify bracket() against sympy differentiation on a random-ish sample."""
    x, Y = sp.symbols('x Y')
    f = 2 * x ** 3 * Y ** 2 - 5 * x ** -1 + 7 * x * Y ** 5
    g = Y + x ** 3 * Y ** 3 - 4 * x ** 2
    direct = sp.expand(sp.diff(f, x) * sp.diff(g, Y) - sp.diff(f, Y) * sp.diff(g, x))
    fd = {(3, 2): sp.Integer(2), (-1, 0): sp.Integer(-5), (1, 5): sp.Integer(7)}
    gd = {(0, 1): sp.Integer(1), (3, 3): sp.Integer(1), (2, 0): sp.Integer(-4)}
    viab = sum(co * x ** a * Y ** b for (a, b), co in bracket(fd, gd).items())
    assert sp.expand(direct - viab) == 0, "bracket formula selftest FAILED"


def wt(m):
    return 2 * m[0] - 3 * m[1]


def knob_slots(max_wt, min_wt=0):
    """All chart tail slots (a,b), 0 <= b <= 3, weight in [min_wt, max_wt],
    excluding c_min/gamma/constant, in (b, a)-lex order (the canonical shield
    order). min_wt=-2 admits the full v > -r admissible lattice (a+b >= 2 is
    automatic outside EXCLUDED)."""
    slots = []
    for b in range(0, 4):
        for a in range(0, (max_wt + 3 * b) // 2 + 1):
            m = (a, b)
            if m not in EXCLUDED and min_wt <= wt(m) <= max_wt and a + b >= 2:
                slots.append(m)
    return slots


def sname(m):
    return f'k{m[0]}_{m[1]}'


class WindowError(RuntimeError):
    pass


def march(max_bucket, max_knob_wt, bmax, amax, minimal=False, with_e2=False,
          inject_e1=True, min_knob_wt=0, verbose=True):
    selftest()
    knobs = {}                      # slot -> symbol (full mode)
    cd = dict(CMIN); cd[(1, 0)] = gam
    if not minimal:
        for m in knob_slots(max_knob_wt, min_knob_wt):
            knobs[m] = sp.Symbol(sname(m))
            cd[m] = knobs[m]
    tails = {m: v for m, v in cd.items() if m not in CMIN}
    w = {(0, 1): E1, (3, 3): E1} if inject_e1 else {}
    if with_e2:
        # E2 * c_min^2 probe. NOTE (2026-07-22 finding): unlike E1*c_min, this
        # is NOT an absorbable injection — a real E2-obstruction fires at
        # bucket -1 (shield (3,2)). The honest kernel elements are g(c) with
        # J(g(c), c) = 0 identically (no march needed); c_min-power seedings
        # test absorbability of the tail-deficit, which holds for E1, not E2.
        # (Seeding E2*c^2 instead would enforce a WRONG equation here: the
        # marcher never brackets seeded w against c_min.)
        for m, co in {(0, 2): 1, (3, 4): 2, (6, 6): 1}.items():
            w[m] = w.get(m, 0) + co * E2
    sub = {}
    gauge_syms = [E1] + ([E2] if with_e2 else [])
    records, shields_seen = [], []

    for beta in range(-4, max_bucket + 1):
        t0 = time.time()
        src = {}
        for m, v in bracket(w, tails).items():
            if wt(m) == beta:
                sv = sp.expand(v.subs(sub))
                if sv != 0:
                    src[m] = sv
        if beta == -4:
            src[(-2, 0)] = src.get((-2, 0), 0) - lam
        mons = [(a, b) for b in range(0, bmax + 1) for a in range(-3, amax + 1)
                if wt((a, b)) == beta + 2 and (a, b) != (0, 0)]
        us = {m: sp.Symbol(f'u{beta}_{m[0]}_{m[1]}') for m in mons}
        Jc = bracket({m: us[m] for m in mons}, CMIN)
        eqm = sorted(set(Jc) | set(src))
        M = sp.Matrix([[sp.expand(Jc.get(m, 0)).coeff(us[mm]) for mm in mons]
                       for m in eqm])
        rec = {"bucket": beta, "n_mons": len(mons), "n_eqm": len(eqm),
               "cokernel_dim": 0, "nonzero_conds": 0, "shield": None,
               "shield_value": None, "core_canonical": None, "cond_canonical": None,
               "gauge_leak": None, "terminal": False}

        for vv in (M.T.nullspace() if len(eqm) else []):
            rec["cokernel_dim"] += 1
            piv = next(i for i in range(len(eqm)) if vv[i] != 0)
            vv = vv / vv[piv]                       # canonical normalization
            rhs_now = sp.Matrix([-sp.expand(sp.expand(src.get(m, 0)).subs(sub))
                                 for m in eqm])
            cond = sp.expand((vv.T * rhs_now)[0])
            if cond == 0:
                continue
            rec["nonzero_conds"] += 1
            leak = [str(gsym) for gsym in gauge_syms if cond.has(gsym)]
            rec["gauge_leak"] = leak or None
            rec["pivot_monomial"] = list(eqm[piv])
            # ---- shield search (pinned rule) ----
            target_slot, target_sym, contrib_cond, target_dsrc = None, None, cond, None
            if minimal:
                for slot in knob_slots(beta + 1, min_wt=-2):
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
                        contrib_cond = sp.expand(cond + delta)
                        break
            else:
                for slot in [m for m in knobs
                             if wt(m) == beta + 1 and knobs[m] not in sub]:
                    ksym = knobs[slot]
                    if cond.has(ksym) and sp.degree(sp.Poly(cond, ksym)) == 1:
                        target_slot, target_sym = slot, ksym
                        break
            fresh_syms = ([] if minimal else
                          [knobs[m] for m in knobs
                           if wt(m) == beta + 1 and knobs[m] not in sub])
            core = sp.factor(cond.subs({f: 0 for f in fresh_syms}))
            rec["core_canonical"] = str(core)
            rec["cond_canonical"] = str(cond)
            if verbose:
                print(f'bucket {beta}: shield={target_slot} '
                      f'core(fresh->0) = {core}'
                      + (f'  GAUGE LEAK {leak}' if leak else ''))
            if target_sym is None:
                rec["terminal"] = True
                records.append(rec)
                raise RuntimeError(
                    f"TERMINAL obstruction at bucket {beta}: {sp.factor(cond)}")
            val = sp.expand(sp.solve(contrib_cond, target_sym)[0].subs(sub))
            rec["shield"], rec["shield_value"] = sname(target_slot), str(val)
            shields_seen.append((beta, target_slot))
            if minimal:
                cd[target_slot] = val               # concrete tail coefficient
                tails[target_slot] = val
                # fold the shield's bracket contribution into THIS bucket's
                # source (in full mode the knob symbol was in src from the
                # start; lazy introduction must add it explicitly)
                for m, v in target_dsrc.items():
                    if wt(m) == beta:
                        src[m] = sp.expand(src.get(m, 0)
                                           + v.subs(target_sym, val))
            else:
                sub[target_sym] = val
        # ---- solve tower (impose-or-halt semantics guaranteed above) ----
        rhs2 = sp.Matrix([-sp.expand(sp.expand(src.get(m, 0)).subs(sub))
                          for m in eqm])
        try:
            sol, params = M.gauss_jordan_solve(rhs2)
        except ValueError as e:
            raise RuntimeError(f"bucket {beta}: tower inconsistent after imposition "
                               f"(silent-drop guard): {e}") from e
        sol = sol.subs({p: 0 for p in params})
        edge_hits = []
        for m, val in zip(mons, list(sol)):
            val = sp.expand(val.subs(sub))
            if val != 0:
                if m[1] >= bmax - 1 or m[0] >= amax - 1:
                    edge_hits.append(m)
                w[m] = sp.expand(w.get(m, 0) + val)
        if edge_hits:
            raise WindowError(f"bucket {beta}: nonzero tower coefficients on "
                              f"window edge {edge_hits}; enlarge --bmax/--amax")
        rec["seconds"] = round(time.time() - t0, 2)
        records.append(rec)
        if verbose:
            print(f'  [bucket {beta}: {len(mons)} unk, {len(eqm)} eq, '
                  f'coker {rec["cokernel_dim"]}, {rec["seconds"]}s]')
    return records, sub, w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--max-bucket', type=int, default=11)
    ap.add_argument('--max-knob-wt', type=int, default=12)
    ap.add_argument('--bmax', type=int, default=15)
    ap.add_argument('--amax', type=int, default=25)
    ap.add_argument('--minimal', action='store_true')
    ap.add_argument('--e2', action='store_true')
    ap.add_argument('--no-e1', action='store_true',
                    help='skip the E1*c_min kernel injection (clean minimal ladder)')
    ap.add_argument('--min-knob-wt', type=int, default=0,
                    help='lower knob-weight bound; -2 = full v > -r admissible lattice')
    ap.add_argument('--json', type=str, default=None)
    args = ap.parse_args()
    t0 = time.time()
    records, sub, w = march(args.max_bucket, args.max_knob_wt, args.bmax,
                            args.amax, minimal=args.minimal, with_e2=args.e2,
                            inject_e1=not args.no_e1, min_knob_wt=args.min_knob_wt)
    meta = {"mode": "minimal" if args.minimal else "full",
            "max_bucket": args.max_bucket, "max_knob_wt": args.max_knob_wt,
            "bmax": args.bmax, "amax": args.amax, "e2": args.e2,
            "total_seconds": round(time.time() - t0, 1)}
    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"meta": meta, "buckets": records}, fh, indent=1)
        print(f"wrote {args.json}")
    print(f"done in {meta['total_seconds']}s")


if __name__ == '__main__':
    main()
