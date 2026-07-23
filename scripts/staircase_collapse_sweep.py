"""Audit-ask 3 (forward direction): every genuine c-staircase satisfies the
fiber-collapse hypothesis of FATEDGE-KILL.

Claim: vertices V_k = (i0 + sum_{s<=k} sigma_s h_s, sum_{s<=k} h_s) of a
staircase polygon (i0 >= 1, slopes sigma_1 > sigma_2 > ... > sigma_E > 0,
heights h_s >= 1, E >= 2 edges) satisfy omega = N*i - alpha*j > 0 for every
vertex except the top (alpha, N) — and hence so does every monomial on or
above the staircase hull.

Proof (recorded here, machine-swept below): with T = sum sigma_s h_s and
P_k the prefix sum, omega(V_k) = i0*(N - j_k) + (N*P_k - T*j_k); the first
term is >= 0 (positive for j_k < N), and P_k/j_k >= T/N because the prefix
weighted average of a strictly decreasing sequence strictly exceeds the
total average for intermediate k. Non-vertex monomials on/above the hull lie
weakly right of the vertex chain, hence strictly right of the chord.
Sweep: random staircases incl. rational slopes (integer lattice points only),
plus randomized above-hull tail monomials.
"""
import random
from fractions import Fraction

random.seed(723)
bad = []
trials = 0
for _ in range(4000):
    E = random.randint(2, 5)
    i0 = random.randint(1, 4)
    hs = [random.randint(1, 3) for _ in range(E)]
    # strictly decreasing positive rational slopes with denominator = h_s
    # (so vertices are lattice points): pick decreasing rationals
    while True:
        nums = sorted({Fraction(random.randint(1, 60), random.randint(1, 4))
                       for _ in range(E)}, reverse=True)
        if len(nums) >= E:
            slopes = nums[:E]
            break
    verts = [(Fraction(i0), Fraction(0))]
    ok_lattice = True
    for s in range(E):
        di = slopes[s] * hs[s]
        vi = verts[-1][0] + di
        vj = verts[-1][1] + hs[s]
        if vi.denominator != 1:
            ok_lattice = False
            break
        verts.append((vi, vj))
    if not ok_lattice:
        continue
    trials += 1
    alpha, N = verts[-1]
    # vertices strictly inside (omega > 0) except top
    for (vi, vj) in verts[:-1]:
        if N * vi - alpha * vj <= 0:
            bad.append(("vertex", verts, (vi, vj)))
    # random above-hull monomials: take a vertex, bump i by 1..4 (stays above)
    for _ in range(6):
        vi, vj = random.choice(verts)
        m = (vi + random.randint(1, 4), vj)
        if m[1] != N or m[0] != alpha:
            if N * m[0] - alpha * m[1] <= 0:
                bad.append(("tail", verts, m))
print(f"{trials} random staircases (incl. rational slopes): "
      f"{'ALL COLLAPSED (omega>0 off top)' if not bad else 'VIOLATIONS'}")
if bad:
    for b in bad[:5]:
        print("  BAD:", b)
    raise SystemExit(1)
