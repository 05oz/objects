"""The smallest error each code cannot see.

A quantum code protects information by spreading it out. An error touching a few qubits is
detected and undone; an error touching enough of the right ones is not, and the information is
silently gone. The number of qubits in the smallest such error is the code's distance, and
determining it is the hard part.

Five bivariate bicycle codes are drawn here, each on the torus its own parity checks define.
The lattice is not assumed: it is recovered by testing which factorisation of the check count
makes every row of H a cyclic shift of the first, which gives 6 x 6, 15 x 3, 9 x 6, 12 x 6 and
12 x 12 for the five. Every qubit is a point on that torus.

Lit on each are the qubits of a minimum-weight logical operator, taken from the witness the
certificates ship. That set is an error the code cannot detect: it commutes with every check,
so no measurement fires, and it changes the encoded state. Nothing smaller does. Those witnesses
are what establish the distance from above, and the matching lower bounds are the LRAT proofs in
the same deposit.

The lit sets have 6, 10, 10, 12 and 18 qubits, on tori of 72, 90, 108, 144 and 288. So the
better code is not the one with the larger lattice but the one whose smallest invisible error is
larger, and that is the thing to look at: the constellation grows while remaining, in every case,
a vanishing fraction of the surface.

They also do not clump. A minimum-weight logical winds; it cannot be squeezed into one region,
because anything local is exactly what the checks are built to catch.

doi:10.5281/zenodo.21831995
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
QDIR = "/Users/kirt/Documents/reserch math/certify-repo/qec-certificates"
CODES = ["bb72", "bb90", "bb108", "bb144", "bb288"]
CELL = 900
W, H = 3 * CELL, int(1.82 * CELL)

def load_H(path):
    t = open(path).read().split()
    m, n = int(t[0]), int(t[1])
    return m, n, [[int(c) for c in t[2 + i]] for i in range(m)]

def recover_torus(H_, n):
    """the lattice, derived: the unique l x m for which every row is row 0 shifted"""
    half = n // 2
    A = [r[:half] for r in H_]
    sup0 = [j for j, v in enumerate(A[0]) if v]
    for l in range(2, half + 1):
        if half % l:
            continue
        m_ = half // l
        if all(sorted(((j // m_ + i // m_) % l) * m_ + (j % m_ + i) % m_ for j in sup0)
               == sorted(j for j, v in enumerate(A[i]) if v) for i in range(1, len(A))):
            return l, m_
    return None

def project(p, ax=1.16, ay=0.36, d=3.6, s=1.0):
    x, y, z = p
    y, z = y * math.cos(ax) - z * math.sin(ax), y * math.sin(ax) + z * math.cos(ax)
    x, z = x * math.cos(ay) + z * math.sin(ay), -x * math.sin(ay) + z * math.cos(ay)
    f = d / (d - z)
    return x * f * s, y * f * s, z

def panel(code, ox, oy):
    mH, n, HX = load_H(os.path.join(QDIR, code, "HX.txt"))
    w = json.load(open(os.path.join(QDIR, code, "witness_X.json")))
    sup = {j for j, v in enumerate(w["x"]) if v}   # |x| is the stated weight
    L, M = recover_torus(HX, n)
    half = n // 2
    r_min = max(0.20, min(0.52, 0.52 * (M / L) ** 0.55))

    def pos(j):
        side, k = divmod(j, half)
        i, c = divmod(k, M)
        u = 2 * math.pi * i / L
        v = 2 * math.pi * c / M + (0.30 if side else -0.30) * (2 * math.pi / M)
        r = 1.0 + r_min * math.cos(v)
        return (r * math.cos(u), r * math.sin(u), r_min * math.sin(v))

    S = CELL * 0.295
    def scr(j):
        X, Y, Z = project(pos(j), s=S)
        return (ox + CELL / 2 + X, oy + CELL / 2 + Y, Z)
    P = [scr(j) for j in range(n)]
    zs = [p[2] for p in P]
    lo, hi = min(zs), max(zs)
    dep = lambda z: (z - lo) / (hi - lo or 1)

    # the torus itself: each qubit joined to its lattice neighbours, so the surface reads.
    # short segments only -- a step that wraps the torus would cut straight through it.
    cand = []
    for j in range(n):
        side, k = divmod(j, half)
        i, c = divmod(k, M)
        for di, dc in ((1, 0), (0, 1)):
            j2 = side * half + ((i + di) % L) * M + (c + dc) % M
            a, b = P[j], P[j2]
            cand.append((math.hypot(a[0] - b[0], a[1] - b[1]), a, b))
    # a step that wraps the torus cuts straight through it; drop those by taking the
    # threshold from this lattice's own median step rather than a fixed length
    med = sorted(d for d, _, _ in cand)[len(cand) // 2]
    seg = [((a[2] + b[2]) / 2, a, b) for d, a, b in cand if d < med * 2.6]
    seg.sort(key=lambda t: t[0])

    o = []
    for z, a, b in seg:
        t = dep(z)
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#2f4a5e" '
                 'stroke-width="%.2f" opacity="%.3f"/>'
                 % (a[0], a[1], b[0], b[1], 0.55 + 0.85 * t, 0.20 + 0.40 * t))

    order = sorted(range(n), key=lambda j: P[j][2])
    for j in order:
        x, y, z = P[j]
        t = dep(z)
        if j in sup:
            o.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="#ffd98a" opacity="%.2f"/>'
                     % (x, y, 15 + 9 * t, 0.10 + 0.09 * t))
            o.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="#fff3d2" opacity="%.2f"/>'
                     % (x, y, 4.2 + 3.0 * t, 0.80 + 0.20 * t))
        else:
            o.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="#6f8ea6" opacity="%.2f"/>'
                     % (x, y, 1.9 + 1.7 * t, 0.32 + 0.44 * t))
    return "".join(o), (code, n, L, M, w["weight"], len(sup))

if __name__ == "__main__":
    body, meta = [], []
    place = [(0, 0), (1, 0), (2, 0), (0.5, 0.86), (1.5, 0.86)]
    for k, c in enumerate(CODES):
        cc, rr = place[k]
        s, m = panel(c, cc * CELL, rr * CELL)
        body.append(s); meta.append(m)
    svg = ('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
           % (W, H, W, H)
           + '<rect width="%d" height="%d" fill="#06070c"/>' % (W, H)
           + "".join(body) + "</svg>")
    open(os.path.join(HERE, "smallest.svg"), "w").write(svg)
    for code, n, L, M, wt, s in meta:
        print(f"  {code:<7} {n:>3} qubits on {L}x{M:<6} minimum-weight logical: {wt} "
              f"({100*s/n:.1f}% of the surface)")
    print("  bytes:", len(svg))
