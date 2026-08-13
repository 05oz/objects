"""The bivariate bicycle family, each on its own torus.

Five codes whose minimum distances this program certified. Each one's torus is
recovered from its own shipped parity-check matrix by testing every factorisation
of the block size against the shift structure of all its rows — nothing is assumed
and nothing is read from a paper:

    [[72 ,12, 6]]   6 x 6      [[90 ,8,10]]  15 x 3      [[108,8,10]]   9 x 6
    [[144,12,12]]  12 x 6      [[288,12,18]] 12 x 12

The tori are not drawn to a chosen shape. The ring and tube radii are set by the
lattice itself, so a 15 x 3 code is a wide ring with a thin tube and a 6 x 6 is
almost round; the object's own proportions decide how it looks. Sizes are relative
to qubit count, so the family is seen at true scale against each other.

Teal is the X checks, rose the Z. Every edge is a geodesic in (u,v): a check's
connection lives on the lattice, so it is drawn along the surface, never through it.

doi:10.5281/zenodo.21831995
"""
import math, os

HERE = os.path.dirname(os.path.abspath(__file__))
QEC = "/Users/kirt/Documents/reserch math/certify-repo/qec-certificates"
W, H = 2600, 1000
SEG = 9

def load(p):
    t = open(p).read().split()
    m, n = int(t[0]), int(t[1])
    return m, n, [[int(c) for c in t[2 + i]] for i in range(m)]

def recover_torus(H_, n):
    """the lattice, derived: the unique l x m for which every row is row 0 shifted"""
    half = n // 2
    A = [r[:half] for r in H_]
    sup0 = [j for j, v in enumerate(A[0]) if v]
    for l in range(2, half + 1):
        if half % l: continue
        m_ = half // l
        if all(sorted(((j // m_ + i // m_) % l) * m_ + (j % m_ + i) % m_ for j in sup0)
               == sorted(j for j, v in enumerate(A[i]) if v) for i in range(1, len(A))):
            return l, m_
    return None

def hsv(h, s, v):
    h = (h % 360) / 60.0; i = int(h) % 6; f = h - int(h)
    p, q, t = v * (1 - s), v * (1 - s * f), v * (1 - s * (1 - f))
    r, g, b = ((v, t, p), (q, v, p), (p, v, t), (p, q, v), (t, p, v), (v, p, q))[i]
    return "#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255))

def wrap(d, P=2 * math.pi):
    while d > P / 2: d -= P
    while d < -P / 2: d += P
    return d

def render(name, cx, cy, scale):
    mX, n, HX = load(f"{QEC}/{name}/HX.txt")
    mZ, _, HZ = load(f"{QEC}/{name}/HZ.txt")
    L, M = recover_torus(HX, n)
    half = n // 2
    # the lattice sets the shape: tube-to-ring ratio follows m/l
    r_maj, r_min = 1.0, max(0.20, min(0.52, 0.52 * (M / L) ** 0.55))

    def uv(j):
        side, k = divmod(j, half)
        i, c = divmod(k, M)
        return (2 * math.pi * (i + (0.26 if side else -0.26)) / L,
                2 * math.pi * (c + (0.26 if side else -0.26)) / M)

    def surf(u, v):
        r = r_maj + r_min * math.cos(v)
        return (r * math.cos(u), r * math.sin(u), r_min * math.sin(v))

    def proj(p, ax=1.14, ay=0.34, d=3.6):
        x, y, z = p
        y, z = y * math.cos(ax) - z * math.sin(ax), y * math.sin(ax) + z * math.cos(ax)
        x, z = x * math.cos(ay) + z * math.sin(ay), -x * math.sin(ay) + z * math.cos(ay)
        f = d / (d - z)
        return cx + x * scale * f, cy + y * scale * f, z

    pts = {j: surf(*uv(j)) for j in range(n)}
    E = []
    for Hm, kind in ((HX, "x"), (HZ, "z")):
        for row in Hm:
            sup = [j for j, v in enumerate(row) if v]
            for a in range(len(sup)):
                for b in range(a + 1, len(sup)):
                    u1, v1 = uv(sup[a]); u2, v2 = uv(sup[b])
                    du, dv = wrap(u2 - u1), wrap(v2 - v1)
                    E.append(([surf(u1 + du * t / SEG, v1 + dv * t / SEG)
                               for t in range(SEG + 1)], kind))

    out = []
    drawn = sorted(((sum(q[2] for q in P) / len(P), P, k)
                    for P, k in ((tuple(proj(p) for p in path), kind) for path, kind in E)),
                   key=lambda t: t[0])
    zs = [d[0] for d in drawn]; lo, hi = min(zs), max(zs)
    for z, P, kind in drawn:
        t = (z - lo) / (hi - lo or 1)
        hue = 196 - 26 * t if kind == "x" else 336 + 14 * t
        out.append('<path d="M%.0f %.0f%s" fill="none" stroke="%s" stroke-width="%.2f" '
                   'opacity="%.3f" stroke-linecap="round"/>'
                   % (P[0][0], P[0][1], "".join("L%.0f %.0f" % (q[0], q[1]) for q in P[1:]),
                      hsv(hue, 0.26 + 0.44 * t, 0.30 + 0.66 * t),
                      0.25 + 1.5 * t ** 1.4 * (scale / 300),
                      0.05 + 0.68 * t ** 2.1))
    qs = sorted((proj(p), j) for j, p in pts.items())
    zq = [q[0][2] for q in qs]; qlo, qhi = min(zq), max(zq)
    for (X, Y, Z), j in qs:
        t = (Z - qlo) / (qhi - qlo or 1)
        out.append('<circle cx="%.0f" cy="%.0f" r="%.2f" fill="%s" opacity="%.3f"/>'
                   % (X, Y, (0.7 + 2.6 * t ** 1.5) * (scale / 300),
                      hsv(38 if j >= half else 190, 0.16 + 0.20 * t, 0.55 + 0.42 * t),
                      0.16 + 0.76 * t ** 1.4))
    return "".join(out), (L, M, n, mX + mZ)

if __name__ == "__main__":
    fam = ["bb72", "bb90", "bb108", "bb144", "bb288"]
    ns = {}
    for f in fam:
        _, n, _ = load(f"{QEC}/{f}/HX.txt"); ns[f] = n
    # true relative scale: linear in sqrt(qubits), so area tracks the code size
    base = 208.0 / math.sqrt(288)
    # place from the true sizes: each torus claims its own projected width, then a gap
    scales = {f: base * math.sqrt(ns[f]) for f in fam}
    widths = {f: scales[f] * 2.9 for f in fam}          # ring + tube + perspective
    gap = (W - sum(widths.values())) / (len(fam) + 1)
    xs, run = [], gap
    for f in fam:
        xs.append(run + widths[f] / 2); run += widths[f] + gap
    body = []
    meta = []
    for f, x in zip(fam, xs):
        s, info = render(f, x, H / 2, scales[f])
        body.append(s); meta.append((f,) + info)
    svg = ('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
           % (W, H, W, H)
           + '<defs><radialGradient id="bg" cx="50%" cy="50%" r="70%">'
             '<stop offset="0" stop-color="#0d1420"/><stop offset="1" stop-color="#04060a"/>'
             '</radialGradient></defs>'
           + '<rect width="%d" height="%d" fill="url(#bg)"/>' % (W, H)
           + "".join(body) + "</svg>")
    open(os.path.join(HERE, "family.svg"), "w").write(svg)
    for name, L, M, n, c in meta:
        print(f"  {name:<7} torus {L:>2} x {M:<3} qubits {n:>3}  checks {c}")
    print("bytes:", len(svg))
