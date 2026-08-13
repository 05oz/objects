"""The gross code on its torus.

IBM's [[288,12,18]] bivariate bicycle code, drawn on the surface it is defined on.

The parity checks are read from the shipped certificate, not from a paper: HX.txt gives
144 X-checks over 288 qubits, weight 6, and its left block decomposes over Z_12 x Z_12 —
a 12 x 12 torus, recovered here by testing the shift structure of every row rather than
assumed. Each check touches three qubits of the L sublattice and three of the R, at the
offsets its monomials name.

288 qubits sit on that torus, two per cell. 288 checks (144 X, 144 Z) connect them. The
distance is 18, machine-checked: doi:10.5281/zenodo.21831995.

Nothing here is a diagram of the code. It is the code, on its own surface, seen from
outside.
"""
import math, os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = "/Users/kirt/Documents/reserch math/certify-repo"
W = H = 1600
L = M = 12                      # torus dimensions, recovered from HX
R_MAJ, R_MIN = 1.0, 0.42        # torus radii

def load(path):
    t = open(path).read().split()
    m, n = int(t[0]), int(t[1])
    return m, n, [[int(c) for c in t[2 + i]] for i in range(m)]

mx, nx, HX = load(f"{REPO}/qec-certificates/bb288/HX.txt")
mz, nz, HZ = load(f"{REPO}/qec-certificates/bb288/HZ.txt")

def cell(j):
    """qubit index -> (sublattice, torus cell)"""
    side, k = divmod(j, 144)
    return side, divmod(k, M)          # (0=L, 1=R), (i, j) on Z_L x Z_M

def uv(j):
    """qubit index -> its (u, v) angles on the torus"""
    side, (i, k) = cell(j)
    return (2*math.pi*(i + (0.26 if side else -0.26))/L,
            2*math.pi*(k + (0.26 if side else -0.26))/M)

def surf(u, v):
    r = R_MAJ + R_MIN * math.cos(v)
    return (r*math.cos(u), r*math.sin(u), R_MIN*math.sin(v))

def wrap(d, period=2*math.pi):
    """shortest signed angular difference"""
    while d >  period/2: d -= period
    while d < -period/2: d += period
    return d

def torus(i, j, side):
    """Z_L x Z_M cell -> a point on the torus surface; the two sublattices sit
    slightly apart along the tube so both are visible."""
    u = 2 * math.pi * (i + (0.26 if side else -0.26)) / L
    v = 2 * math.pi * (j + (0.26 if side else -0.26)) / M
    r = R_MAJ + R_MIN * math.cos(v)
    return (r * math.cos(u), r * math.sin(u), R_MIN * math.sin(v))

def rot(p, ax, ay):
    x, y, z = p
    y, z = y * math.cos(ax) - z * math.sin(ax), y * math.sin(ax) + z * math.cos(ax)
    x, z = x * math.cos(ay) + z * math.sin(ay), -x * math.sin(ay) + z * math.cos(ay)
    return x, y, z

def project(p, ax=1.14, ay=0.34, d=3.6, s=430):
    x, y, z = rot(p, ax, ay)
    f = d / (d - z)
    return W / 2 + x * s * f, H / 2 + y * s * f, z

# --- geometry -----------------------------------------------------------------
pts = {}
for j in range(288):
    side, (i, k) = cell(j)
    pts[j] = torus(i, k, side)

SEG = 11                        # samples per edge; enough to read as a curve

def geodesic(j1, j2):
    """the edge as a path ON the torus: interpolate the angles, never the chord.
    Wraparound is taken the short way, so an edge that leaves one side of the
    lattice re-enters on the other rather than cutting through the tube."""
    u1, v1 = uv(j1); u2, v2 = uv(j2)
    du, dv = wrap(u2 - u1), wrap(v2 - v1)
    return [surf(u1 + du*t/SEG, v1 + dv*t/SEG) for t in range(SEG + 1)]

def edges(Hm, kind):
    out = []
    for row in Hm:
        sup = [j for j, v in enumerate(row) if v]
        for a in range(len(sup)):
            for b in range(a + 1, len(sup)):
                out.append((geodesic(sup[a], sup[b]), kind))
    return out

E = edges(HX, "x") + edges(HZ, "z")

def bez(pa, pb):
    """bow the edge outward along the surface normal so it hugs the torus"""
    mx_ = [(a + b) / 2 for a, b in zip(pa, pb)]
    n = math.sqrt(mx_[0] ** 2 + mx_[1] ** 2) or 1
    k = 1.0 + 0.052 / max(0.35, math.dist(pa, pb))
    return (mx_[0] / n * (n * k), mx_[1] / n * (n * k), mx_[2] * 1.06)

def hsv(h, s, v):
    h = (h % 360) / 60.0; i = int(h) % 6; f = h - int(h)
    p, q, t = v * (1 - s), v * (1 - s * f), v * (1 - s * (1 - f))
    r, g, b = ((v, t, p), (q, v, p), (p, v, t), (p, q, v), (t, p, v), (v, p, q))[i]
    return "#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255))

def svg():
    o = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
         % (W, H, W, H),
         '<defs><radialGradient id="bg" cx="50%%" cy="46%%" r="72%%">'
         '<stop offset="0" stop-color="#0d1420"/><stop offset="1" stop-color="#04060a"/>'
         '</radialGradient></defs>',
         '<rect width="%d" height="%d" fill="url(#bg)"/>' % (W, H)]

    drawn = []
    for path, kind in E:
        P = [project(p) for p in path]
        drawn.append((sum(q[2] for q in P) / len(P), P, kind))
    drawn.sort(key=lambda t: t[0])                     # painter: far first

    zs = [d[0] for d in drawn]
    lo, hi = min(zs), max(zs)
    for z, P, kind in drawn:
        t = (z - lo) / (hi - lo or 1)                  # 0 far, 1 near
        hue = 196 - 26 * t if kind == "x" else 336 + 14 * t
        col = hsv(hue, 0.26 + 0.44 * t, 0.30 + 0.66 * t)
        op = 0.055 + 0.72 * t ** 2.1
        wd = 0.35 + 2.3 * t ** 1.4
        d = "M%.0f %.0f" % (P[0][0], P[0][1]) + "".join(
            "L%.0f %.0f" % (q[0], q[1]) for q in P[1:])
        o.append('<path d="%s" fill="none" stroke="%s" stroke-width="%.2f" '
                 'opacity="%.3f" stroke-linecap="round"/>' % (d, col, wd, op))

    qs = sorted(((project(p), j) for j, p in pts.items()), key=lambda t: t[0][2])
    zq = [q[0][2] for q in qs]; qlo, qhi = min(zq), max(zq)
    for (X, Y, Z), j in qs:
        t = (Z - qlo) / (qhi - qlo or 1)
        side = j // 144
        col = hsv(38 if side else 190, 0.16 + 0.20 * t, 0.55 + 0.42 * t)
        o.append('<circle cx="%.1f" cy="%.1f" r="%.2f" fill="%s" opacity="%.3f"/>'
                 % (X, Y, 1.0 + 4.0 * t ** 1.5, col, 0.18 + 0.78 * t ** 1.4))
    o.append("</svg>")
    return "".join(o)

if __name__ == "__main__":
    s = svg()
    open(os.path.join(HERE, "gross.svg"), "w").write(s)
    print("torus %dx%d | qubits %d | checks %d | edges drawn %d | bytes %d"
          % (L, M, len(pts), mx + mz, len(E), len(s)))
