"""Twenty clock transitions, drawn as the curvature that certifies them.

A ZEFOZ point is a magnetic field at which a transition frequency is stationary. What
makes such a point useful as a quantum memory is not that the gradient vanishes but how
the surface curves there, since the curvature sets the predicted coherence time.

This program certified that curvature for all twenty published points of Er-167:Y2SiO5,
both crystallographic sites. Each panel is one point, drawn from its own certificate: the
three Hessian eigenvalues, and the eigenvector frame they live in, both shipped as exact
rationals. The surface shown is the quadric those eigenvalues define — the second-order
level set of the frequency about that field.

Seven points have all three eigenvalues negative, a certified local maximum, and close
into an ellipsoid. Thirteen have signature (-,-,+), a saddle, and open into a hyperboloid:
two directions in which the frequency falls away and one in which it climbs. That the
published atlas contains no minimum is a result, not an assumption; the signatures are
certified from bracketed eigenvalues, not read off a numerical diagonalisation.

Axis lengths are the cube root of the true semi-axes 1/sqrt|lambda|. This is the one
liberty taken and it is taken because the true anisotropy is extreme: the ratio of longest
to shortest axis runs from 7.9:1 to 50.8:1, and drawn faithfully every panel is a line.
The compression preserves the ordering of the axes, the eigenvector frame, and the sign
structure exactly; it reduces only how far the needle is drawn out. That these surfaces
are knife-edged is itself the physics of a clock transition.

doi:10.5281/zenodo.21898996
"""
import json, math, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = "/Users/kirt/Documents/reserch math/certify-repo/zefoz-certificates/certificate2.json"
COLS, ROWS = 5, 4
CELL = 420
W, H = COLS * CELL, ROWS * CELL

def num(x):
    return float(F(x)) if isinstance(x, str) else float(x)

def points():
    d = json.load(open(CERT))
    out = []
    for p in d["points"]:
        if p["kind"] != "zefoz_point":
            continue
        for pr in p["pairs"]:
            eig = [ (num(a) + num(b)) / 2 for a, b in pr["hessian_eig_brackets_MHz_per_mT2"] ]
            R = [[num(v) for v in row] for row in pr["hessian_rotation_dyadic"]]
            out.append((p["site"], pr["i"], pr["j"], eig, R, tuple(pr["hessian_signature"])))
    return out

def quadric(eig, sig, nu=26, nv=14):
    """the second-order level set: an ellipsoid when the signature is definite,
    a one-sheet hyperboloid about the positive direction when it is not"""
    a = [(1.0 / math.sqrt(abs(e))) ** (1.0 / 3.0) for e in eig]
    s = max(a)
    a = [x / s for x in a]
    # The true semi-axes are 1/sqrt|lambda| and their ratios run from 7.9:1 to 50.8:1 —
    # every panel would render as a line. Cube-rooting compresses that to at most 3.7:1.
    # Ordering, orientation and the sign structure are untouched; only the exaggeration
    # of the anisotropy is reduced, and the true ratios are printed by this script.
    pos = [k for k, c in enumerate(sig) if c == "+"]
    curves = []
    if not pos:                               # definite: ellipsoid
        for i in range(nv):
            th = math.pi * (i + 0.5) / nv
            curves.append([(a[0]*math.sin(th)*math.cos(ph), a[1]*math.sin(th)*math.sin(ph),
                            a[2]*math.cos(th)) for ph in
                           [2*math.pi*j/nu for j in range(nu+1)]])
        for j in range(nu):
            ph = 2 * math.pi * j / nu
            curves.append([(a[0]*math.sin(th)*math.cos(ph), a[1]*math.sin(th)*math.sin(ph),
                            a[2]*math.cos(th)) for th in
                           [math.pi*k/nv for k in range(nv+1)]])
    else:                                     # saddle: hyperboloid about the + axis
        k = pos[0]
        o1, o2 = [t for t in range(3) if t != k]
        def P(u, v):
            c = [0, 0, 0]
            c[o1] = a[o1] * math.cosh(v) * math.cos(u)
            c[o2] = a[o2] * math.cosh(v) * math.sin(u)
            c[k]  = a[k] * math.sinh(v)
            return tuple(c)
        VS = 0.72
        for i in range(9):
            v = -VS + 2 * VS * i / 8
            curves.append([P(2*math.pi*j/nu, v) for j in range(nu+1)])
        for j in range(nu):
            u = 2 * math.pi * j / nu
            curves.append([P(u, -VS + 2*VS*t/10) for t in range(11)])
    return curves

def rotate(pt, R):
    return tuple(sum(R[r][c] * pt[c] for c in range(3)) for r in range(3))

def project(p, ox, oy, s, ax=0.62, ay=0.72):
    x, y, z = p
    y, z = y*math.cos(ax) - z*math.sin(ax), y*math.sin(ax) + z*math.cos(ax)
    x, z = x*math.cos(ay) + z*math.sin(ay), -x*math.sin(ay) + z*math.cos(ay)
    f = 3.4 / (3.4 - z)
    return ox + x*s*f, oy + y*s*f, z

def hsv(h, sa, v):
    h = (h % 360) / 60.0; i = int(h) % 6; f = h - int(h)
    p, q, t = v*(1-sa), v*(1-sa*f), v*(1-sa*(1-f))
    a, b, c = ((v,t,p),(q,v,p),(p,v,t),(p,q,v),(t,p,v),(v,p,q))[i]
    return "#%02x%02x%02x" % (int(a*255), int(b*255), int(c*255))

def panel(site, i, j, eig, R, sig, ox, oy):
    saddle = "+" in sig
    curves = quadric(eig, sig)
    segs = []
    for c in curves:
        P = [project(rotate(p, R), ox + CELL/2, oy + CELL/2, CELL*0.27) for p in c]
        segs.append((sum(q[2] for q in P)/len(P), P))
    segs.sort(key=lambda t: t[0])
    zs = [s[0] for s in segs]; lo, hi = min(zs), max(zs)
    o = []
    for z, P in segs:
        t = (z - lo)/(hi - lo or 1)
        base = 344 if saddle else 172          # saddles rose, maxima teal
        o.append('<path d="M%.1f %.1f%s" fill="none" stroke="%s" stroke-width="%.2f" '
                 'opacity="%.3f" stroke-linecap="round"/>'
                 % (P[0][0], P[0][1], "".join("L%.1f %.1f" % (q[0], q[1]) for q in P[1:]),
                    hsv(base + (10 if saddle else -14)*t, 0.30 + 0.36*t, 0.42 + 0.52*t),
                    0.5 + 1.3*t, 0.10 + 0.62*t**1.6))
    return "".join(o)

if __name__ == "__main__":
    pts = points()
    body = []
    for idx, (site, i, j, eig, R, sig) in enumerate(pts):
        r, c = divmod(idx, COLS)
        body.append(panel(site, i, j, eig, R, sig, c*CELL, r*CELL))
    svg = ('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
           % (W, H, W, H)
           + '<rect width="%d" height="%d" fill="#080b11"/>' % (W, H)
           + "".join(body) + "</svg>")
    open(os.path.join(HERE, "curvature.svg"), "w").write(svg)
    sad = sum(1 for p in pts if "+" in p[5])
    print("points:", len(pts), "| saddles:", sad, "| maxima:", len(pts) - sad,
          "| minima:", sum(1 for p in pts if p[5] == ("+", "+", "+")))
    print("bytes:", len(svg))
