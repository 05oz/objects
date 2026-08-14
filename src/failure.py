"""What a decoder sees, and where the code breaks.

A decoder never sees qubits. It sees detectors that fire, and must decide which faults did it.
That view is a graph, and it is read here out of the certificates rather than drawn from any
picture of the code: one edge per fault mechanism, its endpoints taken from the detector bitmask
the mechanism carries. At distance 3 that is 8 detectors and 23 mechanisms, 15 joining a pair of
detectors and 8 firing a single detector and so running to the boundary. At distance 5 it is 24
detectors and 77 mechanisms, 65 and 12.

Both are placed at the second and third eigenvectors of their own Laplacian. Nothing else.

The distance-3 graph will not lie flat. Three pairs of its detectors land at exactly the same
point, and no spectral coordinate separates them, because the swap 1↔4, 2↔5, 3↔6 is a genuine
automorphism of the graph: those detectors are interchangeable, and a symmetric object has no
asymmetric drawing. They are marked with a double ring rather than nudged apart. At distance 5
the symmetry is gone and all 24 detectors separate cleanly.

On the distance-3 panel every edge is weighted by how much of the code's failure it carries. The
certificate enumerates all 4,823 sets of faults that defeat the decoder, exactly — 55 of weight
two, 690 of weight three, 4,078 of weight four — and each mechanism is drawn thicker and brighter
the more of them it belongs to, from 686 to 1,021. The faint arcs are the 55 weight-two sets,
the minimum-weight failures that hold the distance at 3 and that carry nearly all of the logical
error probability at low physical error rate.

The distance-5 panel carries no such weighting. Its certificate ships per-weight SHA-256 digests
of the uncorrectable sets rather than the sets themselves, so the enumeration cannot be redrawn
from it and no edge weight is claimed.

Four mechanisms flip the logical observable on their own, in red. Every one is a boundary fault,
at detectors 1, 3, 5 and 7. That is the logical operator meeting the edge of the patch, and
nothing here arranged it.

doi:10.5281/zenodo.21895825
"""
import json, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CDIR = "/Users/kirt/Documents/reserch math/certify-repo/wedge-certificates"
PANELS = [("certificate_d3_r1_p1over100.json", "d = 3"),
          ("certificate_d5_r1_p1over100.json", "d = 5")]
PW, H = 1240, 1300
W = PW * len(PANELS)

def load(fn):
    d = json.load(open(os.path.join(CDIR, fn)))
    nd = d["num_detectors"]
    edges = [(i, [k for k in range(nd) if (m["det"] >> k) & 1], m["obs"])
             for i, m in enumerate(d["mechanisms"])]
    sets = {int(w): [tuple(s) for s in v]
            for w, v in d.get("uncorrectable_sets", {}).items() if v}
    return d, nd, edges, sets

def layout(nd, edges, cx, cy, R, stub):
    """eigenvectors 2 and 3, then framed on the bounding box of the whole figure including the
    boundary stubs. Uniform scale in both axes: the found geometry is never stretched to fit."""
    A = np.zeros((nd, nd))
    for _, b, _ in edges:
        if len(b) == 2:
            A[b[0]][b[1]] += 1; A[b[1]][b[0]] += 1
    L = np.diag(A.sum(1)) - A
    _, V = np.linalg.eigh(L)
    x, y = V[:, 1].copy(), V[:, 2].copy()
    mx, my = (x.max() + x.min()) / 2, (y.max() + y.min()) / 2
    x -= mx; y -= my
    # boundary stubs extend outward, so include them when measuring the extent
    ext = 1.0 + stub
    s = max(np.abs(x).max(), np.abs(y).max()) * ext or 1
    return [(cx + x[i] / s * R, cy + y[i] / s * R) for i in range(nd)]

def coincident(P):
    """detectors the spectrum cannot separate — an automorphism orbit, not a drawing fault"""
    out = set()
    for i in range(len(P)):
        for j in range(i + 1, len(P)):
            if math.hypot(P[i][0] - P[j][0], P[i][1] - P[j][1]) < 0.5:
                out.add(i); out.add(j)
    return out

def panel(fn, label, ox):
    d, nd, edges, sets = load(fn)
    cx, cy, R = ox + PW / 2, H * 0.47, PW * 0.40
    STUB = 0.30
    P = layout(nd, edges, cx, cy, R, STUB)
    dup = coincident(P)
    tot = sum(len(v) for v in sets.values())
    part = ([sum(1 for g in sets.values() for s in g if i in s) for i, _, _ in edges]
            if tot else [1] * len(edges))
    lo, hi = min(part), max(part)

    def bpt(p):
        dx, dy = p[0] - cx, p[1] - cy
        n = math.hypot(dx, dy) or 1
        return (p[0] + dx / n * R * STUB, p[1] + dy / n * R * STUB)
    ends = {i: ((P[b[0]], P[b[1]]) if len(b) == 2 else (P[b[0]], bpt(P[b[0]])))
            for i, b, _ in edges}

    o = []
    if sets.get(2):                     # the minimum-weight failures, as arcs between mechanisms
        mid = lambda i: (((ends[i][0][0] + ends[i][1][0]) / 2), ((ends[i][0][1] + ends[i][1][1]) / 2))
        for a, b in sets[2]:
            p, q = mid(a), mid(b)
            mx, my = (p[0] + q[0]) / 2, (p[1] + q[1]) / 2
            o.append('<path d="M%.1f %.1fQ%.1f %.1f %.1f %.1f" fill="none" stroke="#ffd27f" '
                     'stroke-width="1.1" opacity="0.26" stroke-linecap="round"/>'
                     % (p[0], p[1], mx + (cx - mx) * 0.32, my + (cy - my) * 0.32, q[0], q[1]))

    for i, b, obs in edges:
        t = (part[i] - lo) / (hi - lo or 1)
        (a, c) = ends[i]
        col = "#ff6f5e" if obs else ("#7fd4e0" if len(b) == 2 else "#8f9bad")
        wdt = (1.5 + 6.5 * t) if tot else 2.0
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="%.2f" '
                 'opacity="%.2f" stroke-linecap="round"/>'
                 % (a[0], a[1], c[0], c[1], col, wdt, (0.32 + 0.52 * t) if tot else 0.62))
        if len(b) == 1:
            o.append('<circle cx="%.1f" cy="%.1f" r="3.0" fill="%s" opacity="0.7"/>' % (c[0], c[1], col))

    r = 9.0 if nd <= 8 else 6.0
    for i in range(nd):
        o.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="#0b0d13" opacity="0.95"/>'
                 % (P[i][0], P[i][1], r + 2.5))
        o.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="#f3ecdd" opacity="0.96"/>'
                 % (P[i][0], P[i][1], r))
        if i in dup:                    # two interchangeable detectors at one point
            o.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="none" stroke="#f3ecdd" '
                     'stroke-width="1.3" opacity="0.55"/>' % (P[i][0], P[i][1], r + 6.5))

    o.append('<text x="%.1f" y="%.1f" font-family="Georgia,serif" font-size="34" fill="#7f858f" '
             'text-anchor="middle">%s</text>' % (cx, H - 84, label))
    return "".join(o), (nd, len(edges), tot, len(dup))

if __name__ == "__main__":
    body, meta = [], []
    for k, (fn, label) in enumerate(PANELS):
        s, m = panel(fn, label, k * PW)
        body.append(s); meta.append((label,) + m)
    svg = ('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
           % (W, H, W, H)
           + '<rect width="%d" height="%d" fill="#070810"/>' % (W, H)
           + "".join(body) + "</svg>")
    open(os.path.join(HERE, "failure.svg"), "w").write(svg)
    for label, nd, ne, tot, dup in meta:
        print(f"  {label}: {nd} detectors, {ne} mechanisms, "
              f"{tot if tot else 'no'} uncorrectable sets, {dup} coincident detectors")
    print("  bytes:", len(svg))
