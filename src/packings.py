"""Certified packings: cubic graphs cut into paths of three.

Kelmans asked in 1984 whether every 3-connected cubic graph carries the largest
possible number of vertex-disjoint 3-vertex paths. This program verified it for all
6,339,157 such graphs on at most 22 vertices; 43,580 certificates ship, each a graph6
string and the packing found for it.

Every panel here is one certificate, read directly: the graph is decoded from its
graph6 string, and the coloured pieces are the packing that certificate records —
floor(v/3) paths of three vertices, no two sharing a vertex. Grey edges are the ones
the packing does not use. Where v is not divisible by 3, one or two vertices are left
over and stay pale; that is not a failure but the exact bound.

The layouts are not drawn either. Each graph places its own vertices at the second and
third eigenvectors of its Laplacian — the standard spectral embedding — so the shape of
a panel is a property of that graph and of nothing else.

doi:10.5281/zenodo.21897011
"""
import math, os, glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = "/Users/kirt/Documents/reserch math/certify-repo/kelmans-certificates"
COLS, ROWS = 6, 4
CELL = 400
W, H = COLS * CELL, ROWS * CELL

PAL = ["#e0703a", "#5aa9e6", "#6fcf8e", "#d9b23c", "#b07fd0", "#4fc4c4", "#e0648a"]

def graph6(s):
    b = [ord(c) - 63 for c in s]
    n = b[0]
    bits = []
    for x in b[1:]:
        bits += [(x >> k) & 1 for k in range(5, -1, -1)]
    E, i = [], 0
    for j in range(1, n):
        for k in range(j):
            if i < len(bits) and bits[i]:
                E.append((k, j))
            i += 1
    return n, E

def spectral(n, E):
    """the graph's own layout: eigenvectors 2 and 3 of its Laplacian"""
    L = np.zeros((n, n))
    for a, b in E:
        L[a][a] += 1; L[b][b] += 1; L[a][b] -= 1; L[b][a] -= 1
    w, V = np.linalg.eigh(L)
    x, y = V[:, 1].copy(), V[:, 2].copy()
    x -= x.mean(); y -= y.mean()
    s = max(np.abs(x).max(), np.abs(y).max()) or 1
    return x / s, y / s

def certificates(limit):
    out = []
    for f in sorted(glob.glob(f"{CERT}/certs_n*.txt")):
        for line in open(f):
            if not line.startswith("CERT"):
                continue
            parts = line.split("|")
            g = parts[0].split()[1]
            pieces = [tuple(int(v) for v in t.split("-"))
                      for t in parts[1].split() if "-" in t]
            n, E = graph6(g)
            out.append((n, E, pieces))
            break                       # one certificate per order-file per pass
    return out[:limit]

def gather(target):
    """Four panels from each order the deposit certifies, largest orders included:
    24 is where the packing is richest (8 pieces) and where the deposit reports
    search-side completeness only."""
    want = [14, 16, 18, 20, 22, 24]
    per = target // len(want)
    picked = {n: [] for n in want}
    for f in sorted(glob.glob(f"{CERT}/certs_n*.txt")):
        for i, line in enumerate(open(f)):
            if not line.startswith("CERT"):
                continue
            parts = line.split("|")
            n, E = graph6(parts[0].split()[1])
            if n not in picked or len(picked[n]) >= per:
                continue
            if i % 97:                      # spread the sample through the file
                continue
            pieces = [tuple(int(v) for v in t.split("-"))
                      for t in parts[1].split() if "-" in t]
            picked[n].append((n, E, pieces))
        if all(len(v) >= per for v in picked.values()):
            break
    out = [g for n in want for g in picked[n]]
    return out[:target]

def panel(n, E, pieces, ox, oy):
    x, y = spectral(n, E)
    R = CELL * 0.34
    P = [(ox + CELL / 2 + x[i] * R, oy + CELL / 2 + y[i] * R) for i in range(n)]
    used = set()
    for p in pieces:
        for a, b in zip(p, p[1:]):
            used.add((min(a, b), max(a, b)))
    covered = {v for p in pieces for v in p}

    o = []
    for a, b in E:                                   # edges the packing does not use
        if (min(a, b), max(a, b)) not in used:
            o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#39414f" '
                     'stroke-width="1.2" opacity="0.55"/>'
                     % (P[a][0], P[a][1], P[b][0], P[b][1]))
    for k, p in enumerate(pieces):                   # the packing itself
        col = PAL[k % len(PAL)]
        d = "M%.1f %.1f" % P[p[0]] + "".join("L%.1f %.1f" % P[v] for v in p[1:])
        o.append('<path d="%s" fill="none" stroke="%s" stroke-width="4.2" opacity="0.92" '
                 'stroke-linecap="round" stroke-linejoin="round"/>' % (d, col))
    for i in range(n):
        inside = i in covered
        o.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" opacity="%.2f"/>'
                 % (P[i][0], P[i][1], 4.6 if inside else 3.4,
                    "#f2ecdf" if inside else "#6b7484", 0.95 if inside else 0.85))
    return "".join(o)

if __name__ == "__main__":
    panels = gather(COLS * ROWS)
    body = []
    for idx, (n, E, pieces) in enumerate(panels):
        r, c = divmod(idx, COLS)
        body.append(panel(n, E, pieces, c * CELL, r * CELL))
    svg = ('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
           % (W, H, W, H)
           + '<rect width="%d" height="%d" fill="#0a0d13"/>' % (W, H)
           + "".join(body) + "</svg>")
    open(os.path.join(HERE, "packings.svg"), "w").write(svg)
    for n, E, pieces in panels[:6]:
        print(f"  n={n:<3} edges={len(E):<3} pieces={len(pieces)} "
              f"covered={len({v for p in pieces for v in p})}/{n}  floor(n/3)={n//3}")
    print("panels:", len(panels), "| bytes:", len(svg))
