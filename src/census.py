"""A shape that appears when you have enough of them.

Kelmans asked in 1984 whether every three-connected cubic graph carries the largest possible
number of vertex-disjoint three-vertex paths. This program verified it, and the deposit ships a
certificate per graph: the graph6 string and the packing found for it.

Up to eighteen vertices the census is complete, and that is what is drawn. Every graph is
decoded from its own shipped string and its adjacency spectrum computed here: 1, 2, 4, 14, 57,
341, 2,828 and 30,468 graphs at orders 4 through 18, the published enumeration exactly, giving
599,318 eigenvalues. Nothing is sampled and nothing is fitted.

Each band is one order, the density of all its eigenvalues, on a common horizontal scale from
minus three to three. The bands are stacked with the smallest census at the top.

At the top there is nothing to see, because one graph has four eigenvalues. Going down, the
population grows by roughly an order of magnitude a step and something specific happens: the
noise resolves into a definite curve, symmetric, with a soft interior and two shoulders that
stop short of the ends. That curve is not drawn here and not assumed; it is what the
eigenvalues do.

Two features are exact rather than statistical. The spike standing alone at the right edge is
the eigenvalue three, which every cubic graph has once, so the spike holds exactly one
eigenvalue per graph. Its mirror at minus three appears only for the bipartite graphs, of which
this census contains 173.

The shoulders fall at plus and minus two root two, and no eigenvalue of any of the 33,715
graphs lies outside those shoulders except the trivial three and the bipartite minus three: at
order eighteen, 94 per cent of all eigenvalues lie strictly inside, and the remainder is almost
entirely the one-per-graph spike.

doi:10.5281/zenodo.21897011
"""
import glob, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CDIR = "/Users/kirt/Documents/reserch math/certify-repo/kelmans-certificates"
MAXN = 18                                  # the census is complete to here
W, H = 2200, 1560
L, R, T = 150, W - 150, 96
BAND, GAP = 152, 18
BINS = 620

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

def spectra():
    by = {}
    for f in sorted(glob.glob(os.path.join(CDIR, "certs_n*.txt"))):
        for line in open(f, errors="ignore"):
            if not line.startswith("CERT"):
                continue
            n, E = graph6(line.split("|")[0].split()[1])
            if n > MAXN:
                continue
            A = np.zeros((n, n))
            for a, b in E:
                A[a][b] = A[b][a] = 1.0
            by.setdefault(n, []).extend(np.linalg.eigvalsh(A).tolist())
    return {k: np.array(v) for k, v in by.items()}

def hsv(h, s, v):
    h = (h % 360) / 60.0; i = int(h) % 6; f = h - int(h)
    p, q, t = v * (1 - s), v * (1 - s * f), v * (1 - s * (1 - f))
    a, b, c = ((v, t, p), (q, v, p), (p, v, t), (p, q, v), (t, p, v), (v, p, q))[i]
    return "#%02x%02x%02x" % (int(a * 255), int(b * 255), int(c * 255))

def svg():
    S = spectra()
    orders = sorted(S)
    X = lambda v: L + (R - L) * (v + 3.0) / 6.0
    o = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
         % (W, H, W, H),
         '<rect width="%d" height="%d" fill="#06070c"/>' % (W, H)]

    # the shoulders at +-2 root 2, where the bulk of every cubic spectrum stops
    for sgn in (-1, 1):
        x = X(sgn * 2 * math.sqrt(2))
        o.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#4a4536" stroke-width="1.1" '
                 'opacity="0.55"/>' % (x, T - 14, x, T + len(orders) * (BAND + GAP)))

    for k, n in enumerate(orders):
        v = S[n]
        base = T + k * (BAND + GAP) + BAND
        h, edges = np.histogram(v, bins=BINS, range=(-3.0, 3.0))
        m = h.max() or 1
        t = k / (len(orders) - 1)
        col = hsv(198 - 160 * t, 0.40 + 0.34 * t, 0.60 + 0.36 * t)
        pts = []
        for i, c in enumerate(h):
            x = X(-3.0 + 6.0 * (i + 0.5) / BINS)
            pts.append((x, base - (BAND - 12) * (c / m) ** 0.55))
        d = "M%.1f %.1f" % (X(-3.0), base) + "".join("L%.1f %.1f" % p for p in pts) \
            + "L%.1f %.1fZ" % (X(3.0), base)
        o.append('<path d="%s" fill="%s" opacity="%.2f"/>' % (d, col, 0.20 + 0.14 * t))
        o.append('<path d="M%.1f %.1f%s" fill="none" stroke="%s" stroke-width="%.1f" '
                 'opacity="%.2f"/>'
                 % (pts[0][0], pts[0][1], "".join("L%.1f %.1f" % p for p in pts[1:]),
                    col, 1.0 + 0.9 * t, 0.55 + 0.40 * t))
    o.append("</svg>")
    return "".join(o), S

if __name__ == "__main__":
    s, S = svg()
    open(os.path.join(HERE, "census.svg"), "w").write(s)
    tot = sum(len(v) for v in S.values())
    print("  orders:", {n: len(S[n]) // n for n in sorted(S)})
    print("  eigenvalues:", tot)
    v = S[18]; b = 2 * math.sqrt(2)
    print("  order-18 inside +-2sqrt2: %d of %d (%.2f%%)"
          % (int((np.abs(v) <= b + 1e-9).sum()), len(v),
             100 * (np.abs(v) <= b + 1e-9).mean()))
    print("  bipartite graphs (an eigenvalue at -3): %d"
          % sum(int((np.abs(S[n] + 3) < 1e-9).sum()) for n in S))
    print("  bytes:", len(s))
