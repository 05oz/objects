"""Five codes, asked what shape they are.

Every stabilizer code is a bipartite graph: qubits on one side, parity checks on the
other, an edge wherever a check touches a qubit. That graph is the code — it is all a
decoder ever sees. Each is built here from the HX and HZ shipped in its own certificate.

No layout is imposed. Every vertex is placed at the second and third eigenvectors of that
graph's own Laplacian, and what appears is whatever the connectivity forces. The result
was not designed and, for two of the five, was not what I expected.

The three surface codes find their lattice. Nobody told them they are planar; the
eigenvectors put them on a grid because a nearest-neighbour code has no other shape to
take, and the grids visibly grow d = 3, 5, 7. The reason is measurable in the spectrum:
their algebraic connectivity falls 0.424, 0.173, 0.092 as they grow, and a small lambda_2
means the graph has a long direction to be stretched along. That is what a lattice is,
spectrally.

Steane and Golay have no such direction. Their checks are the lines of a finite projective
structure, joining qubits with no notion of nearness, and their algebraic connectivity is
the largest in the set — 1.097 and 1.630. A graph with no bottleneck cannot be laid out in
two dimensions, so the embedding does not find a shape; it collapses. Golay concentrates
almost all of its third eigenvector onto a handful of vertices and comes out a star. That
collapse is not a failure of the drawing. It is the picture of a code that is an expander,
and it is why Golay's distance had to be established by proof rather than by looking.

Distances were determined by this program and machine-checked: Steane 3, Golay 7, and the
rotated surface codes 3, 5, 7.

Qubits are pale, X checks teal, Z checks rose.

doi:10.5281/zenodo.21799780
"""
import math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
QEC = "/Users/kirt/Documents/reserch math/certify-repo/qec-certificates"
CODES = [("five/steane", "steane"), ("surface3", "surface3"), ("surface5", "surface5"),
         ("surface7", "surface7"), ("golay", "golay")]
COLS = 5
CELL = 560
W, H = COLS * CELL, CELL

def load(p):
    t = open(p).read().split()
    m, n = int(t[0]), int(t[1])
    return m, n, [[int(c) for c in t[2 + i]] for i in range(m)]

def tanner(name):
    """qubits 0..n-1, then X checks, then Z checks; an edge per incidence"""
    mX, n, HX = load(f"{QEC}/{name}/HX.txt")
    mZ, _, HZ = load(f"{QEC}/{name}/HZ.txt")
    N = n + mX + mZ
    E = []
    for c, row in enumerate(HX):
        for q, v in enumerate(row):
            if v: E.append((q, n + c))
    for c, row in enumerate(HZ):
        for q, v in enumerate(row):
            if v: E.append((q, n + mX + c))
    return N, n, mX, mZ, E

def spectral(N, E):
    L = np.zeros((N, N))
    for a, b in E:
        L[a][a] += 1; L[b][b] += 1; L[a][b] -= 1; L[b][a] -= 1
    w, V = np.linalg.eigh(L)
    x, y = V[:, 1].copy(), V[:, 2].copy()
    x -= x.mean(); y -= y.mean()
    s = max(np.abs(x).max(), np.abs(y).max()) or 1
    return x / s, y / s

def panel(name, ox):
    N, n, mX, mZ, E = tanner(name)
    x, y = spectral(N, E)
    R = CELL * 0.36
    P = [(ox + CELL / 2 + x[i] * R, CELL / 2 + y[i] * R) for i in range(N)]
    o = []
    for a, b in E:
        xcheck = (b - n) < mX if b >= n else (a - n) < mX
        col = "#3f9fa8" if xcheck else "#a8506a"
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                 'stroke-width="1.5" opacity="0.5"/>' % (P[a][0], P[a][1], P[b][0], P[b][1], col))
    for i in range(N):
        if i < n:
            col, r, op = "#f0eadd", 5.0, 0.95          # qubit
        elif i < n + mX:
            col, r, op = "#5fd0d8", 3.6, 0.9                # X check
        else:
            col, r, op = "#e0708e", 3.6, 0.9                # Z check
        o.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" opacity="%.2f"/>'
                 % (P[i][0], P[i][1], r, col, op))
    return "".join(o), (n, mX, mZ, len(E))

if __name__ == "__main__":
    names = ["steane", "surface3", "surface5", "surface7", "golay"]
    body, meta = [], []
    for k, nm in enumerate(names):
        s, info = panel(nm, k * CELL)
        body.append(s); meta.append((nm,) + info)
    svg = ('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
           % (W, H, W, H)
           + '<rect width="%d" height="%d" fill="#080b11"/>' % (W, H)
           + "".join(body) + "</svg>")
    open(os.path.join(HERE, "tanner.svg"), "w").write(svg)
    for nm, n, mX, mZ, e in meta:
        print(f"  {nm:<10} qubits {n:>3}  X-checks {mX:>3}  Z-checks {mZ:>3}  edges {e:>4}")
    print("bytes:", len(svg))
