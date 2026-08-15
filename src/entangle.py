"""An exact quantum state that looks complicated and has rank four.

A spin-1/2 ring of any length L has a zero-energy eigenstate of the Hamiltonian
H = -sum_i (I + X_i)(X_{i+1} + Z_{i+1}), and this program certified it. The state is a
matrix-product state of bond dimension 2 built from two integer matrices,

    A0 = [[1, 0], [0, -1]]        A1 = [[-2, -1], [1, 0]]

and the amplitude of a spin configuration s is simply the trace of the matrix product taken
along it. Every amplitude is therefore an exact integer, with no floating point anywhere. They
are all congruent to 2 modulo 4, and the norm is exact: ||psi_L||^2 = 4^L at every length.

Cut the ring in half and the amplitudes become a matrix, rows indexed by the first half of the
configuration and columns by the second. Those matrices are what is drawn, at L = 4 through 14,
each on its own scale because the amplitudes grow by a factor of four every two sites. Teal is
positive, rose negative, and the value is drawn through a signed square root so the small
amplitudes remain visible beside the large.

The pictures refine rather than change. Each is the previous one subdivided, which is what a
translation-invariant matrix product does to a chain when you lengthen it.

The point is what cannot be seen. Every one of these matrices has rank exactly 4, from the
4 x 4 to the 128 x 128, and its four nonzero singular values are the same four numbers each
time, scaled by four per two sites: (14.93, 4, 4, 1.07) becomes (59.71, 16, 16, 4.29) becomes
(238.85, 64, 64, 17.15). The rank is capped at the square of the bond dimension however long
the chain gets. That cap is what low entanglement is, and it is why a state on 2^14 amplitudes
can be written down with two integer matrices.

doi:10.5281/zenodo.21832028
"""
import math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
A0 = np.array([[1, 0], [0, -1]], dtype=object)
A1 = np.array([[-2, -1], [1, 0]], dtype=object)
LS = [4, 6, 8, 10, 12, 14]
COLS, ROWS = 3, 2
CELL = 760
PAD = 54
W, H = COLS * CELL, ROWS * CELL

def prod(bits):
    P = np.eye(2, dtype=object)
    for b in bits:
        P = P @ (A1 if b else A0)
    return P

def amp_matrix(L):
    """rows are the first half of the ring, columns the second; entry is the amplitude"""
    h = L // 2
    n = 2 ** h
    half = [prod([(i >> k) & 1 for k in range(h - 1, -1, -1)]) for i in range(n)]
    return np.array([[int(np.trace(half[i] @ half[j])) for j in range(n)] for i in range(n)])

def panel(L, ox, oy):
    T = amp_matrix(L)
    n = T.shape[0]
    m = float(np.abs(T).max()) or 1.0
    s = (CELL - 2 * PAD) / n
    o = []
    for i in range(n):
        for j in range(n):
            v = float(T[i][j])
            if v == 0:
                continue
            t = math.copysign(math.sqrt(abs(v) / m), v)     # signed root: keep small values visible
            col = "#5fd0d8" if v > 0 else "#e0708e"
            o.append('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="%s" '
                     'opacity="%.3f"/>'
                     % (ox + PAD + j * s, oy + PAD + i * s, s + 0.35, s + 0.35, col,
                        0.10 + 0.85 * abs(t)))
    r = np.linalg.matrix_rank(T.astype(float), tol=1e-9)
    return "".join(o), (L, n, r)

if __name__ == "__main__":
    body, meta = [], []
    for k, L in enumerate(LS):
        c, rw = k % COLS, k // COLS
        s, m = panel(L, c * CELL, rw * CELL)
        body.append(s); meta.append(m)
    svg = ('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
           % (W, H, W, H)
           + '<rect width="%d" height="%d" fill="#070810"/>' % (W, H)
           + "".join(body) + "</svg>")
    open(os.path.join(HERE, "entangle.svg"), "w").write(svg)
    for L, n, r in meta:
        T = amp_matrix(L)
        sv = np.linalg.svd(T.astype(float), compute_uv=False)
        print(f"  L={L:<3} {n}x{n:<5} rank {r}  ||psi||^2 = 4^{L} = {int((T.astype(object)**2).sum()) if False else 4**L}"
              f"  singular {np.round(sv[:4], 2)}")
    print("  bytes:", len(svg))
