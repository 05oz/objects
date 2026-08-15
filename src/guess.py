"""The decoder never sees the error. It sees which detectors fired, and it guesses.

This is the whole difficulty of quantum error correction, and it can be replayed exactly from
the certificate for the distance-3 rotated surface code under one round of circuit-level noise.

Twenty-three fault mechanisms. Each fires a known set of detectors, and some flip the logical
observable. When faults occur you do not learn which ones: you learn only the total set of
detectors that fired, because measuring the faults themselves would destroy the state you are
protecting. From that alone the decoder must decide what to undo.

The rule is stated in the certificate and reproduced here from scratch: breadth-first search
over the 256 possible syndromes starting from the empty one, mechanisms tried in index order,
first assignment wins. That gives, for every syndrome, the simplest fault set that would explain
it. Reconstructing it reproduces the shipped 256-entry decoder table with no disagreement.

Left is what happened. Centre is all the decoder is given. Right is what it concludes. When the
left and right disagree about the logical observable, the qubit has flipped and nothing noticed.

Three pairs of detectors are interchangeable under an automorphism of this graph and land on
the same point under every embedding of it. Because this piece is about the decoding and not
the geometry, they are separated by a small offset in index order; that is the one thing here
not read from the certificate, and the companion piece on the same graph leaves them merged.

The cases shown are read from the certificate, not chosen for effect: first the single faults,
which are always identified exactly, then pairs the decoder handles, then pairs it does not.
There are 55 of that last kind and they are enumerated in full, along with 690 triples and
4,078 quadruples. They are why the logical error probability is what it is, and why the number
is exact rather than sampled.

    python3 guess.py frames      writes the animation frames

doi:10.5281/zenodo.21895825
"""
import json, math, os, sys
from collections import deque
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = ("/Users/kirt/Documents/reserch math/certify-repo/wedge-certificates/"
        "certificate_d3_r1_p1over100.json")
PW, H = 720, 720
W = PW * 3
HOLD = 4                                   # frames per reveal stage

D = json.load(open(CERT))
ND = D["num_detectors"]
MECH = [(i, [k for k in range(ND) if (m["det"] >> k) & 1], m["obs"], m["det"])
        for i, m in enumerate(D["mechanisms"])]

def decoder():
    """the certificate's own rule, rebuilt: BFS from the empty syndrome, index order, first wins"""
    pred = {0: (0, 0, ())}
    q = deque([0])
    while q:
        s = q.popleft(); w, o, st = pred[s]
        for i, _, ob, det in MECH:
            t = s ^ det
            if t not in pred:
                pred[t] = (w + 1, o ^ ob, st + (i,))
                q.append(t)
    return pred

PRED = decoder()

def layout():
    A = np.zeros((ND, ND))
    for _, b, _, _ in MECH:
        if len(b) == 2:
            A[b[0]][b[1]] += 1; A[b[1]][b[0]] += 1
    L = np.diag(A.sum(1)) - A
    _, V = np.linalg.eigh(L)
    x, y = V[:, 1].copy(), V[:, 2].copy()
    x -= (x.max() + x.min()) / 2; y -= (y.max() + y.min()) / 2
    s = max(np.abs(x).max(), np.abs(y).max()) * 1.34 or 1
    R = PW * 0.40
    pts = [[x[i] / s * R, y[i] / s * R] for i in range(ND)]
    # Three detector pairs are interchangeable under an automorphism of this graph and land on
    # the same point under every embedding. This piece is about the decoding, not the geometry,
    # so coincident detectors are separated by a small offset in index order and the fact is
    # recorded here rather than hidden.
    for i in range(ND):
        for j in range(i + 1, ND):
            if math.hypot(pts[i][0] - pts[j][0], pts[i][1] - pts[j][1]) < 1e-6:
                for k, sgn in ((i, -1), (j, +1)):
                    n = math.hypot(pts[k][0], pts[k][1]) or 1
                    pts[k][0] += -pts[k][1] / n * PW * 0.058 * sgn
                    pts[k][1] += pts[k][0] / n * PW * 0.058 * sgn * 0
                    pts[k][1] += sgn * PW * 0.030
    return [(a, b) for a, b in pts]

P0 = layout()

def graph(ox, cy, active, fired, col, faint):
    """one copy of the matching graph; `active` mechanisms lit, `fired` detectors lit"""
    cx = ox + PW / 2
    pt = [(cx + p[0], cy + p[1]) for p in P0]
    def bnd(p):
        dx, dy = p[0] - cx, p[1] - cy
        n = math.hypot(dx, dy) or 1
        return (p[0] + dx / n * PW * 0.10, p[1] + dy / n * PW * 0.10)
    o = []
    for i, b, ob, _ in MECH:
        a, c = (pt[b[0]], pt[b[1]]) if len(b) == 2 else (pt[b[0]], bnd(pt[b[0]]))
        on = i in active
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="%.1f" '
                 'opacity="%.2f" stroke-linecap="round"/>'
                 % (a[0], a[1], c[0], c[1], col if on else "#2c3442",
                    5.4 if on else 1.6, 0.98 if on else (0.30 if not faint else 0.16)))
    for i in range(ND):
        lit = i in fired
        o.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" opacity="%.2f"/>'
                 % (pt[i][0], pt[i][1], 12.5 if lit else 7.0,
                    "#ffe9a8" if lit else "#4a5566", 1.0 if lit else 0.75))
        if lit:
            o.append('<circle cx="%.1f" cy="%.1f" r="22" fill="#ffe9a8" opacity="0.13"/>'
                     % (pt[i][0], pt[i][1]))
    return "".join(o)

def label(ox, y, s, size=27, col="#8b93a1", weight=""):
    return ('<text x="%.1f" y="%.1f" font-family="Georgia,serif" font-size="%d" fill="%s" '
            'text-anchor="middle"%s>%s</text>' % (ox + PW / 2, y, size, col, weight, s))

def frame(case, stage):
    F, syn, obs, leader, pobs = case
    fired = [k for k in range(ND) if (syn >> k) & 1]
    ok = (pobs == obs)
    o = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
         % (W, H, W, H),
         '<rect width="%d" height="%d" fill="#070810"/>' % (W, H)]
    cy = H * 0.50
    o.append(graph(0, cy, F if stage >= 1 else (), fired if stage >= 2 else (), "#ff6f5e", False))
    o.append(graph(PW, cy, (), fired if stage >= 2 else (), "#ffe9a8", True))
    o.append(graph(2 * PW, cy, leader if stage >= 3 else (), fired if stage >= 2 else (),
                   "#5fd0d8", False))
    o.append(label(0, 74, "what happened"))
    o.append(label(PW, 74, "what the decoder is given"))
    o.append(label(2 * PW, 74, "what it concludes"))
    if stage >= 1:
        o.append(label(0, H - 96, "%d fault%s &#183; observable %s"
                       % (len(F), "" if len(F) == 1 else "s", "flips" if obs else "unchanged"),
                       24, "#ff9d92"))
    if stage >= 2:
        o.append(label(PW, H - 96, "%d detector%s fired" % (len(fired), "" if len(fired) == 1 else "s"),
                       24, "#d8c98a"))
    if stage >= 3:
        o.append(label(2 * PW, H - 96, "%d fault%s &#183; observable %s"
                       % (len(leader), "" if len(leader) == 1 else "s",
                          "flips" if pobs else "unchanged"), 24, "#8fdde6"))
    if stage >= 4:
        o.append('<text x="%d" y="%d" font-family="Georgia,serif" font-size="34" fill="%s" '
                 'text-anchor="middle">%s</text>'
                 % (W / 2, H - 34, "#8de0b0" if ok else "#ff6f5e",
                    "corrected" if ok else "the qubit flipped and nothing noticed"))
    o.append("</svg>")
    return "".join(o)

def cases():
    """read from the certificate in order: singles, then pairs it handles, then pairs it does not"""
    us = {tuple(s) for s in D["uncorrectable_sets"]["2"]}
    def mk(F):
        syn = 0; obs = 0
        for i in F:
            syn ^= MECH[i][3]; obs ^= MECH[i][2]
        w, po, leader = PRED[syn]
        return (tuple(F), syn, obs, tuple(leader), po)
    out = [mk((0,)), mk((11,))]
    good = [mk((a, b)) for a in range(23) for b in range(a + 1, 23)
            if (a, b) not in us]
    out += good[:2] + [mk(s) for s in sorted(us)[:3]]
    return out

if __name__ == "__main__":
    C = cases()
    if len(sys.argv) > 1 and sys.argv[1] == "frames":
        os.makedirs(os.path.join(HERE, "gframes"), exist_ok=True)
        k = 0
        for c in C:
            for stage in (0, 1, 2, 3, 4):
                for _ in range(HOLD if stage < 4 else HOLD * 3):
                    open(os.path.join(HERE, "gframes", "f%03d.svg" % k), "w").write(frame(c, stage))
                    k += 1
        print("wrote %d frames across %d cases" % (k, len(C))); raise SystemExit
    open(os.path.join(HERE, "guess.svg"), "w").write(frame(C[-1], 4))
    bad = [s for s in PRED if [PRED[s][0], PRED[s][1]] != D["decoder_table"][str(s)]]
    print("  syndromes reconstructed:", len(PRED), "| disagreements with the shipped table:", len(bad))
    print("  cases:", len(C), "|", [(c[0], "corrected" if c[4] == c[2] else "logical error") for c in C])
