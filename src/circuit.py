"""Twelve faults no detector sees.

The circuit is the single-auxiliary colour-code syndrome-extraction circuit of Kishony and
Fowler, with their default schedule, at distance fifteen, run for fifteen rounds with two-qubit
depolarising noise after every CNOT. Its code needs fifteen data errors to flip the logical
qubit. The circuit can be beaten with twelve faults, and the certificate ships them: twelve
fault mechanisms, each named by the detectors it flips, which together flip no detector at all
and flip the logical observable.

That twelve is also the least is proved about an extended code built from the circuit, and
carried to the circuit by two hypotheses that a shipped script checks: eleven or fewer faults
cannot do it. That proof is too large for the archive. It is attached to release v0.16.0 of the
certify repository on GitHub, and REGENERATE.md there rebuilds it.

Everything here is read from the shipped files. Qubits and detectors stand at the coordinates
the circuit itself declares, QUBIT_COORDS and DETECTOR with SHIFT_COORDS applied, with the
vertical axis pointing down as those coordinates do, so the observable's row is at the top. The
coordinates are compressed vertically, and the vertical axis is stretched by the one factor that
makes the triangle's three sides, each holding fifteen data qubits, equal in length. The factor
is computed from the three corner qubits and comes out at two root three. The same factor, and
no other, puts every qubit at the same distance from its nearest neighbour.

The faint points are the 169 data qubits; the brighter ones, along the top edge, are the fifteen
whose final measurements make up the logical observable. The faint rings are the 84 auxiliary
qubits, where the checks are read. Each fault is drawn at the centre of the detectors it flips,
with a line to each of them. The three faults that flip the observable on their own are in red;
the other nine flip only detectors.

Before anything is drawn the script checks the certificate: that the twelve entries are the
listed mechanisms of the shipped mechanism file, that every detector is flipped an even number
of times, and that the observable is flipped an odd number of times. It refuses to draw
otherwise.

Nineteen detectors are touched, and every one of them by exactly two of the twelve faults, so
every one is flipped twice and reads nothing. All nineteen carry the same time coordinate in the
circuit, t = 1, so the drawing is a single plane. The chain starts on the top edge, where the
three red faults sit against the observable's row, and ends one lattice step from the bottom
corner.

doi:10.5281/zenodo.23792064
"""
import json, math, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = "/Users/kirt/Documents/reserch math/certify-repo"
CERTS = os.path.join(REPO, "colorcode-dcirc-certificates", "certs")
CIRCUIT, WITNESS, MECHS = "c_d15_r15.stim", "witness_d15_r15_w12.json", "m_d15_r15.json"

W, MARGIN = 2200, 150
BG, QUBIT, LOGICAL, FAULT, POINT, RED = (
    "#080b11", "#2a3240", "#7b8494", "#8f9bad", "#f3ecdd", "#ff6f5e")


def need(ok, why):
    """A failed check stops the script in every interpreter mode, including python -O."""
    if not ok:
        raise SystemExit("refusing to draw: " + why)


def coords(s):
    return [float(x) for x in re.findall(r"\(([^)]*)\)", s)[0].split(",")]


def read_circuit(path):
    """QUBIT_COORDS, DETECTOR coordinates with SHIFT_COORDS applied, the final data measurement,
    and the data qubits named by OBSERVABLE_INCLUDE."""
    qubit, detector, shift, n_meas, last_m, obs_recs = {}, [], [0.0] * 8, 0, None, []
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        name = s.split()[0].split("(")[0]
        if name == "QUBIT_COORDS":
            qubit[int(s.split(")")[1].split()[0])] = coords(s)
        elif name == "SHIFT_COORDS":
            for i, x in enumerate(coords(s)):
                shift[i] += x
        elif name == "DETECTOR":
            v = coords(s)
            detector.append([v[i] + shift[i] for i in range(len(v))])
        elif name in ("M", "MX", "MY", "MZ", "MR", "MRX"):
            targets = [int(t) for t in s.split()[1:]]
            last_m = (n_meas, targets)
            n_meas += len(targets)
        elif name == "OBSERVABLE_INCLUDE":
            obs_recs += [n_meas + int(t[4:-1]) for t in s.split()[1:]]
    start, data = last_m
    need(all(start <= r < start + len(data) for r in obs_recs),
         "the observable reads something other than the final data measurement")
    return qubit, detector, data, [data[r - start] for r in obs_recs]


def check(witness, mechanisms):
    """The claims the drawing rests on, from the certificate itself."""
    need(all(mechanisms[i] == witness["mechanisms"][j]
             for j, i in enumerate(witness["mechanism_indices"])),
         "the witness entries are not the listed mechanisms")
    need(len(witness["mechanisms"]) == witness["weight"] == 12, "the witness does not hold twelve")
    flips, obs = {}, 0
    for m in witness["mechanisms"]:
        for d in m["dets"]:
            flips[d] = flips.get(d, 0) + 1
        obs += m["obs"].count(0)
    need(all(c % 2 == 0 for c in flips.values()), "a detector would fire")
    need(obs % 2 == 1, "the observable would not flip")
    return flips


def nearest(points, k):
    return [min(math.hypot(x - a, k * (y - b)) for j, (a, b) in enumerate(points) if j != i)
            for i, (x, y) in enumerate(points)]


def main():
    qubit, detector, data, logical = read_circuit(os.path.join(CERTS, CIRCUIT))
    witness = json.load(open(os.path.join(CERTS, WITNESS)))
    mechanisms = json.load(open(os.path.join(CERTS, MECHS)))["mechanisms"]
    flips = check(witness, mechanisms)
    touched = sorted(flips)
    need(all(flips[d] == 2 for d in touched), "a touched detector is not touched exactly twice")
    t_layer = {detector[d][2] for d in touched}
    need(len(t_layer) == 1, "the touched detectors do not share one time coordinate")
    t = t_layer.pop()

    # The vertical factor that makes the triangle equilateral, from its three corner qubits.
    pts = [qubit[q][:2] for q in data]
    y0 = min(p[1] for p in pts)
    edge = [p for p in pts if p[1] == y0]
    left, right = min(edge), max(edge)
    far = max(pts, key=lambda p: p[1])
    side = right[0] - left[0]
    k = math.sqrt((side ** 2 - (far[0] - left[0]) ** 2) / (far[1] - y0) ** 2)
    sides = [side, math.hypot(far[0] - left[0], k * (far[1] - y0)),
             math.hypot(right[0] - far[0], k * (far[1] - y0))]
    need(max(sides) - min(sides) < 1e-9 * side, "the corners do not make an equilateral triangle")
    gaps = nearest([v[:2] for v in qubit.values()], k)
    need(max(gaps) - min(gaps) < 1e-9 * side, "the stretched lattice is not regular")

    s = (W - 2 * MARGIN) / side
    H = int(round(2 * MARGIN + k * (far[1] - y0) * s))
    X = lambda x: MARGIN + (x - left[0]) * s
    Y = lambda y: MARGIN + k * (y - y0) * s          # downward, as the circuit's coordinates run

    dset = set(data)
    sites = sorted(tuple(qubit[q][:2]) for q in qubit if q not in dset)
    need(len(sites) == 84 and set(sites) == {(c[0], c[1]) for c in detector if c[2] == t},
         "the auxiliary qubits are not the places where that layer's checks are read")
    out = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
           % (W, H, W, H), '<rect width="%d" height="%d" fill="%s"/>' % (W, H, BG)]
    for (x, y) in sites:
        out.append('<circle cx="%.1f" cy="%.1f" r="7.5" fill="none" stroke="%s" stroke-width="1.8"/>'
                   % (X(x), Y(y), QUBIT))
    lset = set(logical)
    for q in data:
        x, y = qubit[q][:2]
        out.append('<circle cx="%.1f" cy="%.1f" r="3.8" fill="%s"/>'
                   % (X(x), Y(y), LOGICAL if q in lset else QUBIT))
    faults = []
    for m in witness["mechanisms"]:
        ps = [(X(detector[d][0]), Y(detector[d][1])) for d in m["dets"]]
        c = (sum(p[0] for p in ps) / len(ps), sum(p[1] for p in ps) / len(ps))
        faults.append((c, ps, RED if 0 in m["obs"] else FAULT))
    for colour in (FAULT, RED):                        # red lines last, so they are never covered
        for (cx, cy), ps, col in faults:
            if col == colour:
                for (px, py) in ps:
                    out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                               'stroke-width="3.2" stroke-linecap="round"/>' % (cx, cy, px, py, col))
    for (cx, cy), ps, col in faults:                   # a dark halo, then the node
        out.append('<circle cx="%.1f" cy="%.1f" r="12.5" fill="%s"/>' % (cx, cy, BG))
        out.append('<circle cx="%.1f" cy="%.1f" r="8.5" fill="%s"/>' % (cx, cy, col))
    for d in touched:
        out.append('<circle cx="%.1f" cy="%.1f" r="5.5" fill="%s"/>'
                   % (X(detector[d][0]), Y(detector[d][1]), POINT))
    out.append("</svg>")
    open(os.path.join(HERE, "circuit.svg"), "w").write("\n".join(out) + "\n")
    print("ok: %d detectors, each touched twice, all at t = %g; vertical factor %.12f (two root "
          "three %.12f); nearest-neighbour distance %.9f for all %d qubits; %d logical qubits; "
          "%d red; canvas %d x %d" % (len(touched), t, k, 2 * math.sqrt(3), gaps[0], len(gaps),
                                      len(logical), sum(f[2] == RED for f in faults), W, H))


if __name__ == "__main__":
    main()
