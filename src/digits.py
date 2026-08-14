"""Where the arithmetic stops being true.

Every finite-difference micromagnetic simulator evaluates the same closed form: Newell's
demagnetization tensor for a pair of uniformly magnetized cells. Done in double precision
it dies of catastrophic cancellation, and this program certified exactly where.

Each line is one geometry-and-component combination from the certified table — a cell shape,
an offset direction, a tensor entry — plotted as the number of correct decimal digits the
naive double-precision evaluation still holds, against separation. The digit counts are not
estimated: each is the rigorous comparison of the double against a two-sided rational
enclosure of the true value, computed in outward-rounded interval arithmetic.

The lines fall at about six digits per decade and each one is drawn brighter the longer it
survives. Where a line meets the horizon that geometry has no correct significant figure
left. The crossing is not one number. All fifty cross, and each crossing is located between
two tabulated separations rather than measured at one: they run from between 55 and 100 cells
for an elongated cell's off-diagonal entry to between 1000 and 2000 for a thin film's
out-of-plane entry. The familiar "about 300 cells" is where the cube's own crossings fall,
not a constant. The marker on the horizon is placed by logarithmic interpolation inside each
bracket, which is a drawing convenience and not a measured value.

Past the horizon the lines are drawn dissolving over the next five and a half digits,
following the running minimum. Those are the two rendering decisions here and both are
factual: a value with fewer than zero correct
digits is wrong in sign and order of magnitude, so its continued descent is not information.
The certificate records those numbers; the drawing declines to give them a shape.

doi:10.5281/zenodo.21922469
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = "/Users/kirt/Documents/reserch math/certify-repo/demag-certificates/demag_certificate.json"
W, H = 2400, 1330
LEFT, RIGHT = 205, W - 175
TOP, BOT = 165, H - 165
FLOOR = 5.5                       # digits below the horizon over which a line dissolves

entries = json.load(open(CERT))["entries"]

curves = {}
for e in entries:
    if e["sep_cells"] == 0:
        continue
    d = e["naive_digits"]
    d = d["correct_digits_lo"] if isinstance(d, dict) else d
    curves.setdefault((e["cell"], e["direction"], e["component"]), []).append((e["sep_cells"], d))
for k in curves:
    curves[k].sort()

DMAX = max(d for v in curves.values() for _, d in v)
SEPS = sorted({s for v in curves.values() for s, _ in v})

def X(sep):  return LEFT + (RIGHT - LEFT) * math.log10(sep) / math.log10(SEPS[-1])
def Y(dig):  return TOP + (BOT - TOP) * (DMAX - dig) / (DMAX + FLOOR)

def hsv(h, s, v):
    h = (h % 360) / 60.0; i = int(h) % 6; f = h - int(h)
    p, q, t = v * (1 - s), v * (1 - s * f), v * (1 - s * (1 - f))
    a, b, c = ((v, t, p), (q, v, p), (p, v, t), (p, q, v), (t, p, v), (v, p, q))[i]
    return "#%02x%02x%02x" % (int(a * 255), int(b * 255), int(c * 255))

def crossing(v):
    """the separation at which this geometry loses its last correct digit"""
    for (s1, d1), (s2, d2) in zip(v, v[1:]):
        if d1 > 0 >= d2:
            return s1 * (s2 / s1) ** (d1 / (d1 - d2))
    return None

def split(v):
    """the run that still carries a correct digit, and the run past the horizon"""
    above, below = [], []
    prev = None
    for s, d in v:
        if below:
            below.append((s, d))
        elif d > 0:
            above.append((s, d)); prev = (s, d)
        else:
            s1, d1 = prev
            sc = s1 * (s / s1) ** (d1 / (d1 - d))
            above.append((sc, 0.0))
            below = [(sc, 0.0), (s, d)]
    # Past the horizon the tail follows the running minimum. A count that has reached zero
    # has not recovered when the next separation reports a larger one; that is coincidence
    # in a number already wrong in sign, and drawing it as a rebound would assert otherwise.
    run, mono = 0.0, []
    for s, d in below:
        run = min(run, d)
        mono.append((s, run))
    return above, mono

def line(pts, col, wid, op):
    return ('<path d="M%.1f %.1f%s" fill="none" stroke="%s" stroke-width="%.2f" '
            'opacity="%.3f" stroke-linecap="round" stroke-linejoin="round"/>'
            % (pts[0][0], pts[0][1], "".join("L%.1f %.1f" % p for p in pts[1:]), col, wid, op))

def dissolve(o, below, col, op0):
    """past the horizon, taper to nothing over FLOOR digits — the value has no shape left"""
    for (s1, d1), (s2, d2) in zip(below, below[1:]):
        l1, l2 = math.log10(s1), math.log10(s2)
        for k in range(12):
            u, w = k / 12, (k + 1) / 12
            da, db = d1 + (d2 - d1) * u, d1 + (d2 - d1) * w
            m = (da + db) / 2
            if m <= -FLOOR:
                continue
            f = max(0.0, min(1.0, 1 + m / FLOOR)) ** 2.8
            o.append(line([(X(10 ** (l1 + (l2 - l1) * u)), Y(da)),
                           (X(10 ** (l1 + (l2 - l1) * w)), Y(db))],
                          col, 0.70 + 1.20 * f, op0 * f))

def svg():
    o = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
         % (W, H, W, H),
         '<rect width="%d" height="%d" fill="#07080c"/>' % (W, H)]

    order = sorted(curves.items(), key=lambda kv: crossing(kv[1]))
    n = len(order)
    marks = []
    for i, (key, v) in enumerate(order):
        t = i / (n - 1)
        # one hue family, ember to pale gold: the longer a geometry survives, the brighter
        # it is drawn. Rank sets brightness only; nothing here is a second variable.
        col = hsv(17 + 35 * t, 0.68 - 0.36 * t, 0.57 + 0.40 * t)
        op = 0.60 + 0.26 * t
        above, below = split(v)
        if len(above) > 1:
            o.append(line([(X(s), Y(d)) for s, d in above], col, 1.9, op))
        if below:
            dissolve(o, below, col, op)
            marks.append((X(below[0][0]), col, op))

    # the horizon is a level, so it spans; it only softens at either end rather than
    # stopping against the frame
    a, b = LEFT - 55, max(m[0] for m in marks) + 235
    o.insert(1, '<defs><linearGradient id="h" gradientUnits="userSpaceOnUse" '
                'x1="%.1f" y1="0" x2="%.1f" y2="0">' % (a, b) +
                '<stop offset="0" stop-color="#b3ab9c" stop-opacity="0"/>'
                '<stop offset="0.16" stop-color="#b3ab9c" stop-opacity="0.44"/>'
                '<stop offset="0.87" stop-color="#b3ab9c" stop-opacity="0.44"/>'
                '<stop offset="1" stop-color="#b3ab9c" stop-opacity="0"/></linearGradient></defs>')
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="url(#h)" stroke-width="1.3"/>'
             % (a, Y(0), b, Y(0)))
    for x, col, op in marks:                          # the last correct digit, geometry by geometry
        o.append('<circle cx="%.1f" cy="%.1f" r="13" fill="%s" opacity="%.3f"/>' % (x, Y(0), col, op * 0.12))
        o.append('<circle cx="%.1f" cy="%.1f" r="3.8" fill="%s" opacity="%.3f"/>' % (x, Y(0), col, min(1.0, op + 0.20)))
    o.append("</svg>")
    return "".join(o)

if __name__ == "__main__":
    s = svg()
    open(os.path.join(HERE, "digits.svg"), "w").write(s)
    cs = sorted((crossing(v), k) for k, v in curves.items())
    print("curves:", len(curves), "| separations per curve:", len(SEPS))
    print("digits at one cell: %.2f" % DMAX)
    print("first geometry to lose its last digit: %-38s at %6.0f cells" % (cs[0][1], cs[0][0]))
    print("last  geometry to lose its last digit: %-38s at %6.0f cells" % (cs[-1][1], cs[-1][0]))
    print("bytes:", len(s))
