"""Where the arithmetic stops being true.

Every finite-difference micromagnetic simulator evaluates the same closed form: Newell's
demagnetization tensor for a pair of uniformly magnetized cells. Done in double precision
it dies of catastrophic cancellation, and this program certified exactly where.

Each line is one geometry-and-component combination from the certified table — a cell shape,
an offset direction, a tensor entry — plotted as the number of correct decimal digits the
naive double-precision evaluation still holds, against separation. The digit counts are not
estimated: each is the rigorous comparison of the double against a two-sided rational
enclosure of the true value, computed in outward-rounded interval arithmetic.

The lines fall at about six digits per decade. Where a line crosses zero, that geometry has
no correct significant figure left; below the line the values are not merely inaccurate but
wrong in sign and order of magnitude. The crossing is not one number: it runs from about a
hundred cells for an elongated cell off-diagonal to roughly two thousand for a thin film's
out-of-plane entry. The familiar "about 300 cells" is the cube's crossing, not a constant.

Nothing is drawn but the data. The horizon is zero correct digits.

doi:10.5281/zenodo.21922469
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = "/Users/kirt/Documents/reserch math/certify-repo/demag-certificates/demag_certificate.json"
W, H = 2200, 1250
LEFT, RIGHT = 150, W - 110
TOP, BOT = 120, H - 150

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
DMIN = min(d for v in curves.values() for _, d in v)
SEPS = sorted({s for v in curves.values() for s, _ in v})

def X(sep):  return LEFT + (RIGHT - LEFT) * math.log10(sep) / math.log10(SEPS[-1])
def Y(dig):  return TOP + (BOT - TOP) * (DMAX - dig) / (DMAX - DMIN)

def hsv(h, s, v):
    h = (h % 360) / 60.0; i = int(h) % 6; f = h - int(h)
    p, q, t = v * (1 - s), v * (1 - s * f), v * (1 - s * (1 - f))
    a, b, c = ((v, t, p), (q, v, p), (p, v, t), (p, q, v), (t, p, v), (v, p, q))[i]
    return "#%02x%02x%02x" % (int(a * 255), int(b * 255), int(c * 255))

def crossing(v):
    """the separation at which this geometry loses its last correct digit"""
    for (s1, d1), (s2, d2) in zip(v, v[1:]):
        if d1 > 0 >= d2:
            f = d1 / (d1 - d2)
            return s1 * (s2 / s1) ** f
    return None

def emit(o, pts, col, below):
    d_ = "M%.1f %.1f" % pts[0] + "".join("L%.1f %.1f" % p for p in pts[1:])
    o.append('<path d="%s" fill="none" stroke="%s" stroke-width="%.1f" opacity="%.2f" '
             'stroke-linecap="round" stroke-linejoin="round"/>'
             % (d_, col, 1.0 if below else 2.0, 0.16 if below else 0.72))

def svg():
    o = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
         % (W, H, W, H),
         '<defs><linearGradient id="sky" x1="0" y1="0" x2="0" y2="1">'
         '<stop offset="0" stop-color="#0a1018"/><stop offset="%.2f" stop-color="#121a26"/>'
         '<stop offset="1" stop-color="#05070b"/></linearGradient></defs>'
         % (Y(0) / H),
         '<rect width="%d" height="%d" fill="url(#sky)"/>' % (W, H)]

    # decade rules, faint
    for k in range(0, 5):
        x = X(10 ** k)
        if LEFT <= x <= RIGHT:
            o.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#2a3646" '
                     'stroke-width="1" opacity="0.42"/>' % (x, TOP, x, BOT))

    order = sorted(curves.items(), key=lambda kv: crossing(kv[1]) or 1e9)
    n = len(order)
    for i, (key, v) in enumerate(order):
        c = crossing(v)
        t = i / (n - 1)                      # earliest death cool, latest warm
        col = hsv(196 - 168 * t, 0.50 + 0.22 * t, 0.80 + 0.18 * t)
        # split at the horizon: above it the digit count means something, below it the
        # value is wrong in sign and magnitude and any apparent recovery is coincidence
        run, below = [], None
        for s_, d in v:
            b = d <= 0
            if below is None: below = b
            if b != below:
                run.append((X(s_), Y(d)))
                emit(o, run, col, below)
                run = [(X(s_), Y(d))]; below = b
            else:
                run.append((X(s_), Y(d)))
        if len(run) > 1: emit(o, run, col, below)
        if c:                                 # mark the moment it loses the last digit
            o.append('<circle cx="%.1f" cy="%.1f" r="3.4" fill="%s" opacity="0.95"/>'
                     % (X(c), Y(0), col))

    # the horizon: zero correct digits
    o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#d8cfc0" stroke-width="1.6" '
             'opacity="0.72"/>' % (LEFT - 40, Y(0), RIGHT + 30, Y(0)))
    o.append("</svg>")
    return "".join(o)

if __name__ == "__main__":
    s = svg()
    open(os.path.join(HERE, "digits.svg"), "w").write(s)
    cs = sorted((crossing(v), k) for k, v in curves.items() if crossing(v))
    print("curves:", len(curves), "| separations per curve:", len(SEPS))
    print("digits at one cell: %.2f   worst in table: %.2f" % (DMAX, DMIN))
    print("first geometry to lose its last digit: %-38s at %6.0f cells" % (cs[0][1], cs[0][0]))
    print("last  geometry to lose its last digit: %-38s at %6.0f cells" % (cs[-1][1], cs[-1][0]))
    print("bytes:", len(s))
