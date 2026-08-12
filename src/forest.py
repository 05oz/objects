"""A forest in which every growth rule is a published result.

  topology   the 13 rigid {I3,TT4}-free witnesses of k(3,4) = 21. A branch sitting at
             vertex v forks to out-neighbours of v. 13 objects x 20 roots = 260 trees,
             and because the objects are rigid no two of them repeat.
  internode  TT4-freeness bounds the longest transitive chain at 3, so a lineage runs
             exactly three internodes before it must fork. That is the trees' rhythm.
  angle      QR7 = Cay(Z7,{1,2,4}), forced as the non-neighbourhood of any vertex with
             seven non-neighbours. An arc v->u carries the Cayley difference
             (u - v) mod 7; the branch turns toward the residues {1,2,4} and away from
             the non-residues {3,5,6}. That set is not symmetric, so the trees lean.
  density    A_w, the exact uncorrectable-fault-set spectrum of the distance-5 rotated
             surface code. Normalised log-ratios A_{w+1}/A_w over the weights that carry
             the sub-threshold probability set how many children survive at each depth.
  crown      the 20 certified ZEFOZ points: 13 saddles (-,-,+) open and spread,
             7 maxima (-,-,-) close into a dome.
  distance   Newell. Correct digits = 15.2 - 6.0 log10(n): about six lost per decade,
             none left near 300 cells. That one curve is the haze, the colour, and the
             loss of fine twigs. A tree far enough away has no significant figures left,
             and vanishes.

Only the scatter of the trees across the ground is arbitrary. Every tree's own form is
computed.
"""
import json, glob, math, os

REPO = "/Users/kirt/Documents/reserch math/certify-repo"
HERE = os.path.dirname(os.path.abspath(__file__))
W, H = 1584, 396
HORIZON = 302.0
NEAR, FAR = 1.0, 300.0          # Newell cell separations: the certified regime
D0, SLOPE = 15.2, 6.0           # digits at one cell; digits lost per decade

WIT = [json.load(open(f)) for f in sorted(glob.glob(f"{REPO}/k34add-certificates/w*.json"))]
SPEC = json.load(open(f"{REPO}/wedge2-certificates/certificate_d5_r1_p1over100_exact.json"))
AW = {int(k): int(v) for k, v in SPEC["A_w_full_spectrum"].items() if int(v) > 0}
SIG = [ln.split()[3] for ln in open(os.path.join(HERE, "zefoz_sigs.txt"))]

QR7 = {1, 2, 4}

def out_neighbours(w):
    o = [[] for _ in range(w["N"])]
    for u, v in w["arcs"]:
        o[u].append(v)
    return o

OUT = [out_neighbours(w) for w in WIT]

_ws = sorted(AW)[:11]                       # w = 3..13, where the probability lives
_ratio = [AW[b] / AW[a] for a, b in zip(_ws, _ws[1:])]
_lg = [math.log(r) for r in _ratio]
_lo, _hi = min(_lg), max(_lg)
BRANCH = [int(round(2 + (x - _lo) / (_hi - _lo))) for x in _lg]

def digits(n):
    return max(0.0, D0 - SLOPE * math.log10(n))

# --------------------------------------------------------------------- growth
def grow(segs, x, y, ang, L, wid, v, out, depth, maxd, spread, curve, decay, dd=-1):
    if depth > maxd or L < 0.55 or wid < 0.03:
        return
    px, py, a = x, y, ang
    pts = [(px, py)]
    c = curve * (0.74 ** depth)
    for _ in range(3):                               # TT4-free: three internodes
        a += c * ((v % 7) - 3.0) / 3.0
        px += math.sin(a) * L / 3.0
        py -= math.cos(a) * L / 3.0
        pts.append((px, py))
    segs.append((pts, wid, depth, dd))
    nb = out[v]
    if not nb:
        return
    k = min(BRANCH[min(depth, len(BRANCH) - 1)], len(nb))
    st = (v * 3 + depth * 5) % len(nb)               # walk the object, don't re-read it
    lean = spread * (0.44 + 0.56 * 0.62 ** depth)
    for j in range(k):
        u = nb[(st + j) % len(nb)]
        dd = (u - v) % 7                             # the Cayley difference of the arc
        off = (1.0 if dd in QR7 else -1.0) * ((dd + 1) / 7.0) * lean
        grow(segs, px, py, a + off, L * decay * (0.93 + 0.03 * dd),
             wid * 0.66, u, out, depth + 1, maxd, spread, curve, decay, dd)

def tree(wi, root, height, n):
    d = digits(n)
    maxd = max(2, min(8, int(1 + round(d / 1.70))))
    saddle = SIG[root % 20] == "--+"
    spread = 0.76 if saddle else 0.52
    decay = 0.735 if saddle else 0.705
    curve = 0.075 if saddle else 0.115
    segs = []
    grow(segs, 0.0, 0.0, 0.0, height * 0.44, max(0.45, height * 0.0128),
         root, OUT[wi], 0, maxd, spread, curve, decay)
    return segs

# ---------------------------------------------------------------- composition
def lcg(seed):
    s = seed & 0xFFFFFFFF
    while True:
        s = (1103515245 * s + 12345) & 0x7FFFFFFF
        yield s / 0x7FFFFFFF

GROVES = [-118, 22, 158, 268, 372, 486, 596, 700, 812, 928, 1032, 1146,
          1258, 1372, 1478, 1596, 1706]

def build():
    rnd = lcg(20260812)
    NT = 720
    stand = []
    for i in range(NT):
        u = i / (NT - 1.0)
        n = math.exp(math.log(FAR) * (u ** 0.42))
        wi = i % 13
        root = (i * 7 + i // 13) % 20
        od = len(OUT[wi][root])
        vary = 0.70 + 0.50 * max(0.0, min(1.0, (od - 3.0) / 6.0))
        h = min(278.0, 266.0 * n ** -0.46 * vary)
        by = HORIZON + (H - HORIZON) * n ** -0.55
        g = GROVES[(i * 5 + i // 7) % len(GROVES)]
        x = g + (next(rnd) - 0.5) * 178.0
        # half the stand is grown from the converse digraph D^op, which is also
        # {I3,TT4}-free and sends every residue to a non-residue: the exact mirror
        mir = 1.0 if od % 2 == 0 else -1.0
        stand.append((n, x, by, h, wi, root, mir))
    stand.sort(key=lambda t: -t[0])

    buckets, total = {}, 0
    for n, x, by, h, wi, root, mir in stand:
        t = digits(n) / D0
        if t <= 0.010:
            continue
        op = round(min(0.95, 0.07 + 0.88 * t ** 0.95), 2)
        tq = round(t, 1)
        nbr = OUT[wi][root]
        bq = round(sum(1 for u in nbr if (u - root) % 7 in QR7) / max(1, len(nbr)), 1)
        for pts, wid, dep, dd in tree(wi, root, h, n):
            total += 1
            wb = round(min(3.4, max(0.16, wid)), 1)
            p = pts
            buckets.setdefault((ink(dd, min(dep, 6), tq, bq), op, wb), []).append(
                "M%.1f %.1fQ%.1f %.1f %.1f %.1fT%.1f %.1f" % (
                    x + mir * p[0][0], by + p[0][1], x + mir * p[1][0], by + p[1][1],
                    x + mir * p[2][0], by + p[2][1], x + mir * p[3][0], by + p[3][1]))
    return buckets, total

# Hue carries the Cayley difference of the arc that grew the branch, and the
# tree's own hue comes from the residue fraction of its root's out-arcs. Both are
# rotations, never blends — mixing two opposite hues through RGB only makes grey.
HAZE = (0x11, 0x18, 0x22)          # far trees dissolve into the dark, not out of it
H_COOL, H_WARM = 208.0, 30.0       # b = 0 runs blue, b = 1 runs amber, green between

def hsv(h, sat, val):
    h = (h % 360.0) / 60.0
    i = int(h) % 6
    f = h - int(h)
    p, q, t_ = val * (1 - sat), val * (1 - sat * f), val * (1 - sat * (1 - f))
    r, g, b = ((val, t_, p), (q, val, p), (p, val, t_),
               (p, q, val), (t_, p, val), (val, p, q))[i]
    return (r * 255.0, g * 255.0, b * 255.0)

def ink(dd, dep, t, b):
    """b is the residue fraction of the root's own out-arcs: the tree's hue."""
    b = min(1.0, max(0.0, 0.5 + (b - 0.5) * 2.2))   # push off the middle: fewer pure greens
    h = H_COOL + (H_WARM - H_COOL) * b
    if dd >= 0:                                  # the branch turns its own way off that
        h += (1.0 if dd in QR7 else -1.0) * (9.0 + 6.0 * dd)
    f = min(1.0, dep / 1.30)                     # bark is barely tinted, twigs are not
    hh = h % 360.0
    green = max(0.0, 1.0 - abs(hh - 118.0) / 46.0)   # greens go deep, never neon
    c = hsv(h, 0.13 + 0.44 * f, (0.79 - 0.08 * f) * (1.0 - 0.28 * green))
    sc = t ** 0.70                               # and the haze law washes it back out
    c = [HAZE[i] + (c[i] - HAZE[i]) * sc for i in range(3)]
    return "#%02x%02x%02x" % tuple(int(max(0, min(255, v))) for v in c)

def svg():
    buckets, total = build()
    o = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
         % (W, H, W, H),
         '<defs><linearGradient id="sky" x1="0" y1="0" x2="0" y2="1">'
         '<stop offset="0" stop-color="#05070b"/><stop offset="0.62" stop-color="#0d141d"/>'
         '<stop offset="0.78" stop-color="#16202c"/><stop offset="1" stop-color="#0a0f16"/>'
         '</linearGradient></defs>',
         '<rect width="%d" height="%d" fill="url(#sky)"/>' % (W, H)]

    g, rnd = [], lcg(77)
    for _ in range(430):
        uu = next(rnd)
        n = math.exp(math.log(FAR) * (uu ** 0.45))
        yy = HORIZON + (H - HORIZON) * n ** -0.55
        xx = -60 + next(rnd) * (W + 120)
        dn = digits(n) / D0
        g.append('<path d="M%.0f %.1fh%.0f" stroke="#8fa6bd" stroke-width="%.2f" opacity="%.3f"/>'
                 % (xx, yy, 5 + 92 * dn ** 1.6, 0.45 + 0.9 * dn,
                    round(0.025 + 0.20 * dn, 3)))
    o.append('<g fill="none" stroke-linecap="round">' + "".join(g) + "</g>")

    for key in sorted(buckets, key=lambda k: (k[1], k[2], k[0])):
        col, op, wb = key
        o.append('<g fill="none" stroke="%s" stroke-width="%.2f" opacity="%.2f" '
                 'stroke-linecap="round"><path d="%s"/></g>'
                 % (col, wb, op, "".join(buckets[key])))

    for bx, by, sc in ((452, 54, 7.4), (500, 71, 5.6), (534, 49, 4.3)):
        o.append('<path d="M%.1f %.1fq%.1f -%.1f %.1f 0q%.1f -%.1f %.1f 0" fill="none" '
                 'stroke="#b9c6d4" stroke-width="1.05" opacity="0.5" stroke-linecap="round"/>'
                 % (bx, by, sc * .5, sc * .42, sc, sc * .5, sc * .42, sc))
    o.append("</svg>")
    return "".join(o), total

if __name__ == "__main__":
    s, n = svg()
    open(os.path.join(HERE, "forest.svg"), "w").write(s)
    print("branches %d   bytes %d   branching %s" % (n, len(s), BRANCH[:9]))
