"""Twenty magnetic fields at which a spin forgets more slowly, and the ladders they make.

An erbium-167 ion in yttrium orthosilicate has one electron spin and one nuclear spin of 7/2,
so sixteen levels. Their energies move as the applied magnetic field moves, and at a few
particular fields a transition frequency stops moving to first order. Those are the ZEFOZ
points, and they are where a quantum memory built on this ion holds its coherence longest.

This program certified twenty of them, across both crystallographic sites of the host. Each
column here is one, drawn from its own certificate: the sixteen level energies at that field,
on a scale shared by all twenty so the columns can be compared. The two levels the clock
transition runs between are picked out and joined.

No energy here is a number. Each is a two-sided rational enclosure, and what the certificate
states is that the true level lies inside it. The enclosures are about 2 x 10^-10 MHz wide
against a span of some 100,000 MHz, a relative width near 10^-15, so at any scale a page can
hold they are narrower than the line used to draw them. They are drawn as lines because
nothing wider would be true.

The sixteen levels do not sit evenly. They fall into two bands of eight, and the void across
the middle of the figure is not empty space but the electron Zeeman splitting, which runs from
16,661 to 89,238 MHz across these twenty fields. Each band of eight is the nuclear manifold of
a spin 7/2. Every one of the twenty certified transitions lies inside a single band, none
crosses the gap, and their frequencies run from 736 to 2,220 MHz. So every ZEFOZ point in the
published atlas is a nuclear transition within one electron manifold, which is readable here
straight off the certificate and is why these are the long-lived ones.

The fields are exact rationals in millitesla and are not small, from 191 to 2,475. Columns are
ordered by field strength, which is why the ladders drift. The rule separates site 1 from site 2. Colour carries
one thing only, whether the point is a saddle of the transition frequency or a maximum, and it
follows the convention of the curvature piece, which draws these same twenty objects. Thirteen of the twenty are saddles of the transition frequency and seven are maxima;
none is a minimum, which is a result of the same certificate and not an assumption.

doi:10.5281/zenodo.21898996
"""
import json, math, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = "/Users/kirt/Documents/reserch math/certify-repo/zefoz-certificates/certificate2.json"
W, H = 2300, 1420
TOP, BOT = 150, H - 150

def num(x):
    return float(F(x)) if isinstance(x, str) else float(x)

def points():
    d = json.load(open(CERT))
    out = []
    for p in d["points"]:
        if p.get("kind") != "zefoz_point":
            continue
        lv = [((num(a) + num(b)) / 2, num(b) - num(a))
              for a, b in p["eigenvalue_brackets_MHz"]]
        B = [num(v) for v in p["B_mT"]]
        for pr in p["pairs"]:
            out.append(dict(site=p["site"], levels=lv, B=math.sqrt(sum(v * v for v in B)),
                            i=pr["i"], j=pr["j"],
                            sig=tuple(pr["hessian_signature"])))
    out.sort(key=lambda r: (r["site"], r["B"]))
    return out

def svg():
    P = points()
    n = len(P)
    lo = min(e for r in P for e, _ in r["levels"])
    hi = max(e for r in P for e, _ in r["levels"])
    pad = (hi - lo) * 0.04
    Y = lambda e: BOT - (BOT - TOP) * (e - (lo - pad)) / ((hi + pad) - (lo - pad))
    step = (W - 190) / n
    half = step * 0.40

    o = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
         % (W, H, W, H),
         '<defs><radialGradient id="bg" cx="50%" cy="46%" r="80%">'
         '<stop offset="0" stop-color="#0d1017"/><stop offset="1" stop-color="#05060a"/>'
         '</radialGradient></defs>',
         '<rect width="%d" height="%d" fill="url(#bg)"/>' % (W, H)]

    o.append('<line x1="70" y1="%.1f" x2="%d" y2="%.1f" stroke="#252c37" stroke-width="1" '
             'opacity="0.8"/>' % (Y(0), W - 70, Y(0)))          # zero of energy

    for k in range(1, n):                       # where site 1 ends and site 2 begins
        if P[k]["site"] != P[k - 1]["site"]:
            x = 120 + step * k
            o.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#2a3240" '
                     'stroke-width="1" opacity="0.75"/>' % (x, TOP - 20, x, BOT + 20))

    for k, r in enumerate(P):
        cx = 120 + step * (k + 0.5)
        sad = "+" in r["sig"]
        # the only hue that means anything: saddle or maximum, matching curvature.py,
        # which draws these same twenty points. Everything else is neutral.
        hue = "#e0708e" if sad else "#5fd0d8"
        for idx, (e, _) in enumerate(r["levels"]):
            on = idx in (r["i"], r["j"])
            o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                     'stroke-width="%.1f" opacity="%.2f" stroke-linecap="round"/>'
                     % (cx - half, Y(e), cx + half, Y(e),
                        hue if on else "#7b8494", 3.2 if on else 1.4, 0.95 if on else 0.34))
        ei, ej = r["levels"][r["i"]][0], r["levels"][r["j"]][0]
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                 'stroke-width="1.1" opacity="0.42"/>' % (cx, Y(ei), cx, Y(ej), hue))

    o.append('<text x="120" y="86" font-family="Georgia,serif" font-size="34" fill="#e8e3d6">'
             'Twenty certified ZEFOZ points of ¹⁶⁷Er³⁺:Y₂SiO₅'
             '</text>')
    o.append('<text x="120" y="%d" font-family="Georgia,serif" font-size="24" fill="#767d89">'
             'sixteen levels each &#183; site 1 left of the rule, site 2 right, each ordered by '
             'field strength &#183; saddle rose, maximum teal</text>' % (H - 76))
    o.append("</svg>")
    return "".join(o), P

if __name__ == "__main__":
    s, P = svg()
    open(os.path.join(HERE, "levels.svg"), "w").write(s)
    wid = [w for r in P for _, w in r["levels"]]
    span = max(e for r in P for e, _ in r["levels"]) - min(e for r in P for e, _ in r["levels"])
    sad = sum(1 for r in P if "+" in r["sig"])
    print("  points:", len(P), "| levels each:", len(P[0]["levels"]),
          "| sites:", sorted({r["site"] for r in P}))
    print("  enclosure width: %.3e to %.3e MHz | span %.0f MHz | relative %.1e"
          % (min(wid), max(wid), span, max(wid) / span))
    print("  field strength: %.1f to %.1f mT" % (min(r["B"] for r in P), max(r["B"] for r in P)))
    print("  saddles:", sad, "| maxima:", len(P) - sad, "| minima:",
          sum(1 for r in P if r["sig"] == ("+", "+", "+")))
    print("  bytes:", len(s))
