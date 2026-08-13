"""CFR(5,25): five permutations that cannot all be regular.

A circular Florentine rectangle of order 25 with 5 rows. Each row is a permutation
of Z_25, and the defining condition is that for any two distinct symbols and any
distance, at most one row places the second that far after the first. F_c(25) >= 5
was established by exhibiting this object; the Handbook of Combinatorial Designs,
2nd ed., recorded 4.

Each row is drawn as the closed path that visits the 25 symbols in the order the row
lists them, on 25 points evenly spaced by symbol. Nothing is arranged: the circle is
the symbol set, and the path is the row read left to right.

What the drawing shows is a fact about the object. Two rows are linear maps
j -> cj mod 25, with c = 1 and c = 7, so they close as the regular 25-gon and the
{25/7} star. The remaining three are not linear and cannot be: five rows all of that
form would repeat an ordered pair at some distance and break the condition. Their
irregularity is what makes five rows possible at all. Row 2 keeps a trace of order —
its successive steps are 4, 9, 14, 19, 24, every one congruent to 4 mod 5.

doi:10.5281/zenodo.21831896
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = "/Users/kirt/Documents/reserch math/certify-repo/cfr-certificates/CFR_5_25.json"
S = 1500
C = S / 2
R = S * 0.385

d = json.load(open(CERT))
rows, n = d["rows"], d["n"]

def pt(sym):
    a = 2 * math.pi * sym / n - math.pi / 2
    return C + R * math.cos(a), C + R * math.sin(a)

def linear(r):
    """the multiplier c if the row is j -> cj mod n, else None"""
    for c in range(1, n):
        if all(r[j] == (c * j) % n for j in range(n)):
            return c
    return None

def hsv(h, s, v):
    h = (h % 360) / 60.0; i = int(h) % 6; f = h - int(h)
    p, q, t = v * (1 - s), v * (1 - s * f), v * (1 - s * (1 - f))
    a, b, c = ((v, t, p), (q, v, p), (p, v, t), (p, q, v), (t, p, v), (v, p, q))[i]
    return "#%02x%02x%02x" % (int(a * 255), int(b * 255), int(c * 255))

def path(r, bow):
    """the row as a closed path; chords bow toward the centre so overlapping
    rows stay legible and long chords read as arcs rather than clutter"""
    segs = []
    for k in range(n):
        a, b = pt(r[k]), pt(r[(k + 1) % n])
        mx, my = (a[0] + b[0]) / 2, (a[1] + b[1]) / 2
        vx, vy = C - mx, C - my
        segs.append("M%.1f %.1fQ%.1f %.1f %.1f %.1f"
                    % (a[0], a[1], mx + vx * bow, my + vy * bow, b[0], b[1]))
    return "".join(segs)

def svg():
    o = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
         % (S, S, S, S),
         '<defs><radialGradient id="bg" cx="50%" cy="50%" r="72%">'
         '<stop offset="0" stop-color="#0e1119"/><stop offset="1" stop-color="#05070b"/>'
         '</radialGradient></defs>',
         '<rect width="%d" height="%d" fill="url(#bg)"/>' % (S, S)]

    # irregular rows first, regular ones over them: the stars are the readable frame
    order = sorted(range(len(rows)), key=lambda i: linear(rows[i]) is not None)
    for i in order:
        r = rows[i]
        c = linear(r)
        # hue by row; the two linear rows run cool, the three that cannot be linear run warm
        hue = 200 - 14 * i if c else 44 - 16 * (i - 2)
        col = hsv(hue, 0.30 if c else 0.66, 0.94 if c else 0.95)
        # a linear row's chords are all one length, so a light bow suffices;
        # the irregular rows carry mixed lengths and need more to separate
        bow = 0.10 if c else 0.20 + 0.05 * i
        o.append('<path d="%s" fill="none" stroke="%s" stroke-width="%.1f" opacity="%.2f" '
                 'stroke-linecap="round"/>'
                 % (path(r, bow), col, 2.2 if c else 1.8, 0.70 if c else 0.72))

    for s in range(n):
        x, y = pt(s)
        o.append('<circle cx="%.1f" cy="%.1f" r="5.4" fill="#efe9dc" opacity="0.93"/>' % (x, y))
    o.append("</svg>")
    return "".join(o)

if __name__ == "__main__":
    s = svg()
    open(os.path.join(HERE, "florentine.svg"), "w").write(s)
    for i, r in enumerate(rows):
        c = linear(r)
        steps = sorted({(r[(j + 1) % n] - r[j]) % n for j in range(n)})
        print(f"  row {i}: {'linear c=%d' % c if c else 'not linear'}, "
              f"{len(steps)} distinct steps {steps if len(steps) <= 7 else ''}")
    print("bytes:", len(s))
