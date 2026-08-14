"""The shape of a machine proof, and how much of it the conclusion never uses.

Two refutations, both checked line by line by a standard-library program. Each proves that a
covering number cannot be as large as was asked, and together they determine nu_3(9) = 9 and
nu_3(10) = 12. A refutation is not a narrative: it is a directed acyclic graph. Every lemma
names the clauses it was derived from, and the argument ends when the empty clause appears.

Each lemma is a point. Its horizontal position is its derivation depth, the length of the
longest chain of inferences standing behind it, which runs to 472 in the first proof and 2,225
in the second. Its vertical position is the width of the clause it states, from a single literal
up to 24 and 63. Nothing is arranged; both coordinates are read out of the proof file.

The bright points are the lemmas the empty clause actually depends on, found by walking the
antecedent graph backwards from the contradiction. The faint ones are not. A solver derives what
it needs and a great deal that it does not, and the proportion is the thing worth seeing: of the
4,530 lemmas in the first proof only 1,750 stand behind the conclusion, so 61% of it is
scaffolding. In the second, 7,971 of 12,067 are load-bearing and 34% is not. The checker verifies
all of it regardless, because a proof is checked as written and not as needed.

The two kinds of lemma are not alike, and the difference is sharper than the picture alone
shows. In both proofs the load-bearing lemmas cite about 29 earlier clauses apiece; the discarded
ones cite between 4 and 7. What survives into the argument is what the solver worked hardest to
derive, and cheap inferences are mostly the ones it turns out not to need.

Each bright lemma is joined to its deepest antecedent, one line per lemma rather than the twenty
or so it cites, so what shows is the spine of the derivation and not the full thicket. Everything
converges on a single point at the lower right: width zero, maximum depth, the empty clause.

doi:10.5281/zenodo.21816010
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CDIR = "/Users/kirt/Documents/reserch math/certify-repo/tt3-certificates"
PANELS = [("min9_ge10.lrat", "ν₃(9) = 9", "no 10 blocks on 9 points"),
          ("min10_ge13.lrat", "ν₃(10) = 12", "no 13 blocks on 10 points")]
PW, PH = 1560, 860
W, H = PW, PH * len(PANELS)
L, R, T, B = 130, PW - 90, 96, PH - 104

def parse(fn):
    """lemma id -> (antecedents, clause width), in file order"""
    ants, width, order = {}, {}, []
    for line in open(os.path.join(CDIR, fn), errors="ignore"):
        p = line.split()
        if len(p) < 2 or p[1] == "d":
            continue
        i = int(p[0])
        z = p.index("0", 1)
        width[i] = z - 1
        rest = p[z + 1:]
        z2 = rest.index("0") if "0" in rest else len(rest)
        ants[i] = [abs(int(x)) for x in rest[:z2]]
        order.append(i)
    return ants, width, order

def analyse(ants, width, order):
    empty = [i for i in order if width[i] == 0]
    seen, stack = set(), list(empty)            # walk back from the contradiction
    while stack:
        c = stack.pop()
        if c in seen:
            continue
        seen.add(c)
        stack.extend(ants.get(c, ()))
    depth = {}
    for i in order:
        a = [depth.get(x, 0) for x in ants[i]]
        depth[i] = 1 + max(a) if a else 1
    return depth, {c for c in seen if c in ants}, (empty[0] if empty else None)

def panel(fn, title, sub, oy):
    ants, width, order = parse(fn)
    depth, core, empty = analyse(ants, width, order)
    dmax = max(depth.values()); wmax = max(width.values())
    X = lambda d: L + (R - L) * (d / dmax)
    Y = lambda w: oy + B - (B - T) * (w / wmax)

    o = []
    # the spine: each load-bearing lemma to its deepest antecedent, one line apiece
    for i in order:
        if i not in core or not ants[i]:
            continue
        a = max(ants[i], key=lambda x: depth.get(x, 0))
        if a not in depth:
            continue
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#c98a3a" '
                 'stroke-width="0.55" opacity="0.16"/>'
                 % (X(depth[a]), Y(width[a]), X(depth[i]), Y(width[i])))
    # the lemmas the conclusion never uses
    for i in order:
        if i in core:
            continue
        o.append('<circle cx="%.1f" cy="%.1f" r="1.9" fill="#4f7699" opacity="0.60"/>'
                 % (X(depth[i]), Y(width[i])))
    # the lemmas it does
    for i in order:
        if i not in core:
            continue
        o.append('<circle cx="%.1f" cy="%.1f" r="1.7" fill="#ffc978" opacity="0.66"/>'
                 % (X(depth[i]), Y(width[i])))
    if empty is not None:
        x, y = X(depth[empty]), Y(0)
        o.append('<circle cx="%.1f" cy="%.1f" r="26" fill="#fff0d0" opacity="0.07"/>' % (x, y))
        o.append('<circle cx="%.1f" cy="%.1f" r="13" fill="#fff0d0" opacity="0.13"/>' % (x, y))
        o.append('<circle cx="%.1f" cy="%.1f" r="5.4" fill="#fff4de" opacity="0.97"/>' % (x, y))

    o.append('<text x="%d" y="%.1f" font-family="Georgia,serif" font-size="38" fill="#ece6d8">'
             '%s</text>' % (L, oy + 58, title))
    o.append('<text x="%d" y="%.1f" font-family="Georgia,serif" font-size="25" fill="#79808c">'
             '%s &#183; %d lemmas, %d behind the contradiction, depth to %d</text>'
             % (L + 235, oy + 58, sub, len(order), len(core), dmax))
    return "".join(o), (len(order), len(core), dmax, wmax)

if __name__ == "__main__":
    body, meta = [], []
    for k, (fn, title, sub) in enumerate(PANELS):
        s, m = panel(fn, title, sub, k * PH)
        body.append(s); meta.append((fn,) + m)
    svg = ('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
           % (W, H, W, H)
           + '<rect width="%d" height="%d" fill="#080a10"/>' % (W, H)
           + "".join(body) + "</svg>")
    open(os.path.join(HERE, "proof.svg"), "w").write(svg)
    for fn, n, c, d, w in meta:
        print(f"  {fn}: {n} lemmas, {c} load-bearing ({100*c/n:.1f}%), "
              f"{n-c} scaffolding, depth to {d}, width to {w}")
    print("  bytes:", len(svg))
