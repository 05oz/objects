# Source

The drawings are generated from the certificate files in
[05oz/certify](https://github.com/05oz/certify) — the same JSON that the verification checkers
read. No coordinates are hand-placed, so a picture cannot disagree with the mathematics it is
drawn from.

## `forest.py`

Draws [`../pieces/forest.svg`](../pieces/forest.svg). Six published results decide the trees,
and nothing else does:

| what it decides | where it comes from |
| --- | --- |
| where a branch forks, and to what | the 13 rigid {I₃,TT₄}-free witnesses of k(3,4) = 21, read from `k34add-certificates/w*.json` |
| three internodes to a lineage | those objects contain no TT₄, so the longest transitive chain is 3; the constant is inline |
| which way a branch turns, and its hue | QR₇ = Cay(ℤ₇,{1,2,4}), applied as the Cayley difference (*u* − *v*) mod 7 computed from each arc |
| how many children survive at each depth | A_w, the exact uncorrectable-fault-set spectrum of the distance-5 rotated surface code, read from `wedge2-certificates/certificate_d5_r1_p1over100_exact.json` |
| whether a crown opens or closes | the 20 certified ZEFOZ points of ¹⁶⁷Er³⁺:Y₂SiO₅, read from `zefoz_sigs.txt` alongside this file |
| haze, saturation, and the loss of fine twigs | Newell's demagnetization tensor: 15.2 − 6.0 log₁₀ *n* correct digits, the fitted constants inline as `D0`, `SLOPE` |

`zefoz_sigs.txt` holds one line per certified ZEFOZ point: the first column is the
crystallographic site (1 or 2), the next two are the transition's hyperfine level pair, and
the fourth is the Hessian signature. The point itself — the field vector — lives in the
certificate (`certificate2.json`, field `B_mT`) as exact rationals and is not repeated here. Thirteen read `--+` (saddle) and seven `---`
(maximum).

### Redrawing

`REPO` at the top of the file is an absolute path to a local clone of
[05oz/certify](https://github.com/05oz/certify). Point it at your own clone and run:

```
python3 forest.py
```

It writes `forest.svg` beside itself — 111,132 branches — and the shipped copy is that file
moved to `../pieces/`. The output is deterministic: the same certificates give the same bytes.
The one arbitrary element, the scatter of trees across the ground, is a fixed-seed LCG written
out in the source.

## The spectral panels

[`../pieces/spectral_01.svg`](../pieces/spectral_01.svg) … `spectral_13.svg` are the thirteen
witnesses drawn one to a panel, and `../pieces/spectral_grid.svg` is those same thirteen panels
tiled five across, unchanged but for scale.

The layout is the object's own. The adjacency matrix is made skew-symmetric, its eigenvalues are
then purely imaginary, and a pair of eigenvectors gives the two-dimensional frame the vertices
are plotted in. Colour follows out-degree, on one scale across all thirteen: out-degree 5 is
blue, 6 amber, 7 green, and those are the only out-degrees that occur. An arc takes the colour
of the vertex it leaves. An object with a symmetry would collapse into a circle under this
treatment; these thirteen are rigid, and have none to collapse into.

## What does not ship here

Generators and searchers stay in the program. What is published is the certificate, the checker
that replays it, and — here — the drawing made from it.
