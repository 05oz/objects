# Objects

Pictures of mathematical objects that this program proved exist.

Nothing here is decorative. Every line is the object itself — the vertices, the arcs, the
symmetries are what the theorem says they are, drawn without embellishment. If a picture looks
striking, that is a property of the object, not of the drawing.

Each piece names the result it comes from and the archived record where the object and its
verification certificate live. The mathematics is at
[github.com/05oz/certify](https://github.com/05oz/certify); the program is
[halfounce.io](https://halfounce.io).

---

## Thirteen rigid objects, each placing its own points

![spectral](pieces/spectral_grid.png)

The oriented Ramsey number k(3,4) is 21. On 20 vertices there are graphs containing neither three
mutually non-adjacent vertices nor four in transitive order — and there are at least thirteen of
them, no two alike. Every one is *rigid*: only the identity permutation maps it to itself.

Nothing here was laid out by hand. Each object's adjacency matrix, made skew-symmetric, has
purely imaginary eigenvalues whose eigenvectors give a two-dimensional frame; the vertices sit
where that frame puts them. The layout is computed from the object, not chosen for it. Colour
follows out-degree. A symmetric object would collapse into a circle under this treatment — these
cannot, because they have no symmetry to collapse into, and that is why no two of the thirteen
share a silhouette.

These pictures could not have been made before August 2026, because the objects were not known
to exist.

`doi:10.5281/zenodo.21890619`

## QR₇ — the Paley tournament on seven vertices

![paley](pieces/paley_qr7.png)

Seven points; an arrow from *i* to *j* when *j − i* is a perfect square modulo 7 — that is, one
of {1, 2, 4}. It is the only tournament on seven vertices with no transitively ordered
quadruple, and it appears inside every one of the thirteen witnesses above. The pattern was not
chosen; it is forced by the constraint.

## Cay(ℤ₂₈, {3, 8, 10, 12, 17})

![circulant](pieces/circulant_k63.png)

Twenty-eight points on a circle, each joined to the points 3, 8, 10, 12 and 17 steps ahead. The
connection set is *sum-free* — no two of its elements sum to a third — which is exactly what
makes the resulting tournament free of transitive triangles. It has independence number 5, and
so witnesses that k(6,3) ≥ 29.

The rosette and the void at the centre are not design choices. They are what those five
distances do when drawn.

`doi:10.5281/zenodo.21890619`

---

## Prints

`print/` holds 12 × 12 inch masters at 300 dpi, each with a caption plate stating the object,
the theorem it witnesses, and the DOI of its archived certificate. The `.svg` sources are
resolution-independent and can be enlarged without loss.

## Licence

Images: CC BY 4.0 — use them, credit *Half Ounce Research* and the DOI on the piece.
Source code in `src/`: Apache-2.0.

The objects themselves belong to nobody.
