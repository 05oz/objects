# Source

The drawings are generated from the certificate files in
[05oz/certify](https://github.com/05oz/certify) — the same JSON that the verification
checkers read. No coordinates are hand-placed: vertices go on a circle in index order and the
arcs follow the object's own adjacency, so the picture cannot disagree with the mathematics.

To redraw a piece, read the witness JSON (`N`, `arcs`), place the vertices on a circle, and draw
each arc as a quadratic curve whose control point is pulled toward the centre. The Paley blocks
are the vertices with exactly seven non-neighbours.
