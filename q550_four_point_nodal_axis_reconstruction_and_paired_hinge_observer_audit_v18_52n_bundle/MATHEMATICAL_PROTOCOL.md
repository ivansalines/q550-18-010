# Mathematical protocol — v18.52n

Let q(t) = (Q1,Q2,Q3,Q4,Q5) denote the one and only **physical relative 5D trajectory** per kick. The material clock is inherited unchanged from v18.52l. Distinguish:

1. **Line direction from four observed points:** m ≈ -0.308, angle beta=atan(m) ≈ -17.119°.
2. **Sheet/observer hinge tilt** around that line: alpha=±17.04°.

A line in the observed plane is L={(x,m x+b,0,0,0)}. Its physical unit axial direction and in-plane unit normal are

- `e=(1,m,0,0,0)/sqrt(1+m²)`;
- `n=(-m,1,0,0,0)/sqrt(1+m²)`.

Anchor `a=(0,b,0,0,0)`. Choose a unit vector `h=(0,0,h3,h4,h5)` in the hidden three-dimensional complement to the Q1-Q2 plane. Define

- `u=(q-a)·e`;
- `v=(q-a)·n`;
- `z=q·h`;
- `w_plus=cos(alpha) v+sin(alpha) z`;
- `w_minus=cos(alpha) v-sin(alpha) z`.

The two orthonormal observed planes are `(e, cos(alpha)n ± sin(alpha)h)`. **This changes the physical projection plane around the same candidate sloped hinge.** It does *not* rigidly rotate the existing 2D image. The construction preserves q(t) and material time. `w_plus+w_minus=2cos(alpha)v`, `w_plus-w_minus=2sin(alpha)z` are exact algebraic construction identities; they **must not** be misread as physical twin-lobe evidence.

**Line-fitting fork.** In the absence of the four observed 2D coordinates, `m=-0.308,b=0` is a fixed observational hypothesis (intercept unverified). If four coordinates are separately provided, fit constrained-through-origin OLS, affine OLS and total least squares, with leave-one-out stability. Only `--adopt-four-point-fit` changes the geometry to the fitted affine line.

**Provenance negative control.** The four crossings of historical `Q2=0` have `Q2≈0`; fitting these four points recovers a horizontal line, not the observer's `m≈-0.308`. They cannot substitute for the unprovided observational four points.

**Physical tests.** Observer gate windows are signed intervals satisfying `w_plus(t) w_minus(t)<0`, with exact durations for continuous piecewise-linear interpolation. Observer events are zeros of w±, whose interpolated 5D postures and material times are retained. Compare nearest cross-observer events, sequential same-observer postures, D, perpendicular radius `sqrt(v²+Q3²+Q4²+Q5²)`, axial progress u, projected covariance semiaxis ratio (apparent eccentricity proxy), and reduced angular phase atan2(z,v).

**Limitations.** Every choice of e,h,b,alpha is an observer hypothesis. Orthogonal in-plane rotations alone cannot change covariance semiaxis *eigenvalue ratios*, whereas true tilted projections can. A reduced 3D phase/winding cannot be certified as a topological invariant of an unconstrained open 5D trajectory. The gate pattern does not establish the user's causal shutter-door model. No exact Poincaré return, surface embedding, tube wall, cone, double helix or self-similarity is asserted. Full nonlinear integration is not rerun.
