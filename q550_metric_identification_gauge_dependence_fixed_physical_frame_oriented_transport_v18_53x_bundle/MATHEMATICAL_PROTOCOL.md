# Q550 v18.53x — Mathematical protocol

## Parent geometry
The parent v18.53w has 204 combinations `(seed, kappa_rule, epoch, lambda)`. The parent `v18.53v` provides 141 optimized widths along a fixed-center path for each node; v18.53t supplies the genuine 21D-integrated cumulative bivector `C(xi)`.

`N(a,w) = C(a+3w) - C(a)` is an element of the exterior square `Lambda^2 R^5`, stored in a fixed mass-normalized Cartesian physical 5D frame and represented by 10 components. For an SPD sensitivity metric `H_lambda = Lambda^2 G_lambda`, the historical orthonormal metric-coordinate representative is `B_lambda = H_lambda^(1/2)`, the symmetric positive-definite square root. For a sampled path, the normalized vectors `u_i=B_lambda N_i / ||B_lambda N_i||` yield a 45D oriented geodesic polygon `P`, via the adjacent great-circle arcs `theta_i / sin(theta_i) (u_i wedge u_(i+1))`.

## Metric identification gauge
If `O_lambda` belongs to SO(5) in the mass-normalized Cartesian frame, then `R_lambda=Lambda^2 O_lambda` belongs to SO(10), and `B'_lambda=R_lambda B_lambda` obeys `(B'_lambda)^T B'_lambda=H_lambda` exactly. Thus the metric, all norms, and all same-lambda pairwise intrinsic angles remain unchanged. For each polygon `P`, the corresponding 45D representative changes as `P'=(Lambda^2 R_lambda) P`.

The gauge frames are predefined, use no data fitting, are identity at lambda=0, and are constant along each *individual* xi path:
- `symmetric`: identity.
- `intraplane`: 5D plane (0,1) rotation by `0.10 lambda` radians.
- `mixed`: plane (0,4) rotation by `0.08 lambda` then plane (1,3) rotation by `-0.06 lambda`.
- `countermixed`: opposite signs of the mixed rotations.

Cross-lambda angles to the canonical lambda=0 polygon can change under these pure gauge choices. The same-lambda fixed-center versus moving-center angle must remain unchanged when both polygons are rotated by the same gauge. De-gauging by `(Lambda^2 R_lambda)^T` must reproduce the historical vector to near machine precision.

## Fixed physical frame, without metric identification
For every parent node, reconstruct the physical 10D sampled path from the authenticated v widths and the 21D cubic-Hermite cumulative primitive. Evaluate the 45D polygon using the **same fixed canonical Euclidean metric** `I_10`, regardless of `lambda`. Compare against the historical canonical baseline `P000`. This isolates geometry caused by a change of selected path in a declared fixed Cartesian physical frame. It is not a coordinate-free material holonomy.

For lambda in {0,1,2,3,4}, additionally reconstruct the whole *moving-center* path with the parent width-selection algorithm. Compare its canonical fixed physical frame polygon to `P000` and verify that the parent weighted polygon `P111` is reproduced.

## A deliberately moving gauge along xi
At lambda=4 on the moving-center paths, use a known xi-dependent SO(5) frame whose endpoint values equal identity. Show that a naive chord-angle sum on the moving-gauge vectors changes, even though the underlying physical path does not. Then apply the **explicit pure-gauge connection** `R_i R_(i+1)^T` to the next vector, followed by pullback to the common frame, and recover the original arc length and 45D polygon.

This is a numerical covariance test of an **imposed** moving gauge, not a derivation of the model's unique material connection. No physical 11D return, invariant torsional flip, or material holonomy is certified.

## Gates
The parent w/v/u/t manifests and ZIP CRC must pass. The 204 parent fixed-width weighted polygons must replay with error <= 1e-7, the 60 reconstructed parent moving-center polygons with error <= 2e-7. Gauge de-identification and corrected moving-frame polygon defects must each be <= 5e-11. No new IVP is run in this audit. The entire studied near-return window uses formal clock continuation beyond xi=25.6.
