# v18.53o — Mathematical protocol (English)

Let the isometric 10-coordinate representation of the accumulated 5x5 antisymmetric residual be `C(xi)`. Compute `N=C(e)-C(o)`, with `e=o+3w`. This telescopes three adjacent windows exactly.

- Fixed-width tangent: `D_o N = g(e)-g(o)`, `g=C'`.
- Width tangent: `D_w N = 3g(e)`.
- Two-ended stationarity: `N·g(e)=0` and `N·g(o)=0` follows algebraically from `J^T N=0` for an interior minimum, regardless of how small ||N|| is.
- For `u=N/r`, scalar growth `r' = u·N'` and directional change `u'=(I-uu^T)N'/r`. The antisymmetric 10D orientation form `u ∧ u'` has norm `||u'||` under the standard wedge coordinate norm.
- The corresponding boundary decompositions are `g(e)` and `-g(o)`; their scalar and transverse contributions must be recorded separately, not called energy flux.
- For a locally optimized path `w(o)`, `N'=g(e)(1+3w')-g(o)`, with `w'` estimated from adjacent optimized widths. Check chain-rule rate against point differences away from knots.
- Cosines of oriented differential forms compare the rotation direction across the minimum without assigning a sign to SVD eigenvectors.

A near-constant scalar norm with a changing normalized direction is an observable of this defined 10D curve; it is not a proof of physical holonomy or underlying two-oscillator phase locking.
