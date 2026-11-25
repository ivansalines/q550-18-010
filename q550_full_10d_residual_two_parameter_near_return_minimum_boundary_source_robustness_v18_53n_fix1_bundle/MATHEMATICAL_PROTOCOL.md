# v18.53n — Mathematical protocol

Let `C(ξ)` be the cumulative integrated 5×5 skew-symmetric residual, mapped isometrically to 10 coordinates via `sqrt(2) B_ij` for `i<j`. The sum of any three consecutive windows of common width `w` is

`N(o,w)=C(o+3w)-C(o)`.

With `s=C'(o)` and `e=C'(o+3w)`, the full unprojected residual Jacobian is

`J=[ e-s , 3 e ]`, an exact endpoint-source identity wherever the interpolant is differentiable.

For the scalar objective `f=1/2 N^T N`,

`gradient f = J^T N`.

For an interior minimum: `J^T N=0`. At parameter-box boundaries the KKT projected-gradient conditions replace zero unconstrained gradient. This condition does **not** require `N=0`, and may be satisfied while the physical history has not returned.

Inside open cells of the piecewise-linear cumulative interpolant, `C''=0`, so the local cell Hessian of `f` is `J^T J`; this is NOT a claim about the continuously evolved physical model's actual second variation. For numerical robustness, compare PCHIP and decimated quadrature and physical tight IVPs.

The area element and principal singular values of `J` diagnose parameter sensitivity; a high condition number identifies a poorly located minimum. Distinct observer-plane roots from v18.53m are re-evaluated on the same full 10D norm for comparison. Only full 10D minima are optimized.

Trace domain: `o in [12.8,15.4]`, `w in [12.70,12.91]`, `o+3w<=56.0`. No extrapolation. Use physical trace input authentication before numerical claims.
