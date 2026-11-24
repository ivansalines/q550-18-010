# Q550 v18.51z — Mathematical Protocol

## Scientific question

How can the summed response vector turn opposite to both separately anchored source-channel responses, and which apparent turning reversals can a moving passive observer create without changing the Cartesian result?

The playful title ("Posture, Elegance & Observer-Induced Catwalk") labels the numerical **posture of parameter-sensitivity response vectors**, not the physical posture of a human or a claim of mechanical torque.

## Source genealogy

The exact v18.51y server sendback is bundled and SHA256-locked. It includes authenticated v18.51t Cartesian/source snapshots, v18.51w channel responses and angles, and v18.51x and v18.51y results and scripts. The v18.51y report is checked for `all_passed=true`. This audit does **not** rerun the nonlinear model or Green IVPs.

## Model-derived response vectors

For each seed, take the 2D Cartesian `(Q3,Q4)` responses `A` (from the `v2` component of the model-derived `partial_kappa f`), `B` (from `omega2`), `O` (remaining sources), and `C=A+B+O` (total). The original channel decomposition remains defined only under the historical fixed-support labeling conditions established in v18.51y. There is no independent-actuator claim.

Use only `xi > 0` for angles; the quantitative turning-rate and sign tests use `xi >= 0.01` so near-zero vectors do not create arbitrary angular instabilities.

## Signed turning and source-channel interference

For 2D vectors `u,w`, define `u wedge w = u_x*w_y-u_y*w_x`. The signed angular rate of nonzero `C` is

`Omega_C = (C wedge C') / ||C||^2`.

The numerator decomposes exactly:

`C wedge C' = sum_j (V_j wedge V'_j) + sum_{j<k} (V_j wedge V'_k + V_k wedge V'_j)`

with `V_1=A`, `V_2=B`, and `V_3=O`. The mixed terms contain both amplitude changes and orientation changes. Dividing each by `||C||^2` provides a signed turning *budget*, not an energy or force decomposition. The numerical closure is checked on the authentic 513-sample arrays.

## Passive observer and gauge connection

Apply the **same** proper rotation to all vectors:

`Vtilde = R(-phi(xi)) V`, where

`phi(xi) = (8 deg)*(xi/0.2) + (2.5 deg)*sin(4*pi*xi/0.2)`.

Then `alpha_tilde_j = alpha_j-phi` for each channel and

`d(alpha_tilde_j)/dxi + d(phi)/dxi = d(alpha_j)/dxi`.

The relative angle, dot product, norm, and oriented wedge stay unchanged under this common SO(2) rotation. Absolute component rates and apparent turning direction do not. Deliberately omitting `phi'` is a negative control.

## Physical-clock reparametrization

Use the `tau(xi)` recorded by v18.51w and verify `dtau/dxi > 0`. Then

`d(alpha)/dtau = [d(alpha)/dxi]/[dtau/dxi]`.

A positive clock rate can change angular-rate magnitudes but not their signs. Compare the chain-rule result with direct numerical differentiation in `tau`, allowing finite-difference truncation errors.

## Controls and limits

1. Verify parent ZIP integrity, SHA256, and `all_passed` status.
2. Reconcile W Cartesian channel vectors with independent T snapshot arrays.
3. Check three-channel reconstruction and signed turning numerator closure.
4. Check the local clockwise/clockwise/counterclockwise behavior in both seeds.
5. Verify mutual-angle and area invariance under common passive rotation.
6. Verify the connection term and demonstrate the false turning created if it is omitted.
7. Verify positive `tau` and the derivative chain rule.
8. Use separate synthetic test vectors with known motion, a reflection negative control, and an origin-angle mask.

This is a **local, model-derived response geometry diagnostic** on `xi in [0,0.2]`: it neither certifies a 10D material return nor a constitutive memory law. The pointwise wedge-budget closure is an algebraic identity measured on authenticated numerical trajectories, not an independent nonlinear perturbation.
