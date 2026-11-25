# Q550 v18.53z — local preflight and genuine lightweight replay

Validation was performed locally without running the 204-path full scan or any new nonlinear 11D/21D physical IVPs.

- Authenticated parent y server sendback CRC and SHA-256: PASS (15 sources, 8 outputs).
- Authenticated nested x/w/v/u/t ancestry: PASS (x 14/9, w 13/7, v 14/10, u 14/8, t 13/7).
- Historical v18.53f physical material RHS SHA-pinned in the distributed package.
- Synthetic phase recovery, 5D configuration lifts, half-stencil 11D Jacobian and exterior-square generator tests: PASS.
- Actual Q550 trace replay, original seed / shared kappa, preformal epoch, lambda = 0, 1, 4: PASS, 18 orientation comparisons.
- Parent y minimum-turn angular replay maximum difference in this genuine sample: 3.083e-12 degrees.
- Genuine Q550 gauge checks with two lifts and start/end assignments: 4 PASS.
- Genuine Q550 half-horizon source spin comparison: 1 PASS.
- ZIP packaging regression with deliberately created new directory and output manifest: PASS.

The audit does **not** certify 11D material holonomy or a unique complementary rotation. The full server scan is pending.
