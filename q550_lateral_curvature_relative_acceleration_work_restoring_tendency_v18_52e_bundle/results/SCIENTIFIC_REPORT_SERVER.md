# Q550 v18.52e — Server numerical results

Numerical checks passed: **13/13**.

Eight genuine 10D nonlinear kicks were integrated at fixed kappa and independently compared with the archived v18.52d norm traces.

| Seed and kick | Terminal relative speed ratio | P_rel < 0 fraction | d(|delta Q|²)/dtau > 0 fraction | Qualified inward acceleration fraction |
|---|---:|---:|---:|---:|
| original_negative_1.0e-04 | 0.916481713 | 1.000000 | 1.000000 | 1.000000 |
| original_positive_1.0e-04 | 0.916462573 | 1.000000 | 1.000000 | 1.000000 |
| original_negative_3.0e-04 | 0.916500849 | 1.000000 | 1.000000 | 1.000000 |
| original_positive_3.0e-04 | 0.916443429 | 1.000000 | 1.000000 | 1.000000 |
| autonomous_negative_1.0e-04 | 0.916495456 | 1.000000 | 1.000000 | 1.000000 |
| autonomous_positive_1.0e-04 | 0.916476306 | 1.000000 | 1.000000 | 1.000000 |
| autonomous_negative_3.0e-04 | 0.916514600 | 1.000000 | 1.000000 | 1.000000 |
| autonomous_positive_3.0e-04 | 0.916457152 | 1.000000 | 1.000000 | 1.000000 |

**Do not infer recovery solely from negative relative work.**
The power is a Euclidean relative-kinematics proxy, not physical dissipated energy; the full position separation and clock must also be examined.

Scientific classification and numerical tolerances are in `report.json`.
