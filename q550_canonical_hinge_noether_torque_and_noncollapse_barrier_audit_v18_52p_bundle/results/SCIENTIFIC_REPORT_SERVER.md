# Q550 v18.52p — Scientific server report

**Stage:** `FULL_RESOLUTION_NATIVE_10D_AND_RELATIVE_5D_KINEMATIC_AUDIT`

## Exact canonical coordinate identity
Fixed observation-based nodal slope `m=-0.308`; passive nodal angle `theta=-17.118832336 deg` (negative). Visual hinge tilt is a DIFFERENT angle, `±17.04 deg`.
Symplectic block matrix error `1.873e-17`; max coordinate-transform J error `2.118e-22`. These are algebraic identities, NOT tests of dynamical conservation.

## Finite-horizon physical relative kinetic observations
| Physical branch | J min | J max | sign flips | relative drift range/median | min r2 | min r5 | open path turns |
|---|---:|---:|---:|---:|---:|---:|---:|
| original_negative_1.0e-04 | -5.68755e-08 | -1.45668e-08 | 0 | 1.31 | 0.000177762 | 0.000177869 | -1.67488 |
| original_negative_3.0e-04 | -5.11866e-07 | -1.31091e-07 | 0 | 1.31 | 0.000533369 | 0.000533698 | -1.67487 |
| original_positive_1.0e-04 | -5.6877e-08 | -1.4568e-08 | 0 | 1.31 | 0.000177735 | 0.000177838 | -1.67489 |
| original_positive_3.0e-04 | -5.11906e-07 | -1.31123e-07 | 0 | 1.31 | 0.000533122 | 0.000533425 | -1.6749 |
| autonomous_negative_1.0e-04 | -5.68557e-08 | -1.45628e-08 | 0 | 1.31 | 0.000177707 | 0.000177811 | -1.67487 |
| autonomous_negative_3.0e-04 | -5.11687e-07 | -1.31055e-07 | 0 | 1.31 | 0.000533204 | 0.000533523 | -1.67486 |
| autonomous_positive_1.0e-04 | -5.68572e-08 | -1.4564e-08 | 0 | 1.31 | 0.00017768 | 0.00017778 | -1.67488 |
| autonomous_positive_3.0e-04 | -5.11728e-07 | -1.31086e-07 | 0 | 1.31 | 0.000532957 | 0.000533251 | -1.67489 |

## Absolute canonical kinetic proxy versus relative two-flow quantity
Native original 10D states available: `True`.
From exact historical kinetic embedding, the first polar pair gives `J1_kinetic=r1^2*w1`, and second pair `J2_kinetic=r2^2*w2`. Both are reconstructed and numerically differentiated along EACH physical material clock. The original action / true canonical momentum in presence of forcing is not supplied, so the proxy is NOT a certified Noether charge.
| Native physical trajectory | J1 min | J1 max | J1 relative range | J2 sign flips |
|---|---:|---:|---:|---:|
| original__reference__Y_tau | -1.44657 | -1.4365 | 0.0069863 | 7 |
| autonomous__reference__Y_tau | -1.44658 | -1.4365 | 0.0069887 | 7 |
| original_negative_1.0e-04__kick__Y_tau | -1.44654 | -1.43647 | 0.0069863 | 7 |
| original_negative_3.0e-04__kick__Y_tau | -1.44647 | -1.4364 | 0.0069863 | 7 |
| original_positive_1.0e-04__kick__Y_tau | -1.44661 | -1.43654 | 0.0069863 | 7 |
| original_positive_3.0e-04__kick__Y_tau | -1.44668 | -1.43661 | 0.0069862 | 7 |
| autonomous_negative_1.0e-04__kick__Y_tau | -1.44654 | -1.43647 | 0.0069887 | 7 |
| autonomous_negative_3.0e-04__kick__Y_tau | -1.44647 | -1.4364 | 0.0069887 | 7 |
| autonomous_positive_1.0e-04__kick__Y_tau | -1.44661 | -1.43654 | 0.0069886 | 7 |
| autonomous_positive_3.0e-04__kick__Y_tau | -1.44668 | -1.4366 | 0.0069886 | 7 |

The quantity `J_kick-J_reference` is not the same as `(Q_kick-Q_reference) x (V_kick-V_reference)`; an explicit matched-time comparison is written to CSV.

## Noncollapse: necessary conditions rather than a false theorem
In the fixed physical 2D plane, `T2 = 1/2 (dr/dtau)^2 + Jkin^2/(2r^2)` and `rddot = Jkin^2/r^3 + (q_plane·a_plane)/r` are exact kinematic identities. Observed `J_min>0` and finite sampled speed yield a posteriori `r >= min|J|/max|v2|` on the sampled domain only; no global speed bound, exact Noether conservation or continuous-time zero-exclusion proof is present.
Open-path unwrapped angular turns are not integral topological charges.

## Scientific disposition
**Noether charge:** not certified. **Topological charge:** not certified. **Global impossibility of collapse:** not certified. A physically meaningful finite-horizon geometric obstruction may still be present and is measured here.
