# Q550 v18.51r — Scientific report

**Long-Horizon Paired-κ₀ Transverse Transport & Return-Section Search Audit**

## Local replay decision

The authenticated, **unchanged** nonlinear 10D Q550 solver was propagated to `xi=0.2` for two historical seed preparations, each at the two distinct historical κ₀ values and a midpoint parameter-curvature control, at three numerical solver settings: **18 authentic physical IVPs**. The run reached the prescribed endpoint on all branches; **17/17 built-in scientific consistency checks passed**. 513 saved sections per parameter branch extend the v18.51q observation horizon by a factor of four.

## Cross-κ₀ separation and phase

| Measurement | Original seed | Autonomous seed |
|---|---:|---:|
| Cartesian position difference at `xi=0.2` | 3.116924118e-06 | 3.116893446e-06 |
| Cartesian tau-velocity difference at `xi=0.2` | 1.797780806e-05 | 1.797753742e-05 |
| Local xi-equivalent lag at endpoint | -3.379001667e-07 | -3.377970403e-07 |
| Transverse fraction of position difference | 0.987907131 | 0.987914543 |
| Ratio to conservative numerical discrepancy proxy | 3.91684e+07 | 3.91681e+07 |
| Sections resolved above 20 times proxy | 512/513 | 512/513 |

The two κ₀ branches retain a small but numerically resolved position separation, about **15.35 times** the q-audit separation at `xi=0.05`, while around 98.8% of the difference remains perpendicular to the reference physical Cartesian velocity. This is a *local tangent projection in the selected Euclidean kinetic coordinates*, not a metric-independent torsion or global shadowing theorem.

A separately optimized **single constant phase delay**, held over the fixed window `[0.05,0.195]`, only decreases the RMS Cartesian position mismatch by approximately **0.85%**, for either seed. Consequently a constant shift along the same trajectory explains very little of the observed pair separation over this window.

## Clock and return search

The parent nonlinear clock identity is retained; physical Q endpoint at `xi=0.05` matches the authenticated q saved state within `7.06e-15` in Euclidean position norm.

**No sign change of `dx/dxi` nor nontrivial return to initial `x` was detected on `[0,0.2]`** for the six high-accuracy state branches. This is a fixed-grid projected search, not a certified exclusion between samples, nor a search over the full cycle. The reference physical x in the original seed decreases from approximately `0.504505889` to `0.504113124`. There is no measured full-state return in this interval, and no return-closure certificate is claimed.

## Authenticity, numerical limits and interpretation

Source ancestry is authenticated by rechecking the complete SHA-256 manifests for the v18.51q parent and the 19 historical ancestors re-established by the parent q executable. The same-seed parameter intervention and distinct historical seed preparations are preserved. Three DOP853 budgets provide **numerical discrepancy probes**, not certified trajectory-error enclosures.

The continuation offers a useful answer to the motivating 'shadow' analogy: the two distinct-κ₀ trajectories remain close in the kinetic Cartesian visualization but develop a resolvable, predominantly transverse offset that a local or global constant time shift barely removes. It remains unknown whether the offset saturates, rotates, changes its scaling, or closes at a later physical return section. An extended full return/reversal study with validated error control would be a separate experiment.

The `r` in the version is an audit version letter; it is not the radius `R` of the original cylindrical toroidal field geometry and is not the right-hand cocycle sheet `R=M-σG`.
