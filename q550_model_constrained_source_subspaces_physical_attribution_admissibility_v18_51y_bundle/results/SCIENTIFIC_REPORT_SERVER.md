# Q550 v18.51y | Server Scientific Report

Status: **PASS** (18/18)

**Question:** Which source assignments are physically anchored in the historical coordinate chart, and what changes under a consistent coordinate gauge?

The model defines one partial-kappa source vector; its v2 and omega2 coordinate supports (zero-based indices 5, 7) are fixed. Under strict disjoint support the decomposition is unique. The X mixing changes supports, not underlying source physics.

A time-varying chart T(xi) must transform P, J, b, and response together and retain Tdot*T^-1. Omitting projector transport or connection creates chart artifacts.

| Seed | mean off-column J v2 | mean off-column J omega2 | X critical eta | support leak 1 at critical | support leak 2 at critical | max covariant generator defect |
---|---:|---:|---:|---:|---:|---:|
| original | 0.93481159 | 0.88711042 | 0.060231477 | 0.038631655 | 0.10565738 | 3.55e-15 |
| autonomous | 0.93485928 | 0.88717169 | 0.060251848 | 0.038654909 | 0.10566961 | 3.55e-15 |

**Limitation:** No independently adjustable source knobs, unconstrained physical-metric invariance, or material-memory law has been established.
