# Q550 v18.51y — Scientific Decision Rules

1. **PASS source authenticity:** the T, W and X snapshots and parental report must agree with their SHA-256 manifests. When present, the full X ancestry must also pass nested ZIP CRC and all X source/output hash checks.
2. **PASS strict support identity:** at every historic source sample, `P5 b + P7 b + Pother b = b`, with nonzero model entries at slots 5 and 7. Any tested nonzero symmetric X mixing must leak into the other labelled support; equality of total forcing alone is not admissibility.
3. **PASS dynamic coupling observation:** the J columns at 5 and 7 must have finite nonzero off-support components. Their Euclidean norms are diagnostics, not invariant physical energies.
4. **PASS coordinate covariance:** with `T(xi)` and its derivative, transported projectors select the same physical source, `T^-1(T S)=S`, and the moving tangent generator satisfies the exact connection identity within floating-point accuracy.
5. **PASS negative controls:** frozen index projectors under a rotating chart must develop false leakage, and dropping `T' T^-1` must produce a detectable nonzero defect.
6. **NO causal overclaim:** no independently variable forcing channels, action-level source decomposition, nontrivial physical metric invariance, return-event, or constitutive memory certificate is asserted.

Expected numerical tolerances: source and tangent covariance max errors <=2e-11; negative-control frozen index leakage >0.25; both seed families checked; 257 source/J samples and 513 stored response samples per seed.
