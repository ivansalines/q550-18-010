# Scientific decision rules

1. **Numerical integrity** passes only when source hashes and CRC match, identities of independently rebuilt five-dimensional vector scalars match archived values, spline/quadrature closure is within recorded tolerances, and all vectors are finite.
2. **Sampled inward insufficiency** means every sampled `C=P-I` is strictly positive. Report exact minimum `C/P` and maximum `R`; do not call this a certified continuous ODE bound.
3. **Sampled convexity threshold** occurs when `R>=1` (equivalently `C<=0`) at a sample. Record location and sign; do not infer position return.
4. **Sampled halt of outward position separation** occurs when `G<=0` at a **nontrivial** sample after the initial coincidence. Record separately. A sign reversal in `G` is stronger than a convexity threshold.
5. **Observed position return** requires actual reduction in `||q||` after a preceding increase, not merely inward acceleration.
6. Report all cases, including scientific failure, independently of test/pass technical status.
7. `F_needed` is a frozen-state scalar counterfactual only, **never** a valid off-orbit force or recovery certificate.
8. The physical effect of extending observation beyond the prior final time remains **unmeasured**; no extrapolated threshold or asymptotic stability claim is allowed.
