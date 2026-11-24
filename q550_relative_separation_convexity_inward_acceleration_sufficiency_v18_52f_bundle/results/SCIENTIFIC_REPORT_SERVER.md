# Q550 v18.52f — Server replay scientific report

Numerical checks: **48/48**.
Authentic physical trajectories: **8**, sampled sections **4104**.
Maximum sampled inward-to-positive ratio: **0.141399428** (threshold=1).
Final required inward factor range (frozen-state diagnostic only): **7.072164–7.081873**.
Minimum normalized convexity across all sampled states: **0.858600572**.
Scientific classification: **NO_SAMPLED_CONVEXITY_REVERSAL_OR_POSITION_RETURN**.

| Kick | R at final | Minimum C/P | Required factor at final | First sampled C<=0 | First sampled G<=0 |
|---|---:|---:|---:|---|---|
| autonomous_negative_1.0e-04 | 0.1412572 | 0.8587428 | 7.07929 | None | None |
| autonomous_negative_3.0e-04 | 0.1412056 | 0.8587944 | 7.08187 | None | None |
| autonomous_positive_1.0e-04 | 0.1413088 | 0.8586912 | 7.07670 | None | None |
| autonomous_positive_3.0e-04 | 0.1413606 | 0.8586394 | 7.07411 | None | None |
| original_negative_1.0e-04 | 0.1412961 | 0.8587039 | 7.07734 | None | None |
| original_negative_3.0e-04 | 0.1412445 | 0.8587555 | 7.07992 | None | None |
| original_positive_1.0e-04 | 0.1413477 | 0.8586523 | 7.07475 | None | None |
| original_positive_3.0e-04 | 0.1413994 | 0.8586006 | 7.07216 | None | None |

A loss of separation convexity (C<=0) would **not by itself** imply separation return (G<0).
No rigorous continuous-time sign enclosure, global dynamics, or asymptotic recovery certificate is claimed.
The parent nonlinear trajectories were not reintegrated in this audit.
