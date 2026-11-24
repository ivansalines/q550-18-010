# v18.52g scientific decision rules

1. Source provenance passes only when the exact f/e/d archives, reference code and output hashes match, and the exact v18.51s source archive SHA is validated before model staging.
2. Historical continuity must independently match saved nonlinear trajectories, clocks and physical kick vectors on the original event-clock interval 0..0.2; failure means no interpretation of the new horizon.
3. All eight original/autonomous, ±, 1e-4/3e-4 finite kicks retain the original fixed-κ field and identical seed-specific initial data.
4. Material clock must stay positive along each branch; all scalar diagnostics use independently inverted matched τ and record phase mismatch.
5. Inward acceleration I>0 alone is never described as proven restoration. C<=0 alone is never reported as a decrease in position distance.
6. An observed local maximum in D requires G to cross positive→negative in a valid numerical bracket. Tangencies or isolated grid zeros are reported as unresolved.
7. Event roots must be re-bracketed on a finer independently integrated model when `--refine` is used. Even this is a numerical convergence test, not a continuous-time certificate.
8. If no event is detected, conclude only that **no event was detected** in the scanned finite numerical horizon, not that no event exists between samples or at later time.
9. Keep all failed checks and incomplete branches visible; do not turn a solver failure into a no-event claim.
10. Server results, generated figures, source scripts, input ZIPs, SHA-256 manifests and full descriptive file/folder names must remain together in the output bundle.
