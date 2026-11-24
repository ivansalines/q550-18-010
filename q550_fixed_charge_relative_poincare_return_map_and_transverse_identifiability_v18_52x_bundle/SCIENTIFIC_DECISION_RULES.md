# Scientific decision rules — v18.52x

- **FAIL CLOSED** for absent SHA-pinned v18.51s ancestry or v18.52l native states when `--require-parent --require-native` is used; do not silently substitute a guessed phase or fabricate an initial 10D state.
- Original phase must be reconstructed through the authenticated historical `ctx['theta']`, never the eleventh native column.
- Late source clock must be strictly negative, with `cubic_clock_sign == -1` at representative points; otherwise a single-period stroboscopic interpretation is rejected.
- Preserve only the true counterrotation group action; physical angular velocities must not be gauge-removed arbitrarily.
- **Every individual orbit's exact J** must be conserved numerically (max observed drift <1e-9); verify independently before comparing returned states.
- **Charge offsets across kicks may be nonzero.** Report them; do not force/renormalize actual kicked trajectories onto the reference leaf.
- At least eight independent same-charge directions are required to infer the full 8D transverse Jacobian without an analytic tangent solver. Four kick points, even if independent, are insufficient. Refined duplicates do not raise rank.
- One phase winding in each stored history measures a single application of the return map, not a repeated discrete orbit, not a limit cycle and not a Floquet spectrum.
- A relative periodic orbit and attraction must not be declared unless appropriate exact/repeated closure and independently stable transverse spectrum are established.
- The raw coordinate norm combines different physical units. It is a section mismatch diagnostic only, never a canonically normalized physical distance.
- Verify user-provided large inputs against pinned SHA-256 before using; log failures in `results/FAILURE_DIAGNOSTIC.txt` and package the report.
