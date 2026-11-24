# Q550 v18.51w — Scientific Decision Rules

The audit must not claim a new dynamical integration. Acceptance requires:

- `SOURCE_MANIFEST.json` full verification, followed by source and output hashes within the authenticated v18.51v server sendback.
- The V physical second-plane lab-frame channel vectors agree with independent T Cartesian vectors at a numerical tolerance of `1e-9`.
- The SO(2) signed-area and dot-product invariance defects must be below `1e-12`; the reflection negative control must reverse the area.
- Every available critical observer axis must produce a projected zero within absolute tolerance `1e-12`, and transform modulo `pi` under the prescribed observer rotations.
- The normalized dot/area identity and source-channel additivity must pass within `1e-12` or better.
- The initial zero-response section must be masked as *undefined*, not treated as 90 degrees or any preferred direction.
- Both seed families must be included and their terminal mutual angles must agree within `0.001 degree`.
- Any inferred temporal trend or sign persistence is restricted to the **sampled** 513 sections of `xi in [0,0.20]`. No theorem of global monotonicity is inferred.

Interpretation must distinguish (i) mutual channel angle, (ii) magnitude-dependent oriented area, (iii) relative projection-zero axes, and (iv) observed cancellation percentage. Invariance of the first two under proper frame changes cannot be extended to observer-independence of channel selection, the signed area under reflections, or material-memory certification.
