# Q550 v18.51t — Scientific decision rules

1. The v18.51s server sendback must pass archival source/output SHA-256 and CRC checks; its 16/16 full-mode tests must hold. If not, **do not** run or promote a later inference.
2. Source split is `b = b_v2 + b_omega2 + b_other` in the explicit historical 10D state chart. Unknown `other` contributions are retained as a measured third channel, not forced to zero.
3. The three independent forced IVPs must add to the unsplit IVP at every section, and their Green-integral counterparts must add as well.
4. Green's representation must agree with direct integration on Simpson-closed even nodes. Check conditioning of `Phi`; a numerically unstable Green integral is a failed test, not a physical discovery.
5. Source derivative step refinement, sampled `J,b` grid refinement and independent v18.51s endpoint and v18.51r physical κ-pair comparisons must all pass.
6. Declare signed radial/azimuthal source transfers and cancellation, including cross-channel entries. Do **not** infer response dominance from the instantaneous norm of a source alone.
7. Passing the audit establishes on-trajectory first-order parameter-source response in the fixed historical coordinate chart over `0≤xi≤0.2`. It does not establish parameter-independent force modularity, material torque, topological winding, full 10D return, global shadowing, or interval-validated error bounds.
