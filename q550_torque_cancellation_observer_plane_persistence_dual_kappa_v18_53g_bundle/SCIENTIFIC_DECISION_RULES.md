# Pre-registered scientific decision rules

1. Confirm exact SHA-256 of native NPZ, historical parent ZIP, and previous v18.53f authenticated reports. Fail closed.
2. At shared κ and zero tilt, integrated radial and phase torques must reproduce v18.53f independently.
3. Validate the material RHS against the original parent at selected states for **all three** κ values; reject if the field changes by implementation error.
4. Projected `cancellation_score` must be >= 0.90 and components have opposing signs for each tested plane to qualify as a strongly cancelling plane. Report **every** failure; do not take best plane or best angle.
5. Coordinate-free `bivector_cancel_score` >= 0.90 and `bivector_cosine` <= -0.90 are separate, stricter observer-independent conditions. They can fail even if a selected plane passes.
6. All seven angles and both historical κ prescriptions need to pass for the label `PERSISTENT_ACROSS_PREDECLARED_PROBES`. Failures must be reported by angle, κ, seed, and sample.
7. These κ variants are **same-state parameter counterfactuals only**. Even an all-pass label does not certify persistence on new dynamically propagated orbits.
8. The exact counterrotation charge is a distinct quantity; hinge torque is not a certified conserved charge.
9. If the authenticated parent or native states are missing, emit `PARTIAL_MISSING_AUTHENTIC_SOURCE_OR_NATIVE_STATES` and a sendback; never mark source-verified persistence.
