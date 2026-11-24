# Q550 v18.51v — scientific decision rules

**PASS** means all numerical diagnostics meet their declared tolerances using the saved authenticated physical trajectories and Green kernels. It does not mean a material memory law or a complete physical return has been proven.

**Hard integrity constraints:** U/T/S/R nested ancestry and both cache copies must match the source SHA256 lineage; all reports in the ancestry must have passed their scientific checks. The parameter midpoint, two distinct historical kappas and positive tau clocks must be present.

**Physical-clock reparameterization:** expect `max|d tau/dxi-c|<5e-9`, `max|K_tau*c-K_xi|<1e-13`, causal upper triangle `<1e-15`, and independent nonuniform tau quadrature error `<1e-8` in absolute source-response units.

**Fixed-tau sensitivity:** expect midpoint model sensitivity versus nonlinear fixed-xi secant `<2e-7` absolute, nonlinear same-tau interpolation versus chain-rule clock correction `<2e-7`, and clock-corrected midpoint model sensitivity versus direct nonlinear same-tau comparison `<2e-7`. The effective clock correction is measured, not assumed zero.

**Moving-frame consistency:** for both independent passive observer schedules expect reconstructed 2D physical vectors, physical Green kernels and vector norms to agree within `1e-12`, source channel additivity `<1e-9`, and connection-aware derivative reconstruction `<2e-6` under actual saved trajectories.

**Interpretation checks:** if the radial/azimuthal cancellation fraction changes under a frame rotation, record the change as a projection effect, not as a physical alteration of the model. Do not assign the entire clock-shift sensitivity to an individual source without a separate decomposition of the parameter-dependent clock.

Failure of a numerical threshold must be reported as a diagnostic failure, not silently relaxed. If a metric or nonuniform quadrature test fails, examine time-discretization and floating-point conditioning before making physical claims. If the two seed families behave differently, retain their distinct results.
