# Scientific decision rules — v18.52c

A numerical PASS means successful data authentication, algebraic reconstruction, finite-difference diagnostics, and any requested physical midpoint probes. It does **not** certify continuous physical positivity.

Possible scientific outcomes:

- `PHYSICAL_MIDPOINTS_POSITIVE__CONDITIONAL_CONTINUUM_BOUND_ONLY`: All 24 independent nonlinear/tangent midpoint evaluations retain positive signed margins. Smoothness between samples remains an assumption.
- `SAMPLED_PHYSICAL_SIGN_LOSS_FOUND`: At least one newly integrated branch has a nonpositive physical margin.
- `PHYSICAL_MIDPOINT_PROBES_PARTIAL_OR_UNRESOLVED`: Solver, model reconstruction, or tangent propagation limits prevented full physical assessment.
- `PHYSICAL_MIDPOINT_COVERAGE_INCOMPLETE__CONDITIONAL_BOUND_ONLY`: The requested probe grid was intentionally restricted or omitted.

Separately, report the exact PCHIP extrema as **surrogate-only**, the observed parameter/time secant envelopes, and the conditional positivity margin under an explicitly stated (unproven) derivative majorant. A positive conditional margin cannot be described as a theorem about the nonlinear ODE.

A fully rigorous physical continuous certificate would additionally require validated (e.g. interval) enclosures of the forward nonlinear flow, its first kappa sensitivity and relevant time derivatives, regular denominators in the angular budget, and rigorous global numerical errors across the entire kappa–xi rectangle.
