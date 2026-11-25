# v18.53q — Decision Rules

- A run only starts if p and nested o CRC, SHA-256 source/output manifests, expected scientific decisions and 10 parent IVPs are verified. This version performs **zero new physical IVPs**.
- Source match at a clock knot is an interpolation constraint, not a physical time-convergence certificate. Report C/g inconsistency rather than enforcing agreement.
- Require all reconstructed residual norms to be positive. Report global width candidates and boundary minima; do not silently relabel a branch jump as continuous dynamics.
- An interior positive-curvature stationary point may support implicit w'(o). Denominator-zero/negative-curvature cases are excluded rather than assigned an artificial tangent.
- Geodesic length and integrated oriented form are separately compared. Their agreement cannot certify physical holonomy.
- If the three reconstructions give different path lengths or derivative patterns, classify this as interpolation sensitivity until a new physical source density test is completed.
- If they agree, classify only numerical robustness *under these interpolants*, not the model's continuum limit.
- Formal clock continuation after xi=25.6 remains explicitly marked. No full 11D material recurrence or Q550/coupled-oscillator identity is claimed.
