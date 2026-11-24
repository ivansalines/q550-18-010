# Decision rules for Q550 v18.52z

**Fail-closed authenticity:** SHA-256 of v18.51s full historical sendback, the v18.52l `extended_horizon_10d_states.npz`, and the embedded v18.52y server candidate report must match the manifest. Verify source period and first archived native replay before advancing Newton.

**Root not accepted:** If the one-turn chart closure norm is >=1e-8, or either direct 2-/3-turn closure threshold fails, or semigroup comparison, charge bound, independent-solver agreement, or exact leaf check fails, the result is a **non-certified near miss**. Do not report any Floquet multipliers or energy cycle identity as an observed closed cycle.

**Root accepted but DP unresolved:** If closure passes but three finite-difference step sizes or tighter-solver Jacobian do not converge, label `RELATIVE_FIXED_POINT_REACHED__FULL_FLOQUET_UNCERTIFIED`; raw Jacobian comparisons can be reported but NOT a certified Floquet spectrum.

**Conditional Floquet certification:** Only if root and all Jacobian consistency checks pass, record the eight eigensolutions, spectral radius, and conditional stability within a fixed-J leaf. Distinguish attracting, neutral/undetermined, or unstable by margin; do not call weakly neutral modes attracting.

**Energy:** If cycle is accepted, report E, drive work, x-work, radial dissipation, x-dissipation and both quadrature residuals. It is not proof of global conservation or historical full forward-return closure.

**Fallback:** If historical parent or native NPZ is missing, emit a report of missing inputs and (if strict flags were passed) nonzero exit plus valid server sendback. Never synthesize authentic Q550 outputs.
