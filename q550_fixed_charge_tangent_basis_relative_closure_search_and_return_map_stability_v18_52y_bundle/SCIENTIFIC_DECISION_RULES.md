# Scientific decision rules — Q550 v18.52y

| Gate | PASS requires | If it fails |
| --- | --- | --- |
| Authenticity | SHA match for v18.51s and v18.52l, known source parameters, fixed kappa | `AUTHENTIC_AUDIT_FAILURE`; no Q550 certification |
| Source stroboscopy | late sign exactly -1, external phase changes by -2π in one source period and -4π in two | reject periodic-map interpretation |
| Exact fixed-J geometry | rank 8 and analytic/FD tangent match; ∇J·B negligible | reject fixed-J derivative |
| First return replay | DOP853 source continuation reproduces independent native first return | reject further interpretation |
| Gauge quotient | counterrotated initial gauge gives same reduced one-turn state | reject quotient |
| Return repeat | direct P² matches composition of P twice | reject autonomous late stroboscopic map |
| Full derivative | 8 independent central probes on the same J leaf, two FD steps | report only where well conditioned; no Floquet claim |
| Relative fixed point | independently recomputed first and two-turn residuals below thresholds | no cycle certificate, even if phase nearly repeats |
| Conditional Floquet | true fixed point plus two-step Jacobian convergence | do not report cycle multipliers |
| Conditional attracting cycle | all conditional multipliers strictly inside unit disk with margin | no attraction claim |
| Conservative dynamics | independent global work–dissipation evidence, not only charge J | never inferred by this audit |

**Result categories:** `AUTHENTIC_ANCESTOR_MISSING__SYNTHETIC_ONLY`, `AUTHENTIC_AUDIT_FAILURE`, `AUTHENTIC_FIXED_CHARGE_8D_TANGENT_MAP_AND_TWO_TURN_PROPAGATION__NO_RELATIVE_CYCLE_CERTIFICATE`, `RELATIVE_FIXED_POINT_CANDIDATE__FLOQUET_UNRESOLVED`, and `RELATIVE_FIXED_POINT_AND_LOCAL_FIXED_J_RETURN_SPECTRUM_CONDITIONALLY_CERTIFIED`.

The user runs IVPs on their server; local ChatGPT work comprises only source inspection, symbolic/linear tests, smoke packaging, and file manifests.
