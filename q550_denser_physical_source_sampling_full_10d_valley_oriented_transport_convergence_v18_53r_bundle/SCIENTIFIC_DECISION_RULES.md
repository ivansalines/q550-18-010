# v18.53r — Scientific decision rules

1. Reject unauthenticated source archives, mismatched exact parent SHA hashes, unverified original external Q550 numerical inputs, missing physical modules or incomplete solver runs. Synthetic tests cannot replace physical integrations.
2. Record 12 new IVPs by default (six 2× and six 4×) with exact identical historical state, parameters and start/end ξ. The 1× parent is an inherited authenticated result. Verify tightened DOP853 maximum step, finite states, finite physical source and bounded charge drift by the original parent tolerance.
3. Require replays at original knots to agree with authenticated original state/primitive/source within measured tolerances; report maxima. Do not claim independent physical realizations simply because a new t_eval is used.
4. For the width valley, report cases where an optimizer hits bounds, stationarity root bracketing fails, singular implicit derivatives or different minimizer branches occur. For each case/level/interpolant, retain actual coordinates, vector directions and oriented 45D integrals.
5. Evaluate numerical convergence across 1×/2×/4×. An apparently stable minimum and orientation is only **numerical transport robustness**. `NO_MATERIAL_HOLONOMY` remains an explicit scope limit. If computed residual <1e-5, flag it for independent follow-up instead of asserting material closure.
6. All tail points ξ>25.6 are a formal model-clock continuation. Preserve full provenance and distinguish genuine 11D source reintegration from post-processing-only work.
