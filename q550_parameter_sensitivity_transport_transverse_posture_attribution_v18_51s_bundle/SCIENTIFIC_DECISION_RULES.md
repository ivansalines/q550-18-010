# Scientific decision rules — v18.51s

1. Authenticate v18.51r's source and output manifests; inherited scientific controls must all pass. Do not rewrite the underlying authentic nonlinear generator.
2. In the cache-only baseline, use the exact same initial seed for both κ values and do not call the finite pair a variational derivative.
3. In the full model stage, source `J_Y` from the historic analytic `structural_jet`, not from differentiation of the saved trajectory alone. Build `partial_kappa f` at fixed state and fixed nonautonomous drive phase, checking two κ stencil widths.
4. Integrate the nonhomogeneous sensitivity ODE with zero initial derivative; require comparisons with the independently integrated high/low κ trajectories and with a finer jet sampling grid. Failure is a result, not a reason to silently loosen tolerances.
5. Record radial/azimuthal and tangential/transverse splits separately; a transverse Euclidean residue is not a material torque measurement.
6. Report limitations of ordinary float64, spline interpolation, three numerical solver budgets inherited from r and the proximity of the historical κ pair.
7. The long-horizon physical return, hysteretic memory, observer-independent material structure and globally valid κ law are NOT tested or certified by this audit.
