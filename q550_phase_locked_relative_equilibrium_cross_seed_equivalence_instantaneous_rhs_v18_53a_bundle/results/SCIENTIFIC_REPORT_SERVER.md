# Q550 v18.53a — phase-locked relative equilibrium audit

**AUTHENTIC_INSTANTANEOUS_PHASE_LOCKED_RELATIVE_EQUILIBRIUM_AND_CROSS_SEED_MATCH_WITH_GENERATOR_FLOQUET_AGREEMENT**

Instantaneous source RHS is authenticated against the v18.51s ancestor, and the candidate is SHA-pinned to v18.52z.

The local algebraic root uses Delta=pi, x=1/2, fixed J, and a late negative external clock. The parent map and this root are independent computations.

## original

- Radial equilibrium r1=1.7247943535035, r2=0.11823231263973; omega1=-0.48446203570555, omega2=-0.015537964294452
- Candidate instantaneous RHS residual 1.00257691038e-09
- Independently rooted RHS residual 4.10782519111e-15
- Candidate delta-minus-pi residual 1.3900882896e-09
- Relative Floquet-vs-exp(local generator) max error 4.37882117622e-05
- Instantaneous algebraic equilibrium certified: True
- Local spectral comparison passed: True

## autonomous

- Radial equilibrium r1=1.7247943535035, r2=0.11823231263973; omega1=-0.48446203570555, omega2=-0.015537964294452
- Candidate instantaneous RHS residual 1.00888997157e-09
- Independently rooted RHS residual 2.69784194984e-14
- Candidate delta-minus-pi residual 1.3962948804e-09
- Relative Floquet-vs-exp(local generator) max error 4.37890485275e-05
- Instantaneous algebraic equilibrium certified: True
- Local spectral comparison passed: True

## Cross-seed equivalence
{
  "J_difference": 1.532107773982716e-14,
  "intrinsic_locked_chart_difference": 3.802412560270631e-15,
  "candidate_phase_offset_invariance_defect": -6.206701819166938e-12,
  "phase_external_difference": -2.103025815713977e-05,
  "candidate_alpha_difference": -2.103026436384159e-05,
  "equivalent_on_fixed_charge_and_phase_frame": true
}

## Strict interpretation
The result is conditional on the *constant negative externally prescribed return drive*, exact fixed J and quotient of counter-rotation. It is neither a full historical forward-return cycle nor a global conservative Hamiltonian statement.
