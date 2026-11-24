# Q550 v18.51v — server scientific report

**Status:** PASS
**Classification:** CLOCK_REPARAMETERIZED_GREEN_AND_FRAME_COVARIANCE_NUMERICAL_DIAGNOSTIC_NOT_MATERIAL_MEMORY_CERTIFICATE

Produced by the server run; the pre-distributed scientific report is only a protocol.

## Fixed physical-clock sensitivity versus frozen-clock kernel

### original

- Physical elapsed tau: 0.334818563828944
- Clock start/end: 1.67363341943275 / 1.67493737885126
- Terminal fixed-xi Cartesian sensitivity norm: 0.0464822319008
- Terminal fixed-tau clock correction norm: 1.02344802736e-07
- Terminal correction / fixed-xi norm: 2.20180483059e-06
- Nonuniform tau quadrature discrepancy: 7.925e-14
- Direct matched-tau versus corrected formula: 6.064e-12
- Connection derivative defect (maximum across frames): 5.455e-07
- Physical source ages at final observation (radial omega2): {"mass": 0.0027233045337023347, "median_xi": 0.055859375, "median_tau": 0.0934903743865276, "mean_lag_xi": 0.13792457228278054, "mean_lag_tau": 0.23091967711526717, "first_half_xi_mass_fraction": 0.8024869006703055, "first_half_tau_mass_fraction": 0.8024869006703055}

### autonomous

- Physical elapsed tau: 0.334818495558058
- Clock start/end: 1.67363306400597 / 1.67493706702195
- Terminal fixed-xi Cartesian sensitivity norm: 0.0464817744998
- Terminal fixed-tau clock correction norm: 1.02374263965e-07
- Terminal correction / fixed-xi norm: 2.20246032056e-06
- Nonuniform tau quadrature discrepancy: 7.925e-14
- Direct matched-tau versus corrected formula: 6.625e-12
- Connection derivative defect (maximum across frames): 5.454e-07
- Physical source ages at final observation (radial omega2): {"mass": 0.002724162368311101, "median_xi": 0.055859375, "median_tau": 0.09349035458762876, "mean_lag_xi": 0.13792578614335532, "mean_lag_tau": 0.2309216625689108, "first_half_xi_mass_fraction": 0.8024958472955942, "first_half_tau_mass_fraction": 0.8024958472955942}

## Interpretation boundaries

- The frozen-history Green density K_tau=K_xi/c(source) is a reparameterization of the same fixed-xi tangent response; alone it does not imply invariance of full kappa sensitivity at fixed physical tau.
- The parameter-dependent physical clock requires a separate clock-shift correction. This is measured using already-integrated distinct-kappa trajectories and checked with direct tau-matched interpolation.
- Clock-shift attribution to individual source channels is not unique without an additional physically justified clock-source decomposition; it is retained as a separate term.
- The full 2D Cartesian response and its Euclidean norm survive a passive rotating-frame change. Radial and azimuthal component magnitudes and cancellation percentages are chart dependent.
- The observer rotation is a known passive O(2) coordinate transformation and does not establish arbitrary state-dependent gauge invariance of physical laws.
- All source forcing is an on-trajectory first variation at a single midpoint kappa0; not independently controlled physical actuators.
- No complete physical return section, interval certificate, persistent shadowing law, or constitutive material memory is certified.
- Time conversions are tied to the actual positive physical clock of the two seed trajectories over xi in [0,0.2].
