# Q550 v18.53c — Dual-κ₀ Fold-Locus Continuation & Minimum-Radius Envelope Audit

**AUTHENTIC_DUAL_KAPPA_FOLD_LOCUS_AND_SAMPLED_CRITICAL_RADIUS_ENVELOPE__NOT_GLOBAL_MINIMUM_OR_INTERVAL_PROOF**

## Dual-kappa radius at the *critical relative equilibrium*

Original kappa = 0.1666666666666667, r2_fold = 0.305717064033169, s_fold = 0.926743480465969
Autonomous kappa = 0.1665996104, r2_fold = 0.305658386253367, s_fold = 0.926942190287349

**Crucial:** the two fold radii are evaluated at distinct signed drives. They do NOT form a pair at a single simultaneous drive setting.
At common drive s=-1: original r2 = 0.118241999589576, autonomous r2 = 0.118222627469916.

Sampled min of min(r1,r2): 0.305658386253367 at k=0.1665996104
Sampled max of min(r1,r2): 0.30571706403317 at k=0.166666666666667
Monotonic on sampled grid: True, derivative range: [0.8750421282616185, 0.8750636810848281]

## Scientific limitations

- This is a narrow envelope over prescribed kappa rules at a **fixed charge** and **counterfactual constant signed drive fold**. It is not a global lower bound on dynamical radii.
- Endpoints bound *sampled* extrema only if the scan is monotonic on the sampled grid. A rigorous continuous enclosure requires interval arithmetic or a validated implicit derivative sign.
- The two original/autonomous kappa values are two parameter prescriptions, not two independent seed branches at the shared anchor.
- A local saddle-node normal-form diagnostic is not a global bifurcation theorem, and the historical forward-return ramp is not simulated.
