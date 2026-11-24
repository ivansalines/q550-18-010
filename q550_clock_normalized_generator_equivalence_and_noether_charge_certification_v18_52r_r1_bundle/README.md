# Q550 — v18.52r-r1: Clock-Normalized Generator Equivalence & Noether-Charge Certification

**Packaging correction v18.52r-r1:** explicitly creates `results/` upon fresh ZIP extraction. Supersedes the untested v18.52r packaging; the scientific protocol is unchanged.

The parent **v18.52q-r1** demonstrated an exact angular charge for archived Q550 equations and
measured a normalized disagreement 0.865616264... between the archived **support-clock** RHS
and the authenticated `structural_jet` original. The v18.52r checks whether this disagreement
is due to a **constant positive clock scaling**, independent of state, seed, phase, and kappa.

## Server (no ODE integrations)

```bash
cd /root/elab1/hysteresis
unzip q550_clock_normalized_generator_equivalence_and_noether_charge_certification_v18_52r_r1_bundle.zip
cd q550_clock_normalized_generator_equivalence_and_noether_charge_certification_v18_52r_r1_bundle
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh \
  2>&1 | tee ../v18_52r_server.log
```

The authenticated archive `q550_parameter_sensitivity_transport_transverse_posture_attribution_v18_51s_bundle_server_sendback.zip`
shall be in the parent workspace (as recovered and used by v18.52q-r1). The authentic `v18.52l`
folder containing `results/extended_horizon_10d_states.npz` must also be in the parent workspace.
Both are SHA-256 pinned. No missing reference is silently substituted. Explicit paths:
`bash run_on_server.sh --parent-v1851s /absolute/path/to/v18_51s...zip --states /path/to/extended_horizon_10d_states.npz`.

### Decision boundary

- `structural_jet` is evaluated on randomly generated **off-orbit** states for two authenticated
  seed contexts, 48 probes each. Its 10D field is compared to the archived 10D support-clock RHS.
- The best-fit clock multiplier is calculated pointwise and must be **positive, constant** and
  reproduce all ten components. Its value is fixed by a training subset and checked on holdout states.
- This demonstrates **numerical vector-field equivalence up to a constant positive time change**.
  The best-fit number is *not itself* an analytic clock law; original source and clock metadata are
  exported to check its provenance.
- Angular counterrotation charge `C=r1^2 omega1-r2^2 omega2` is checked by exact formal algebra (SymPy when installed, Fraction-valued exact monomial algebra otherwise),
  authenticated native trajectories, arbitrary states, and deliberate falsification synthetic tests.
- The full forced/dissipative model is **not** asserted globally Hamiltonian. Angular-sector
  Noether–d'Alembert follows conditionally because the angular kinetic and coupling potential are
  invariant under counterrotation, while the stated dissipative/forced terms have no angular work.
- No global no-collapse theorem or quantized topological charge is claimed.

Generated results: `results/report.json`, `results/SCIENTIFIC_REPORT_SERVER.md`,
`results/exact_symbolic_noether.json`, `results/clock_normalized_off_orbit_engine_comparison.csv`,
`results/native_charge_per_series.csv`, `results/original_structural_jet_source_excerpt.txt`
and `SHA256SUMS_OUTPUT.json`. The bundle packages the scripts and source genealogy in the sendback.
