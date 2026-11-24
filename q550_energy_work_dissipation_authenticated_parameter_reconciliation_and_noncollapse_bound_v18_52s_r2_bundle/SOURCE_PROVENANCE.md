# Provenance and genealogy

1. `reference/test_support_q550_oriented_coupled_bulge_cone_audit_v1.py` is an unmodified archived model RHS from earlier Q550 work (Library source).
2. `reference/v18_52r_clock_audit.py` is retained as a source snapshot; v18.52r-r1 certified clock normalization of modern `structural_jet` against this archived RHS.
3. v18.52l 10D NPZ is not bundled (101.6 MB); expected SHA-256 `92f54032256777febbffc9fd37c8889cab0a92a5913e802ba053c14b80b3cf12` is enforced.
4. v18.52s audit scripts were created specifically to calculate model energy and the conservative upper comparison, without solving new ODEs.
5. Names/parameters and analyses refer to the *archived, documented* Q550 model. New externally forced terms or modified source terms must cause re-derivation, not silent reuse.

6. The current v18.52l NPZ stores augmented material time as column 11, **not** global phase. v18.51s original phase context is required to evaluate native model energy; reference/parent_v18_51t_run_audit.py is preserved to reconstruct it from authenticated source.
7. v18.52g used a common kappa equal to the mean of historical original/autonomous kappa values; native energy ledger must load that common kappa from original ancestor NPZ, rather than assuming a separate seed kappa.

8. Authenticated v18.52l `report.json` records `source_lineage.fixed_model_kappa=0.1666331385421536`; the reconstructed midpoint is compared to this pinned value before evaluating the native ledger.

### v18.52s-r2 parameter correction

The seven audited parameters (six changed and gamma_rho unchanged) were deduced from the exact v18.52s failure report on 2026-11-08 (`engine.params - v18.52s.model_params`). They are not inferred from a fitted energy curve. They are re-authenticated from the exact SHA-256 parent archive (`v1851s_sha256` in `inputs_config.json`) at runtime. Any discrepancy fails closed. The previously quoted v18.52s quantitative finite-time bound with `e_reserve=0.26` did not apply to the authenticated Q550 original engine.
