# Q550 v18.52s-r2 — Parameter lineage repair

**Observed failure:** `Original model parameters differ from historical config` after the successful source manifest and synthetic tests.

**Root cause:** the original v18.52s used unverified historical parameter values in `inputs_config.json` and `energy_math.Params`. The authentic original generator differs in six of the seven parameters. This was scientifically important, not a mere missing Python dependency.

**Changes:** audited all seven parameters and corrected the six mismatches from the exact error diagnostic; restored the exact same authentic-engine SHA-256 requirement; wrote a machine-readable `results/parameter_audit.json` before native evaluation; updated parameter defaults, tests, reports, symbolic scope, and a-priori bound coefficient. Previous v18.52s numeric bound should not be cited as applicable to the original Q550.

**SymPy:** optional, not the cause. Algebraic identities remain documented and the synthetic/off-orbit checks still run without SymPy.

**Remaining checks:** authoritative native phase lineage and 20-series energy budget must pass on the user's server. No long IVP is rerun.
