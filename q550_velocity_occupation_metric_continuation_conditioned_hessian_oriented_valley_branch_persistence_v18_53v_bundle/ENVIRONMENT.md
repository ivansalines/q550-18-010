# Reproduction environment

- Recommended Python 3.13; NumPy, SciPy, Matplotlib; Python's standard `zipfile`, `pathlib`, `csv`, `hashlib`, `json`.
- `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` used for predictable numerical behavior.
- The original historical 5D kinetic metric helper is authenticated by SHA256 before import; no solver implementation is reconstructed or altered.
- All model results derive from immutable historical v18.53t traces embedded under v18.53u. Re-running server v18.53v is post-processing only.
- No numerical thresholds are silently tightened after observing server results; see `inputs_config.json` and `SCIENTIFIC_DECISION_RULES.md`.
