# Environment

Tested locally with Python 3.13 or current interpreter, NumPy, SciPy, and Matplotlib. Server usage expects `python3.13`. `OPENBLAS_NUM_THREADS=1` recommended for deterministic, resource-controlled postprocessing. All operations are offline and access only the authenticated parent archive bundled in `inputs/`. No package downloads and no physical 11D/21D IVP integrations are performed.
