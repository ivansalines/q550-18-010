# Executable source genealogy

The exact source data for v18.52g are obtained from the user-supplied **v18.52f server sendback**:

`v18.52f ZIP -> embedded v18.52e ZIP -> embedded v18.52d ZIP -> historical v18.51r material-clock 10D seed states and physical kick vectors`.

The v18.52g ZIP physically contains all of the above, including the actual parent Python audit scripts (separately copied for easier inspection) and the complete user-supplied f archive. `scripts/run_audit.py` hashes f/e/d, checks their source/results manifests, checks reference scripts are byte-identical, then uses the independently SHA-256-checked v18.51s archive, whose authenticated ancestry reconstructs the original nonlinear engine.

The exact external v18.51s archive is required but is **not** in the f archive: this is an explicitly reported prerequisite, not a silent omission. Any alternate implementation must be verified to be functionally equivalent with the exact nonlinear RHS, clock and embedding and cannot silently replace an unknown source with a drawing or reduced model.

The historical window is `xi<=0.2` (513 sections), whereas new server integrations extend to a user-selected `xi_end`. Results are generated only on the server. All produced files have SHA-256 digests collected in `OUTPUT_MANIFEST.json`.
