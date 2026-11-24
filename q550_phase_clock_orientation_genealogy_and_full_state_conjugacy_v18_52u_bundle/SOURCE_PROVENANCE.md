# Executable genealogy

Inputs: (1) v18.51s original SHA256-authenticated server sendback including v18.51r and its ancestral engine stage; (2) v18.52l native 10D states, whose 11th column is **material time**; (3) v18.52t audited source and server report bundled under `reference/` for continuity.

Original v18.51s archive is intentionally not included: it must be supplied from the user's server and must match the SHA256 in inputs_config.json. The code prints paths, digests and source snippets for staged original theta builder.

No attempt to alter source/observations, fit phase offsets, or rerun heavy ODE integrations. Scripts and tests in English; batch mode uses Python 3.13.
