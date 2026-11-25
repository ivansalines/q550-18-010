# Source provenance and executable genealogy

The source of every genuine Q550 integration is the SHA-pinned historical **v18.51s** parent generator, resolved and staged via `reference/v18_53f_original_run_audit.py` and `reference/parent_v18_51t_run_audit.py`. These exact reference scripts are included (not reconstructed by guesswork). The source hash expected for the v18.51s sendback is `553775d6f46c613b6be384fb39abef81a5ab429a3188f691c125200a185d2c8b`.

Native v18.52l physical 10D states are external (~97 MB), with SHA-256 `92f54032256777febbffc9fd37c8889cab0a92a5913e802ba053c14b80b3cf12`. The package does NOT duplicate large inputs; the server is expected to have them under the historical Q550 workspace. The authenticated parent and native content are verified **before** staging/running the new IVPs.

The v18.53h **fix1** authenticated sendback, which completed 10 genuine IVPs and resolved an earlier `clock_error` NameError, is preserved by hash-pinned copies of its report, complete output manifest, per-run ledger and principal plane diagnostics in `reference/`. The exact v18.53h fix1 script is archived as `reference/v18_53h_fix1_run_audit.py`. The new v18.53i script imports the historical verified RHS, geometry, initial-state reconciliation and numerical solver (not a rewritten source implementation). This is **functional/executable genealogy**.

The historical support field, kinetic 5D embedding and the v18.53g source-output hash lineage are copied from v18.53h fix1 and rechecked before computation. `SOURCE_MANIFEST.json` hashes package source files; `OUTPUT_MANIFEST.json` hashes sources plus results; sendback ZIP contains all files excluding intermediate `__pycache__`.

Scientific scope: same authenticated initial state and forcing per seed, eleven κ values, 26 actual IVPs including four high-accuracy refinements, default ξ to 51.2. Since native archival history ends at ξ=25.6, the rest is an explicitly labelled **extension of the original model clock prescription**, not a new observation or a reconstruction of the entire ramp. Do not conflate counterrotating conserved J with a putative conserved hinge moment.
