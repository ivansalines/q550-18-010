# Q550 v18.53x — Local preflight validation

The preflight is a combination of exact provenance checks, synthetic identities, and a **limited authenticated replay**. No new Q550 11D or 21D integration was run locally.

- Exact nested ancestry verification: W source 13, output 7; V source 14, output 10; U source 14, output 8; T source 13, output 7. All CRC and SHA-256 checks passed.
- Exact historical 5D metric helper recovered and pinned byte-for-byte.
- SO(5)-induced 10D and 45D orthogonality, constant-gauge polygon covariance, de-gauging and moving-gauge parallel alignment tests: PASS.
- Packaging regression using a deliberately minimal synthetic result tree: PASS, with correct pathlib entry naming.
- Authenticated 3-node smoke replay (lambda 0, 1, 4, original/shared/preformal): PASS.
- Full 141-point moving-center test at lambda=4, original/shared/preformal: archived weighted 45D polygon replay error `5.8442900945e-14`.
- Authenticated lambda=4 original/shared/preformal cross-lambda 45D angles in chosen gauges: symmetric `6.5594643613` degrees, intraplane `15.1564004229` degrees, mixed `28.1128761036` degrees, countermixed `23.0055453678` degrees. Canonical physical-frame path angle `0.6438971214` degrees. These are *single-case local replay observations* and not the full 204-node server result.

The full audit (204 parent nodes, 816 gauge rows, 204 physical-frame rows, 60 moving-center physical-frame curves and 24 time-dependent-gauge checks) is reserved for the user's server. No continuum/matter holonomy conclusion was made by this validation.
