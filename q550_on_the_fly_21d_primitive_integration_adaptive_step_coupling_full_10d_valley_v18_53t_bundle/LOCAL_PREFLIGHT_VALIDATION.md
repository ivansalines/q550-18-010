# Local preflight — v18.53t

- Historical v18.53s server sendback SHA256: `3bace60f8b02cfdcb1b423965fce375cc15b0d82452c2c43ee67e3c28362bb9a`.
- Genuine s parent: 12 SHA256 source entries, 9 result entries, no CRC failures.
- Embedded r parent: 12 source entries and 10 result entries, no CRC failures.
- Executable ancestry checked down to v18.53h via the inherited `require_ancestry` logic; intermediate n/m checkpoints authenticated.
- Synthetic test: physical RHS independence from C, clock/source conversion, analytic constant-source integral, directed geodesic arclength.
- Additional historical-source test compares new source evaluation to `B_U+B_C` from the genuine historical force implementation on synthetic 11D states; PASS.
- Packaging self-test: original type-mismatch bug (`str`/`PosixPath`) is not present; CRC and SHA verification pass on a synthetic fixture.
- No new physical 11D, 21D or 10D integration executed locally.

Actual server results are not included here; numerical conclusions require the genuine server sendback.
