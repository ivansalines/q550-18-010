# v18.53q — Local validation

- Recovered genuine v18.53p server sendback; outer CRC and 14/14 source SHA-256 plus 9/9 output SHA-256 checks passed.
- Recovered nested genuine v18.53o server sendback; CRC and 14/14 source plus 11/11 output SHA-256 checks passed; scientific decision and inherited physical 10 IVPs confirmed.
- Synthetic tests for cubic source reconstruction, bounded stationary-width minimization, angular geometry, and an interior helical near-return all passed.
- One-case / 11-origin replay for three reconstruction methods passed on genuine Q550 10D primitives; no new physical IVP executed.
- Source-integrated piecewise-linear source recovers the stored primitive with about 6.48e-7 maximum 10D knot deviation in inspected cases; this is a numerical consistency diagnostic, not an alternative trajectory.
- The full six-case 141/281/561-origin scan is reserved for the user server. No scientific result of v18.53q is claimed prior to that run.
- Sendback packaging syntax uses pathlib Path composition to avoid the prior string/PosixPath exception.
