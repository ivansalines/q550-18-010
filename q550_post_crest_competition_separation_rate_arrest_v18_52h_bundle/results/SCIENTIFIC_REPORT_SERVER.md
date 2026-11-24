# Q550 v18.52h — Scientific server sendback

Exact parent g archive SHA-256: `3a22ed6254608dba1b623653a9879a16381c17141f7054e54976ea5b4a42535e`.
Original nonlinear model resumed at xi=0.8 to xi=1.6.
Numerical checks: 40/40; all passed: **True**.

| Kick | New C zeros | G zeros | D maxima | Return R<1 | Terminal G | Classification |
|---|---:|---:|---:|---:|---:|---|
| original_negative_1.0e-04 | 1 | 1 | 1 | 1 | -3.67186296e-09 | CONFIRMED_NUMERICAL_D_MAXIMUM |
| original_positive_1.0e-04 | 1 | 1 | 1 | 1 | -3.67475709e-09 | CONFIRMED_NUMERICAL_D_MAXIMUM |
| original_negative_3.0e-04 | 1 | 1 | 1 | 1 | -3.3021011e-08 | CONFIRMED_NUMERICAL_D_MAXIMUM |
| original_positive_3.0e-04 | 1 | 1 | 1 | 1 | -3.30991515e-08 | CONFIRMED_NUMERICAL_D_MAXIMUM |
| autonomous_negative_1.0e-04 | 1 | 1 | 1 | 1 | -3.67035399e-09 | CONFIRMED_NUMERICAL_D_MAXIMUM |
| autonomous_positive_1.0e-04 | 1 | 1 | 1 | 1 | -3.67324482e-09 | CONFIRMED_NUMERICAL_D_MAXIMUM |
| autonomous_negative_3.0e-04 | 1 | 1 | 1 | 1 | -3.30074604e-08 | CONFIRMED_NUMERICAL_D_MAXIMUM |
| autonomous_positive_3.0e-04 | 1 | 1 | 1 | 1 | -3.30855117e-08 | CONFIRMED_NUMERICAL_D_MAXIMUM |

A local maximum of squared separation D requires a positive-to-negative zero of G. A new C zero, R crest or reduced positive G alone does not establish return.
The run is a numerical finite-horizon audit with independent refinements, not a rigorous interval enclosure.
If any numerical check fails, do not promote the physical interpretation.
