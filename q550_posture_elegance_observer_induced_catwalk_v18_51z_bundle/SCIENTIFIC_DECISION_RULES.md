# v18.51z — Decision Rules

- **PASS local posture opposition** only if both anchored channel responses have strictly negative signed angular rates, while their total has strictly positive rate, over all qualified samples for both historical seeds.
- **PASS signed-turning reconstruction** only if the exact self/mixed wedge budget reproduces the total within its declared numerical tolerance.
- **PASS observer covariance** only if common SO(2) rotations preserve mutual angle and oriented area, and including the frame connection restores absolute angular rates.
- **PASS negative controls** only if omission of the observer connection creates a measurable false turning and the deliberate observer profile produces extra apparent turning reversals.
- **PASS physical-clock test** only for a strictly positive `dtau/dxi`, with numerical agreement between direct `tau` differentiation and the chain-rule result within the stated finite-difference tolerance.
- **LIMITED RESULT:** No change of physical source definition, actuator, force, energy, material-return surface or memory constitutive law is certified by these tests.

Preflight is valid for the bundled y sendback and the exact declared local xi range only. Treat near-zero angular quantities as undefined.
