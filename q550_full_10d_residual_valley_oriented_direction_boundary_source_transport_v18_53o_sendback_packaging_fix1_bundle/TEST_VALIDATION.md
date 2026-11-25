# Packaging fix1 — Local validation (English)

The recovery Python module was syntax-checked and tested against an isolated **synthetic packaging fixture** with the original v18.53o source manifest and original 22 MB embedded n-fix1 ancestor, but synthetic results substituted solely for packaging validation.

Verified behavior:
- The synthetic fixture contains 10 reported IVPs, four reported refinements, 846 + 846 orientation rows, 6 orthogonality rows, 12 summary rows and 6 physical case rows; values are synthetic and are **not scientific results**.
- A deliberately broken output archive already present on disk is replaced after validation.
- Original source hashes (14/14) and synthetic fixture output hashes pass.
- Resulting archive passes CRC, exact member inventory and both manifest checks.
- No call to the Q550 solver was made.
- The actual user's physical results are not accessible in this environment, so their recovered archive must be verified on the server and sent back for interpretation.
