# v18.53v fix1 — Scientific and operational changelog

**Cause:** the final historical-parent arclength regression compared different origin-centered curves, not different numerical solutions to the same minimization. The historical v18.53u scan used the v18.53t center, but v18.53v re-centered the curve at the optimized v18.53u canonical minimum. The v18.53v regression collapsed all minimum and arclength fields into one absolute threshold (2e-4), leading to a false `PARENT_U_REGRESSION_BREACH`.

**Correction:** no threshold relaxation. Preserve historical archived arclength, independently verify it from authenticated U NPZ, and recompute the historical U width-search on the v18.53v curve's exact observation interval. Then gate comparable arclengths at 1e-5 rad and minimum fields at the original 2e-4 tolerance. Store raw unmatched and matched comparisons separately.

**No rerun:** 204 metric nodes and the 28,764 curve rows are used as already computed; no new Q550 physical IVP. This patch only recomputes the 36 historical reference curves and packages the existing results.
