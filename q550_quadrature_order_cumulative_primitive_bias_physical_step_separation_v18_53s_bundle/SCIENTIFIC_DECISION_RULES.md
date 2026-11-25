# Scientific decision rules — v18.53s

1. **Fail-fast integrity:** verify parent v18.53r ZIP CRC, all source/output hashes, scientific decision and genuine 12 physical IVPs. Preflight must independently authenticate nested model ancestry. Reject wrong or changed input hashes.
2. **No hidden physical work:** preflight/smoke/quadrature-only modes perform no new 11D integrations. Default server full run performs six genuine tighter IVPs using authenticated model functions and initial branches.
3. **Grid separation:** fail if 1x/2x/4x grids are not nested. Simpson must be applied only to equal-width subpanel pairs. Reject non-finite or non-increasing xi grids.
4. **Exact split:** the 2x/4x cumulative trapezoid identity must close to <= 5e-13 in Euclidean 10D norm.
5. **Interpretation:** a held-source quadrature result is not an independent physical history. The physical-step replay, though genuine, is a numerical intervention only. Neither single-run bias nor Richardson contrasts are validated continuum error bars.
6. **Full success:** six tighter 11D replays, historical SHA checks, finite source, charge and RHS gates, correct output shapes, finite valley metrics, output manifests and ZIP CRC. Report quantitative sensitivities instead of forcing a favorable convergence threshold.
7. **No closure:** the residual and area are derived from Euclidean 10D bivector projections, beyond xi=25.6 formal clock continuation. No full-state recurrence, spin, material memory or material holonomy certificate follows.
