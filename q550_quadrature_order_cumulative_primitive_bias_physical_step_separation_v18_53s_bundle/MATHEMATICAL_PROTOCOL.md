# Mathematical protocol — v18.53s

Let `g(ξ)` denote the authenticated 10D source and `C(ξ)` its numerically archived cumulative primitive. The 10D residual between endpoints is `N(a,b) = C(b)-C(a)`; previous audits used `b = a+3w`.

For a fixed **4x** source history `g4`, evaluate cumulative composite trapezoid integrals on the `1x`, `2x`, and `4x` nested grids. Also compute composite Simpson on `2x` and `4x`. The `1x` grid is generally nonuniform, but within each inherited base interval the 2x and 4x subdivisions are equidistant; Simpson is applied only within those blocks. No global uniform-spacing assumption is made.

Define `T_j[g]` as the cumulative trapezoid at base-grid nodes using density `j`, and `S_j[g]` analogously for Simpson. For smooth source, test the heuristic order ratio `||T1-T2||/||T2-T4||`, as well as the contrasts `S2-S4`, `(4 T4 - T2)/3`, and `(16 S4 - S2)/15`. These are **not** certified error bounds or exact continuum values.

**Exact numerical bookkeeping:**

`T4[g4] - T2[g2] = (T4[g4]-T2[g4|2x]) + (T2[g4|2x]-T2[g2])`.

The first term changes only quadrature resolution on one fixed 4x numerical source; the second holds the 2x quadrature nodes/rule fixed and changes source values. This identity holds to rounding precision and carries no causal claim beyond its definitions.

For a separate physical-step discriminator, solve the **same original** Q550 11D ODE at the **same 4x xi nodes** using DOP853 with half the previous max physical step and tighter tolerances. Compare the resulting 11D states, 10D source, and 4x quadrature on the identical sample grid. These finite differences probe step/tolerance sensitivity; they do not isolate the true solution without a convergence or validated enclosure.

For valley geometry, use the common 281 origins and original width search band. Construct `Primitive(x,C,g,'hermite')` for each available primitive: archived 4x, composite trapezoid 4x, exact integral of a PCHIP interpolant of the sampled source 4x, and (full execution) the same source-PCHIP integral of the tighter IVP. Record directed geodesic arclength of normalized `N`, the center speed `||u'||`, and the 45D oriented area from `u wedge du`.

All measured endpoints are subject to ξ>25.6 formal-model-clock limitation. Do not interpret `u wedge du` as material holonomy. Original comments, documents and labels are intentionally in English.
