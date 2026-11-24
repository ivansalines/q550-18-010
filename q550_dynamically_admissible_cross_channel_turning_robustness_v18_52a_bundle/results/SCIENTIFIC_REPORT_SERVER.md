# Q550 v18.52a: Dynamically Admissible Cross-Channel Turning Robustness Audit

**Numerical audit:** 36/36 PASS
**Scientific decision:** `SUPPORTED_IN_AUTHENTIC_THREE_KAPPA_BRANCH_BAND`

## Model-defined physical perturbation

- Three already integrated kappa0 branches: 0.166599610417641, 0.166633138542154, 0.166666666666667.
- Each branch uses its own archived nonlinear Y(t), current structural Jacobian J(t), and fresh parameter source partial_kappa f(t).
- Independent channel responses are forced tangent IVP solutions; no channel mixing is used.

## Endpoint signed turning budget (rad / xi)

| Seed | Branch | Mixed | Self v2+omega2 | Other terms | Total | Min margin |
|---|---|---:|---:|---:|---:|---:|
| original | autonomous | +0.71113699 | -0.48010764 | -0.00021678945 | +0.23081256 | +0.23034383 |
| original | midpoint | +0.71112478 | -0.48009981 | -0.00021683349 | +0.23080813 | +0.2303402 |
| original | original | +0.71111256 | -0.48009198 | -0.00021687753 | +0.2308037 | +0.23033656 |
| autonomous | autonomous | +0.71130195 | -0.48013891 | -0.00021676818 | +0.23094628 | +0.23047548 |
| autonomous | midpoint | +0.71128973 | -0.48013108 | -0.00021681222 | +0.23094184 | +0.23047185 |
| autonomous | original | +0.71127751 | -0.48012325 | -0.00021685625 | +0.23093741 | +0.23046821 |

## Scientific scope

- Only the historical original-to-autonomous kappa0 band and the existing two seed preparations are tested; no general admissible-perturbation ball is certified.
- The previously integrated nonlinear trajectories are reused; only their Jacobian, parameter source and forced channel responses are recomputed.
- The signed turning decomposition is coordinate/channel dependent; source supports are anchored to the original model equations.
- The model does not establish a physical material-return section on xi in [0,0.2].
- Discrete sample signs and ordinary numerical convergence are not rigorous interval bounds on sign preservation between samples.
- The scientific robust/nonrobust decision is recorded separately from numerical audit validity.
