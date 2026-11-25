# Mathematical protocol — Q550 v18.53i

All equations and symbols in this audit are intentionally documented in English.

## 1. Source and parameterization

Authentic source: v18.51s nested generator recovered by SHA-pinned v18.53h fix1 `stage_original`, with native 10D initial conditions independently SHA-authenticated from v18.52l. Physical state `Y=(ln r1,v1,theta1,omega1,ln r2,v2,theta2,omega2,z,vz)`; `x=sigmoid(z)`; appended native coordinate `tau` is material time, NOT external phase. On an IVP in the source parameter ξ, `d tau/d xi = omega_engine/x` (checked on native source). The externally forced angular phase is `Theta_ext(xi)=context.theta(context.event+xi)`, and the signed external clock uses the original source's `cubic_clock_sign` prescription. The source's 10D ξ-RHS and 5D material acceleration are mutually checked. The model source parameter κ enters phase coupling and active force channels as documented in `reference/test_support_q550_oriented_coupled_bulge_cone_audit_v1.py` and `scripts/torque_math.py`.

There are three historical κ prescriptions: `κ_original=1/6`, `κ_autonomous≈0.16659961041764054`, and `κ_shared≈0.1666331385421536`. For sensitivity, evaluate a declared fixed sweep of 11 κ values: all three plus shared ±{0.25,1,3,10} times Δκ with Δκ=κ_original−κ_autonomous. κ is constant within each IVP. No distinct-κ pre-12.8 ramp is implied.

## 2. The true late-horizon IVP

Start at one of two historical native states at ξ=12.8 (original/autonomous seed family). Integrate all κ choices with the **same initial 10D Y and tau** and the **same external phase/drive function** within each family, using the independently staged original source RHS and DOP853. The archive covers [12.8,25.6]. Continue the source forcing to default ξ=51.2. The interval after 25.6 is explicitly **continued-model protocol**, not an archive-observed historical ramp.

Require exact v18.53h fix1 source and output SHA checks; require original native material clock agreement; reproduce original-v18.53h singular strength and 5D cancellation at ξ=25.6 on identical quadrature grid. Require source numerical residual, charge conservation and tighter integration on four separate κ IVPs. Use a baseline shared-κ native path replay test.

## 3. Residual torque as a full 5D skew tensor

`Q=(r1 cos theta1,r1 sin theta1,r2 cos theta2,r2 sin theta2,sqrt(m_x)x)` and material 5D acceleration force channels from source. Define skew tensor

`B(Q,F) = Q F^T - F Q^T`.

For any interval [a,b], set

`U[a,b]=∫_a^b B(Q,F_radial_potential) d tau`,
`C[a,b]=∫_a^b B(Q,F_phase_coupling) d tau`,
`R[a,b]=U[a,b]+C[a,b]`.

Cancel score `S=1-||R||_F/(||U||_F+||C||_F)` when denominator is nonzero. Signed integrated cancellation does NOT imply instantaneous opposition or full six-force cancellation.

Any real 5×5 skew tensor has two rank-two blocks with strengths σ1≥σ2≥0 and a null direction. Real Schur decomposition provides orthonormal rank-two plane projectors `P1,P2` that do not depend on sign choices of Schur basis vectors. Near σ1≈σ2, the individual planes lose identifiability: report degeneracy without claiming a physical flip. The combined nonzero four-dimensional subspace and the full bivector remain defined.

## 4. Parametric differential response

For a preset half-step `h=m Δκ`, `m∈{0.25,1,3,10}`:

`D_h R(κ)=[R(κ+h)−R(κ−h)]/(2h)`,
`D_h P_j=[P_j(κ+h)−P_j(κ−h)]/(2h)`.

Report `||D_h R||_F`, `||D_h P_j||_F`, strengths and cancellation derivatives, and `||R(κ+h)−2R(κ)+R(κ−h)||_F/||R(κ+h)−R(κ−h)||_F`. Compare the derivative magnitudes across h to diagnose conditioning. This is a finite-difference numerical diagnostic, not a differentiability theorem; beware numerical conditioning for very small h and ambiguous projectors if σ1 and σ2 nearly coincide. In addition, compare exact original/autonomous historical κ variants independently with shared baseline.

## 5. Temporal persistence

Report both disjoint windows [12.8,25.6],[25.6,38.4],[38.4,51.2] and cumulative windows [12.8,25.6],[12.8,38.4],[12.8,51.2]. For each truly integrated path, compare the consecutive-segment projectors and whole skew tensors; do NOT infer a continuous local plane rotation merely from the cumulative-integral comparisons. For every κ and fixed time interval, compare the full 5D bivector cosine and the two individual principal-plane maximum angles to shared κ. The threshold `cos>0.99` and `each max angle<10 degrees` is a **predeclared strong-alignment screening threshold** only; a failure is a valid scientific result, not a numerical source-fidelity failure.

## 6. Limits

No global charge-independent stability, conserved hinge Noether charge, complete historical ramp replay, formal interval-arithmetic proof, absolute-trajectory closure, or unique observer-plane symmetry is asserted. Any extended-window instability or source-clock extrapolation issue must be recorded, not silently replaced with an apparently stable path. The archive time end is ξ=25.6, regardless of the numerical continuation endpoint.
