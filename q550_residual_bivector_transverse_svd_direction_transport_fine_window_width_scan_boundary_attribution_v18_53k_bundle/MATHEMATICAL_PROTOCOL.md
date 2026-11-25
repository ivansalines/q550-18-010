# v18.53k — mathematical protocol (predeclared)

Let R_i(o,w)=∫[o+(i-1)w,o+iw] (B_U(ξ)+B_C(ξ))dτ. Each R is a real 5×5 skew matrix, interpreted in the Frobenius-isometric upper-triangle coordinates ρ(R)=√2*(R01,R02,R03,R04,R12,R13,R14,R23,R24,R34), in R^10.

The 3×10 row matrix M=[ρ(R₁);ρ(R₂);ρ(R₃)] has SVD UΣVᵀ, with σ₁≥σ₂≥σ₃≥0. Define the total residual N=R₁+R₂+R₃ in its 10D coordinates n. Set v₃=the third right singular vector, and **T=(n·v₃)v₃**, **P=n−T**. Neither T nor the projector v₃v₃ᵀ depends on the arbitrary sign chosen for v₃. Use `||T||²/||n||²` to quantify the share in the third singular direction, whenever n≠0; use ||P||, ||T|| and gaps to prevent ratio-only artifacts.

For each independently sampled origin, continuously align signs of v₃ from anchor w=12.8 toward each side. Only publish a *conditional signed coefficient* c₃=n·v₃(aligned) when σ₃/σ₁>10⁻⁸, (σ₂−σ₃)/σ₂>0.02, and adjacent |v₃(w)·v₃(w−δ)|≥0.25. A sampled sign change is *not* a root proof or a physical flip. The projective angle arccos|v₃(w)·v₃(12.8)| is sign-independent and measures rotation of the third **10D bivector-coordinate line**.

Under a fixed trajectory, the width derivative obeys the exact endpoint identity:

    d/dw [R₁+R₂+R₃] = 3 (B_U+B_C)(o+3w) (dτ/dξ)(o+3w).

The right-hand side is *boundary evaluation*, not a new material interaction. The derivative is checked against finite differences of the 10D total residual. This is a fundamental-theorem-of-calculus check of the window procedure, not an argument for dynamical memory. For individual Ri, the endpoints contribute the corresponding weighted torque rates at o+w, o+2w, o+3w.

**Controls:** authentic j source and output SHA comparison; original 72 coarse rows and baseline replays; tau-clock/source check; DOP853 integrator tightening; every-other-grid quadrature sensitivity; direct offgrid quadrature; 5D skew isometry, O(10) orthogonal invariance, row permutation invariance; equilateral planar triangle plus synthetic small common transverse bias (which deliberately creates high transverse share without any exotic motion); continuous-gauge signed-crossing synthetic test.

**Claims explicitly forbidden:** a physical third 5D axis from v₃, an exact period, closed physical orbit, observer-independent/material holonomy, singular-vector sign change as a physical reversal, or proof of causation from an integrated window decomposition.
