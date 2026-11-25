# v18.53g — Mathematical protocol

Let the archived material states be `Y(tau)` and let the authentic external phase be `Theta(tau)` recovered from v18.51s for each historical seed. The parent generator is verified at κ_shared.

For each fixed κ ∈ {κ_aut, κ_shared, κ_orig}, compute the six 5D acceleration terms from the *same* `(Y,Theta)`:

`A_U=-U'(r_i)e_{r_i}`; `A_C=-κ*x*r_j*cos(Delta)e_{r_i}+κ*x*r_j*sin(Delta)e_{theta_i}`; four additional channels as in authenticated v18.53f. These vectors are source decomposition terms, not unique Hamiltonian interactions.

For observer axes `(n_beta,h)` with `n_beta=cos(beta)n+sin(beta)e`, `h` fixed and `(e,n,h)` orthonormal, define the physically fixed-plane integrated channel torques

`P_beta=∫ wedge(Q,A_U;n_beta,h) dtau`, `C_beta=∫ wedge(Q,A_C;n_beta,h) dtau`.

`S_beta=1-|P_beta+C_beta|/(|P_beta|+|C_beta|)` lies in [0,1] (set to null when its denominator vanishes). `S_beta` concerns two integrated channels only; separate other four channels and the final torque. Track pointwise opposite-sign fraction and separately integrated magnitudes, to avoid confusing temporal cancellation with channel cancellation.

For observer independence, let `M_F=∫ (Q⊗F-F⊗Q) dtau`, an antisymmetric 5x5 torque bivector. Then `P_beta=n_beta^T M_U h` and `C_beta=n_beta^T M_C h`. The Euclidean Frobenius inner-product cosine `cos(M_U,M_C)` and vector norm score `S_5D=1-||M_U+M_C||_F/(||M_U||_F+||M_C||_F)` are invariant under *simultaneous* global orthogonal coordinate changes, but they are not guaranteed large by one projected cancellation.

RHS verification: evaluate original authenticated parent `rhs_from_original` at each κ, keeping archived Y and Theta fixed; convert from xi to material time by x/omega, compare all 10 derivatives to the independent documented material RHS. Replacing κ does not make an archived orbit a solution for that κ; only shared κ is dynamically authentic on these states.

Independent genealogy: check against the historical v18.53f signed integrated radial and coupling torques at κ_shared, beta=0 for all 20 series. If this fails, block classification.

The frozen-external-state clock sign is not a direct seventh physical force: differing clock histories can alter subsequent phases and motion. No self-consistent paired-κ trajectories and no Noether hinge symmetry are claimed.
