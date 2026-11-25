# v18.53r — Mathematical protocol

Let the original Q550 11D state satisfy z'(ξ)=F(z,ξ;κ₀) with the authenticated v18.53o equations, forces, and initial branches. Extract the model force bivector B(ξ) and the physical clock t(ξ). Define g(ξ)=B(ξ)t'(ξ), with ten independent bivector coordinates, and C(ξ)=∫g(ξ)dξ using the archived trapezoidal accumulation convention. For any three equal contiguous windows, their total residual is N(a,w)=C(a+3w)-C(a).

The valley solves min_w ½||N(a,w)||² in a prescribed, fixed width band. The interior stationary-width equation is N·g(end)=0. Compare the unit directions u=N/||N|| and their oriented forms u∧du/dξ, whose 45 components lie in Λ²R¹⁰. Direct polygonal arc length is Σ 2 arctan2(||uᵢ₊₁-uᵢ||,||uᵢ₊₁+uᵢ||); tangent arc length uses ∫||du/dξ||dξ.

Densification must refine both the grid at which physical forces are *evaluated* and the maximum internal DOP853 step. An unchanged adaptive solver with only denser t_eval would not demonstrate independent time-step convergence. The baseline 1× samples are the original authenticated archive, and 2× and 4× are new physical integrations from the same state. Preserve exact overlapping knot coordinates. Compute differences of physical states, integrated primitives, sample-source derivatives, widths, norms, angular speed, direct length, oriented area. Do not infer hard error bounds from successive differences.

An observed minimum with nonzero full-10D residual is not a closed state in R¹¹. A projected or formal-clock near-return is not a material return or holonomy certificate.
