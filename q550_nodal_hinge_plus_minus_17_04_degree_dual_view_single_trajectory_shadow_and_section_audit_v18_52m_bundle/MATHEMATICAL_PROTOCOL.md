# Mathematical protocol — v18.52m

Let `q(tau) in R5` be the *single* physically embedded Q550 relative displacement at matched material time `tau`. The two orthographic view maps use a fixed hinge `e1` and one hidden unit vector `h`:

```
P_+(q) = (q.e1, q.(cos(theta)*e2 + sin(theta)*h))
P_-(q) = (q.e1, q.(cos(theta)*e2 - sin(theta)*h))
theta = 17.04 degrees (pre-registered candidate)
```

The measured gate condition is `w_+(tau)*w_-(tau)<0`; its duration fraction depends on the physical time ordering of `q(tau)` and the chosen observer. A section event is a sign change of `w_+` or `w_-`, linearly interpolated between neighboring authenticated data samples; the resulting full 5D `q` vectors are retained for comparison.

Intrinsic *hinge-perpendicular* radius is `sqrt(q2^2+q3^2+q4^2+q5^2)`, invariant under a choice of orthonormal directions spanning the physical complement of the candidate Q1 hinge. It depends on the candidate hinge, and is **not** an invariant of the model under arbitrary 5D transformations.

A signed phase diagnostic is computed from `arg(q2 + i(h.dot(q3:5)))`, which is dependent on the chosen 2D transverse plane. Phase turns **cannot** be promoted to a homotopy or knot invariant of a 5D trajectory. In particular, an axis in R5 has a 4D perpendicular complement, not the 2D complement needed for ordinary 3D linking number intuition.

Scale-normalized tube-versus-cone proxies: variation coefficient of the full hinge-perpendicular radius, and the coefficient of determination of its affine fit against `|q1|`. Either proxy can be suggestive or misleading; neither alone identifies an actual tubular/conical manifold.

A projected identity at the **same time** is algebraically guaranteed:

```
w_+ + w_- = 2 cos(theta) q2
w_+ - w_- = 2 sin(theta) (h.dot(q3:5))
```

These two equations are **not a discovery**. Nontrivial observations are time-ordered section structure, control-angle sensitivity, robustness across eight trajectories and fixed hidden observer choices, and 5D posture mismatch between *distinct times*.

Section-candidate observables include the coefficient of variation of `D` and the 5D hinge-perpendicular radius at successive tilted crossings, and consecutive section return pairs. These quantities are **not** asserted to be dynamically conserved. The same-observer successive-pair table reports 5D same-posture residual `||q_b-q_a||/max(||q_a||,||q_b||)` and antipodal residual `||q_b+q_a||/max(||q_a||,||q_b||)`.
