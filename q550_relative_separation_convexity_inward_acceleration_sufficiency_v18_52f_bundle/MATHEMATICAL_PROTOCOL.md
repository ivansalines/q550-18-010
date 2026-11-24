# Mathematical protocol — v18.52f

Let `q=Q_kick-Q_reference`, `v=V_kick-V_reference`, and `a=a_kick-a_reference` be physical Cartesian differences at **matched material time tau**. The chart is the specified five-dimensional kinetic Cartesian embedding, not a new coordinate-independent geometry.

Define `D=0.5 q·q`, `G=q·v`, `P=v·v`, `I=-q·a`, and `C=P-I`. The exact kinematic identities are `D'=G` and `G'=C`. Here `I>0` means an inward projection of acceleration, `P>=0` is the outward kinetic-kinematic term, and `C` is the convexity of squared position separation.

The diagnostic sufficiency ratio is `R=I/P` only when `P>0`; `R>=1` implies `C<=0`. A convexity reversal **does not imply** decreasing position distance. The latter requires `G<0` (away from coincident positions), and even that is not global state recovery.

The instantaneous amplification factor `F_needed=P/I=1/R` (only when `I>0`) answers the arithmetic question: by what factor would the inward term have to increase **if q and v were frozen** to achieve `C=0`? It is NOT a physical prediction, modified acceleration field, or bifurcation threshold.

Compare source NPZ archived `relative_position_rate`, `position_acceleration_dot`, `second_position_separation_derivative` with independently computed values. Verify `G'` via cubic spline on interior samples; integrate `C` to predict the change in `G`, and `G` to reconstruct `D`. Document spline/quadrature errors; these checks are not interval enclosures.

An additional shape-preserving scalar PCHIP interpolation check is allowed **only as a check of the chosen interpolant**, never as certification of the underlying ODE between samples. Report every sign reversal in the sampled data separately from a genuine event located by the continuous dynamical flow.

The audit covers the exact eight finite, lateral, initial velocity kicks, two fixed physical seed families, two signs and two relative amplitudes (1e-4 and 3e-4), each over its authenticated 513-sample material-time grid. The external driving phase can differ slightly at fixed material time; this remains a finite-clock comparison and no stronger causal recovery statement follows from it.
