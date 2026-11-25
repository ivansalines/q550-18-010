# Mathematical Protocol — v18.53q

Let C(xi) be the authenticated 10D cumulative bivector source primitive, and g_i be the sampled source dC/dxi inferred from archived 10D bivector and physical clock rate. Define N(o,w)=C(o+3w)-C(o); r=||N||; u=N/r in the Euclidean coordinate metric used previously.

Reconstructions: PCHIP interpolates C; cubic Hermite interpolates C and enforces C'(xi_i)=g_i at each knot; source-integrated quadratic uses piecewise-linear g and integrates it, starting at the authentic C(xi_0). Distinct nodal C values in the third method are a diagnostic of archived quadrature consistency, not source correction.

Set e=o+3w and f=||N||^2/2. With g=C' and h=C'', the interior stationarity equation is S=N·g(e)=0; f_w=3S; f_ww=9(||g(e)||²+N·h(e)). Its implicit origin derivative at regular stationary points is

w'(o) = -[(g(e)-g(o))·g(e)+N·h(e)]/[3(||g(e)||²+N·h(e))].

The denominator can be small or nonsmooth at branch switches, so the implicit derivative is masked and every branch candidate is retained. Bound minima are not assigned an implicit derivative.

Along a smooth width path, dN/do = g(e)(1+3w')-g(o), and du/do=(I-uu^T)(dN/do)/r. The polygonal geodesic length is Σ 2 atan2(||u_{i+1}-u_i||,||u_{i+1}+u_i||). The 45D oriented segment-area is Σ (θ/sinθ)(u_i ∧ u_{i+1}); compare to ∫u∧u' do. This is oriented geometry in a 10D coordinate space, not material holonomy.

For branch switches or solver nondifferentiability, prefer direct normalized-vector polygonal arclength. Do not infer a delta-like physical angular impulse from a derivative spike in a nonsmooth argmin curve. All comparisons are made on the same archived physical samples.
