# v18.52j mathematical protocol

At common material time tau, let q(tau)=Q_kick-Q_reference in the model-derived
Cartesian R^5 physical embedding. Set v=dq/dtau, a=dv/dtau, then

    D = 1/2 <q,q>
    G = D' = <q,v>
    C = D'' = <v,v> + <q,a> = P-I
    P = <v,v> ,  I = -<q,a>

The reference and kicked histories have separate clocks tau_ref(xi),
tau_kick(xi). Match tau via inverse clock maps, retaining xi_ref and xi_kick.
The global drive phase is xi, NOT xi-3.2. Continue the authentic 11D states
at xi=3.2, then integrate reference and kick independently with DOP853.

For q nonzero, radial velocity is v_r=<q/||q||,v>; transverse speed is
||v-v_r(q/||q||)||, and angular speed is v_t/||q||. Compute the arccos of the
unit 5D q direction against two anchors: the *fixed* historical h seam at
xi=1.6 (retrieved from i trace), and the new j seam at xi=3.2. The fixed
anchor's angle is bounded in [0,180] deg and does NOT measure winding number;
integrated angular speed additionally supplies path length on S^4, not holonomy.

G crossings classify local extrema; C crossings classify local curvature
reversals. Compare each new D-minimum q to the last i D-minimum q using the
actual physical 5D vectors, with norm, direction and absolute difference.
Failure of any comparison or insufficient orientation norm must not be
interpreted as configuration closure.
