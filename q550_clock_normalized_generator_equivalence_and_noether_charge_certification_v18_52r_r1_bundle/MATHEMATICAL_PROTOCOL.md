# Mathematical protocol — v18.52r

Let the independent state be `(l1,v1,theta1,omega1,l2,v2,theta2,omega2,z,vz,Theta)`, where
`r_i=exp(l_i)` and `x=sigmoid(z)`. The angular coupling phase is
`Delta=theta1+theta2-phi0-Theta`, `ke=kappa0*x`, with

```
omega1dot=ke*(r2/r1)*sin(Delta)-2*v1*omega1
omega2dot=ke*(r1/r2)*sin(Delta)-2*v2*omega2
l1dot=v1; l2dot=v2
```

For `C=r1^2*omega1-r2^2*omega2`:
`dC/dt = [2*r1^2*v1*omega1+r1^2*omega1dot] - [2*r2^2*v2*omega2+r2^2*omega2dot] = 0`.
Under any regular clock `t=t(s)` and `dt/ds>0`, `dC/ds=(dC/dt)(dt/ds)=0`.

The angular Lagrangian can be written locally
`L_ang=0.5*(r1^2*omega1^2+r2^2*omega2^2)-kappa0*x*r1*r2*cos(Delta)`.
It is invariant under `theta1->theta1+eps`, `theta2->theta2-eps`; its momentum map is
`C=p_theta1-p_theta2`. Other radial and auxiliary forces in the full model do not enter
this angular virtual variation. This supports a *sectoral Noether–d'Alembert* classification,
not a globally conservative full Hamiltonian statement.

For the original engine's `F_o` and archived support field `F_s`, measure
`lambda_i = (F_o dot F_s)/(F_s dot F_s)` and the defect
`||F_o-lambda_i F_s||/max(1,||F_o||,||lambda_i F_s||)`.
Do NOT assume equality from the old 0.865616... relative discrepancy. Verify that the
`lambda_i` is positive, state-independent, common to ten components, and transferable
across seed families and held-out arbitrary states. Independently inspect the clock law
in the original source before calling the proportionality an exact clock-identity.

Conservation does not imply topological integer winding, and only conditionally
excludes simultaneous `r1=r2=0` when a global bound on kinetic energy is separately established.
