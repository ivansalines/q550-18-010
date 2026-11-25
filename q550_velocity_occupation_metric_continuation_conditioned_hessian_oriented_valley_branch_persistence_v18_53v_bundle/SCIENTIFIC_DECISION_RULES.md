# Scientific decision rules — Q550 v18.53v

1. **Authenticate source**: verify parent v18.53u bytewise SHA-256, all 14/8 source-output manifest items, embedded v18.53t bytewise SHA and all 13/7 manifest entries, six expected 8961-point 11D/10D physical traces, and historical metric helper identity. Fail on any mismatch.
2. **Canonical control**: both epoch families equal the original Euclidean kinetic metric at lambda=0; verify archived v18.53u lambda=0,1,3 results within declared tolerance. Treat identical lambda=0 cases as the same metric, not two independent experimental observations.
3. **Continuation screen**: the minimizer is internal, gradient below `1e-9`, and the *scaled* 2x2 Hessian is positive with minimal eigenvalue above `1e-8`. Check independent multi-start competitors at every other node and key historical lambdas.
4. **Sensitivity screen**: compare analytical \(dp/d\lambda\) and \(d\|N\|_H/d\lambda\) to centered inter-node derivatives. Discrepancies are documented rather than silently used to discard cases.
5. **Oriented geometry**: compare 141-point minimized valley, length and 45D polygon in a clearly specified metric-square-root convention. Distinguish changes of numerical metric from physical rotations.
6. **Classify conservatively**: A smooth, conditioned branch in the two selected families supports *numeric metric-conditioned continuation* only. A competing lower minimum, loss of Hessian positivity, boundary minimum or stationarity failure triggers a warning classification.
7. **Limitations**: far endpoints lie in formal-clock continuation beyond xi=25.6. No equality of all metrics, no exact continuum certificate, no return of full 11D state, no material memory and no material holonomy are claimed.
