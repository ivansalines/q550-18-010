# Mathematical protocol — Q550 v18.52h

All observations use the exact 5D physical Cartesian kinetic embedding `H(Y)=(Q,V)` from the authenticated v18.49h source, carried forward unchanged by the g generator. Comparison is at the same independently integrated **material time tau**, with no assumption that two trajectories have the same drive-event clock xi at that instant.

For `q(tau) = Qkick - Qref`, `v = Vkick - Vref` and `a = akick - aref`, define

`D = 1/2 q·q`, `G = q·v`, `P=v·v`, `I=−q·a`, and `C = P−I`.

The identities are `dD/dtau = G` and `dG/dtau = C`. If `P > 0`, `R=I/P` and `C/P=1−R`.

Historical facts accepted from the SHA-authenticated g server report: eight `C` positive→negative crossings around `tau≈0.788`; no detected `G` crossing through xi=0.8; a subsequent **sampled** crest of `R` around tau≈1.239. h restarts exactly at xi=0.8 using both archived coarse and fine 10D+clock terminal states, without retuning kappa or forcing. It follows the next post-crest chronology: `R` re-cresting, a potential `R=1` downcrossing (equivalently `C` negative→positive), potential arrest `G=0`, and a confirmed local D maximum only for `G` positive→negative. A positive `G` with negative `C` indicates continued outward separation at a decreasing rate, not inward return.

Roots are refined using Brent on two nonlinear dense DOP853 flows after inverting the physical clock. Distinct runs using the archived fine and coarse g endpoints assess numerical robustness. The event grid does not certify sign between sample brackets; a nonlinear solver success is a numerical result, not a rigorous enclosure of all roots. An excessively small P suppresses the ratio diagnostic; no R event is inferred then.
