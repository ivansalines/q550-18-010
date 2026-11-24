# Q550 v18.52d — Pre-registered scientific report

**Status: server physics pending.**

Authenticated v18.51r paired 10D trajectories already support a lightweight physical-clock-matched kinematic comparison. The archived-velocity derivative indicates that near xi=0.20, the difference in acceleration between the original and autonomous kappa0 choices is predominantly lateral, approximately 98.5% in Cartesian norm for both seed families. This result is a trajectory derivative estimate, not independent nonlinear RHS acceleration.

The v18.52c sendback, used only as supporting provenance, completed 24 midpoint physical branches and found sampled positive signed turning margins, with a conditional rather than rigorous continuum bound.

The server stage will independently compute the nonlinear acceleration on each stored state, compare the two methods and test *actual finite initial velocity kicks* at common fixed kappa0. The classification of physical recovery must remain open until the kick trajectories are available. The mere magnitude of a lateral acceleration is not evidence of recovery.

The full report is generated in `results/SCIENTIFIC_REPORT_SERVER.md` and `results/report.json` by the server.
