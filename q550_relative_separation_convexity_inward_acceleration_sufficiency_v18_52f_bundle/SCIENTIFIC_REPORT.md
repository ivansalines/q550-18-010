# Q550 v18.52f — Local numerical audit and server replay protocol

The authenticated v18.52e physical-vector archive was inspected locally.
The new direct-vector convexity diagnostic passed **48/48** numerical checks on eight physical kick traces (4,104 sample sections); the synthetic adversarial test passed **16/16**.

Across the sampled horizon, the maximum ratio `-delta Q·delta a / ||delta V||²` is approximately **0.14139943**, less than the convexity reversal threshold `1`. The instantaneous frozen-state inward-factor requirement at the final samples is approximately **7.07–7.08**. The sampled position separation rate stays positive after the initial coincident position. The physical position return is **not observed**.

These results are sampled-domain tests, not continuous interval certificates and not a physical model of amplified force. The same executable scripts and authenticated parent archive are supplied to the user's server for reproducible replay; its results will be written to `results/SCIENTIFIC_REPORT_SERVER.md`. No new nonlinear integrations are requested by this audit.
