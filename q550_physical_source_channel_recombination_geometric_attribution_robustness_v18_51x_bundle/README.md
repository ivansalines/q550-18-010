# Q550 — v18.51x

**Physical Source-Channel Recombination & Geometric Attribution Robustness Audit**

This audit uses the actual Q550 Green-response data from the authenticated server sendbacks v18.51t and v18.51w, held in small, SHA256-locked source snapshots. The full bundle also includes the *entire authenticated v18.51w server sendback* for explicit ancestry verification. No nonlinear or tangent ODE replay is needed in v18.51x.

## How to run (Python 3.13 recommended)

```bash
cd /root/elab1/hysteresis
unzip q550_physical_source_channel_recombination_geometric_attribution_robustness_v18_51x_bundle.zip
cd q550_physical_source_channel_recombination_geometric_attribution_robustness_v18_51x_bundle
set -o pipefail
OPENBLAS_NUM_THREADS=1 PYTHON_BIN=python3.13 bash run_on_server.sh \
  2>&1 | tee ../v18_51x_server.log
```

The light bundle contains all required numerical inputs but not the 143 MB W ancestral sendback. It runs identically on these recorded inputs and explicitly marks the full ancestral ZIP as not present. The full bundle checks the parent against CRC, source hashes, result hashes and the included W NPZ snapshot.

## Outputs

- `results/report.json`: complete scientific check ledger and per-seed metrics
- `results/SCIENTIFIC_REPORT_SERVER.md`: concise reproducible account
- `results/channel_recombination_endpoint_ledger.csv`: per-seed, per-mixing results with condition numbers and recomposition defects
- `results/channel_recombination_section_ledger.csv`: five representative mixing families across all 513 sections
- `results/channel_recombination_robustness_arrays.npz`: actual vector responses, critical mixing trajectories and endpoint scans
- `results/channel_recombination_robustness.png`: endpoint recombination angle, compensation, signed area, and critical source mixing over sections
- `results/sectional_channel_mix_response.png`: along-section angle and cancellation under selected source mixes
- `SHA256SUMS_OUTPUT.json`: output hashes; `SOURCE_MANIFEST.json`: immutable input/source hashes

## Meaning and limitations

Changing physical observation axes does not alter vector angles, while **changing the two source-channel definitions** can alter their mutual angle and projected cancellation even when the original 10D physical forcing and its response are fixed. This is a source-attribution sensitivity test, not a different material experiment. It does not invalidate the physical support of v2 and omega2 channels at the baseline, nor prove material memory.

The two seeds should agree closely on the endpoint robustness threshold near 6% (reallocation of each source into the other). Invertible redefinitions are tracked by matrix conditioning; the rank-deficient eta=0.5 case is excluded. Comments and scientific documents are in English.
