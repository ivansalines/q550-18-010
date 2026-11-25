# Q550 v18.53h — Residual-Torque Principal Planes & Dual-κ Dynamic Persistence

## Scientific decision
`AUTHENTIC_SOURCE_SELF_CONSISTENT_DUAL_KAPPA_LATE_WINDOW_DYNAMICAL_GEOMETRY__NOT_FULL_RAMP_OR_NOETHER_HINGE`

## Actual source availability
Parent: SHA256_AUTHENTICATED; Native: SHA256_AUTHENTICATED

## Scientific scope
The frozen-state v18.53g provenance is authenticated, but does not certify distinct trajectories.
On full execution, two independent native 10D IVPs per initial seed are run at κ-original and κ-autonomous, plus a shared-κ baseline replay; a stricter replay checks convergence.
The external historical phase is preserved from the authenticated context; the initial state at xi=12.8 is common between κ variants. This is NOT a complete re-integration of pre-xi=12.8 historical ramp.
Integrated skew bivectors are decomposed into two oriented principal planes. Principal-plane tracking near spectral degeneracy must be interpreted only at the combined-subspace level.
Bivector compensation score uses signed integrals and is not time-local cancellation. No globally conserved hinge charge or new bifurcation certificate is asserted.

## Metrics
{
  "n_traces": 6,
  "n_windows": 30,
  "max_charge_drift": 4.218847493575595e-15,
  "max_source_rel": 1.5491080640472887e-15,
  "max_replay_error": 9.059419880941277e-14,
  "max_refinement_error": 5.445643935786393e-14,
  "seed_pair_comparison": [
    {
      "seed": "original",
      "kappa_original": 0.16666666666666666,
      "kappa_autonomous": 0.16659961041764054,
      "kappa_gap": 6.705624902611507e-05,
      "max_5d_position_separation": 9.843153409885998e-05,
      "final_5d_position_separation": 7.629183124827695e-05,
      "max_native_10d_scaled_separation": 0.000549062148139566,
      "endpoint_5d_separation": 7.629183124827695e-05,
      "original_cancellation": 0.795235196502094,
      "autonomous_cancellation": 0.7951636409867101,
      "full_bivector_cosine": 0.999999985747583,
      "relative_residual_change": 0.00017195164131739464,
      "principal_plane_pairing_ill_conditioned": false,
      "planes": [
        {
          "mode": 1,
          "sigma_a": 0.744653479399364,
          "sigma_b": 0.7446305443319089,
          "principal_angles_deg": [
            0.0021312437406026214,
            0.007009551366331099
          ]
        },
        {
          "mode": 2,
          "sigma_a": 0.03605303052083836,
          "sigma_b": 0.03602398840064498,
          "principal_angles_deg": [
            0.0021571864700219564,
            0.12496844349984725
          ]
        }
      ],
      "sigma_a": [
        0.744653479399364,
        0.03605303052083836
      ],
      "sigma_b": [
        0.7446305443319089,
        0.03602398840064498
      ]
    },
    {
      "seed": "autonomous",
      "kappa_original": 0.16666666666666666,
      "kappa_autonomous": 0.16659961041764054,
      "kappa_gap": 6.705624902611507e-05,
      "max_5d_position_separation": 9.843679141152279e-05,
      "final_5d_position_separation": 7.623196569613361e-05,
      "max_native_10d_scaled_separation": 0.0005491763047584947,
      "endpoint_5d_separation": 7.623196569613361e-05,
      "original_cancellation": 0.7952323782427452,
      "autonomous_cancellation": 0.7951608274152475,
      "full_bivector_cosine": 0.9999999857122299,
      "relative_residual_change": 0.00017216269053833976,
      "principal_plane_pairing_ill_conditioned": false,
      "planes": [
        {
          "mode": 1,
          "sigma_a": 0.7446582092498326,
          "sigma_b": 0.7446352700351462,
          "principal_angles_deg": [
            0.002134415674161976,
            0.007014727888946693
          ]
        },
        {
          "mode": 2,
          "sigma_a": 0.03615070085906837,
          "sigma_b": 0.036121359650680834,
          "principal_angles_deg": [
            0.0021605026793364815,
            0.12474034134210484
          ]
        }
      ],
      "sigma_a": [
        0.7446582092498326,
        0.03615070085906837
      ],
      "sigma_b": [
        0.7446352700351462,
        0.036121359650680834
      ]
    }
  ],
  "sampled_all_geometric_comparisons": true
}
