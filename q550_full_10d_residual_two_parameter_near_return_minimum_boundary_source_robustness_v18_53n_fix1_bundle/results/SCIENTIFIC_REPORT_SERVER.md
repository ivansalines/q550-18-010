# Q550 v18.53n-fix1 — Full 10D Residual Minimum Audit

**Scientific decision:** `AUTHENTIC_FULL_10D_NONZERO_NEAR_RETURN_MINIMUM__NO_MATERIAL_CLOSURE`

Actual Q550 physical IVPs: 10

Full bivector near-return is not a full material/physical orbit return.

```json
{
  "n_grid_rows": 13674,
  "n_multistart_runs": 132,
  "n_distinct_local_minima": 16,
  "min_full_net_norm": 0.0786780691436737,
  "best_case": {
    "seed": "autonomous",
    "kappa_rule": "original",
    "origin_xi": 14.327570426977106,
    "width_xi": 12.799954323176044,
    "end_xi": 52.72743339650523
  },
  "min_parent_root_replay_error": 2.7755575615628914e-17,
  "best_full_vs_m_projected_gains": [
    {
      "seed": "original",
      "kappa_rule": "shared",
      "gain": 0.00202700268976648
    },
    {
      "seed": "original",
      "kappa_rule": "original",
      "gain": 0.0020483883343749842
    },
    {
      "seed": "original",
      "kappa_rule": "autonomous",
      "gain": 0.002005801030294213
    },
    {
      "seed": "autonomous",
      "kappa_rule": "shared",
      "gain": 0.0020322591153932312
    },
    {
      "seed": "autonomous",
      "kappa_rule": "original",
      "gain": 0.0020537112178807826
    },
    {
      "seed": "autonomous",
      "kappa_rule": "autonomous",
      "gain": 0.002010991132429918
    }
  ],
  "max_pchip_reoptimized_origin_width_scaled_shift": 0.0006704579964469901,
  "max_quadrature_minimum_parameter_shift": 0.00014358022866227678,
  "max_refined_minimum_parameter_shift": 3.786092813687234e-12,
  "all_bests_in_formal_continuation": true
}
```

Error: none
