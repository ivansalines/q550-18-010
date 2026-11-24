# Q550 v18.51x | Source Provenance

- v18.51t original nonlinear-model-based Green transport was run on the server. Its archived `green_channel_transport_states.npz` is included as a SHA256-locked input with source Jacobians, full parameter forcing, four independent 10D forced and Green response columns, and Cartesian/local projections.
- v18.51w channel-vector angle and critical-angle audit was run on the server. Its archived `channel_vector_geometric_invariants.npz`, source script, and report are included in the v18.51x input/reference snapshots.
- The **full** v18.51x archive includes the verbatim v18.51w server sendback, which itself embeds the v18.51v genealogy. The **light** archive has the authenticated numerical inputs, their fixed SHA256 hashes, and parent scripts/reports, but not the complete archival ancestry.
- v18.51x is algebraic analysis of these verified numerical inputs, not an additional physical integration. Source mixing is applied to the archived 10D forcing and independent forced/Green responses consistently.
- All source and output manifests are regenerated for this version and checked with standalone scripts. No parent result is modified in place.
