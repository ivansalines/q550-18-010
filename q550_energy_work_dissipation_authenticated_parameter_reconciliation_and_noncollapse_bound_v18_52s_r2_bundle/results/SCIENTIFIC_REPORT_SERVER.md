# Q550 v18.52s — Energy–Work–Dissipation & Noncollapse Scientific Report

## Classification
ANALYTIC_PROOF_WITH_PARTIAL_OR_FAILED_TRACE_SUPPORT

## Exact physical-time identity
E= 1/2 Σ(rdot_i²+r_i²ω_i²)+1/2*m_x*xdot²+ΣU(r_i)+V_B(x)+κ0*x*r1*r2*cos(Δ).
dE/dt = -γρ Σrdot_i²-γx xdot² + r1*r2*(μx*sinΔ+κ0*cosΔ)*xdot + κ0*x²*r1*r2*sinΔ.
Δ=θ1+θ2−φ0−Θ; Θdot=x. The last term is a phase-drive work rate, NOT dissipation.

## Rotating frame
H=E - (x/2)(J1+J2), Hdot= -γρ Σrdot_i²-γx xdot² +xdot*[r1*r2*(μx*sinΔ+κ0*cosΔ)-(J1+J2)/2].

## Finite-time theorem vs uniform global theorem
For |κ0|<2, e_reserve>0, B>0, mx>0 and nonnegative damping, let G=E+86.4.
An explicit lower bound gives G ≥ (1/10)Σr_i^4+(1-|κ0|/2)Σr_i²+(e²/2)Σr_i^(-2)+(mx/2)xdot²+B ln4+other nonnegative kinetic terms.
Thus Gdot≤K G, K=sqrt(μx²+κ0²)*(2.5+1/mx)+|κ0|/[2(1-|κ0|/2)].
Therefore on any finite regular physical-time interval G(t)≤G(0)exp(Kt), and r_i(t)≥|e|/sqrt(2G(0)exp(Kt)).
For C=J1−J2 constant, r1²+r2²≥C²/[2G(0)exp(Kt)].
**This excludes finite-time collapse of either absolute radius, conditional on continuation and the stated model domain; it does NOT yield a time-independent lower bound as t→∞.**
No theorem concerning collapse of relative 5D separations follows from this radial statement.

## Independently evaluated RHS probes
{
  "n": 160,
  "max_abs_E_balance": 3.4416913763379853e-15,
  "max_abs_H_balance": 3.552713678800501e-15,
  "max_lower_bound_violation": 0.0,
  "max_dJdiff": 1.5543122344752192e-15,
  "pass_all": true
}

## Native traces
{
  "status": "AUTHENTIC_NATIVE_AND_PHASE_VERIFIED",
  "file_sha256": "92f54032256777febbffc9fd37c8889cab0a92a5913e802ba053c14b80b3cf12",
  "phase_source_sha256": "553775d6f46c613b6be384fb39abef81a5ab429a3188f691c125200a185d2c8b",
  "n_series": 20,
  "min_count": 65537,
  "max_energy_ledger_relative_error": 2.3009612010876923e-11,
  "max_phase_rate_minus_x": 1.0058303739589887,
  "max_charge_range": 3.197442310920451e-14,
  "max_energy_upper_bound_violation": 0.0,
  "ledger_pass": true,
  "phase_clock_pass": false,
  "bound_pass": true,
  "plots_generated": [
    "original",
    "autonomous"
  ]
}

## External phase provenance
{
  "status": "AUTHENTIC_ORIGINAL_PHASE_RECONSTRUCTED",
  "parent_sha256": "553775d6f46c613b6be384fb39abef81a5ab429a3188f691c125200a185d2c8b",
  "original_engine_slope": 0.1343837359721647,
  "kappa_original": 0.16666666666666666,
  "kappa_autonomous": 0.16659961041764054,
  "kappa_shared": 0.1666331385421536,
  "xi_interval": [
    12.8,
    25.6
  ],
  "phase_interval": {
    "original": [
      -9.837654285725517,
      -20.645435602643357
    ],
    "autonomous": [
      -9.83763325546736,
      -20.6454145723852
    ]
  },
  "theta_is_external": true,
  "Y_tau_column_11_is_material_time": true,
  "parameter_authentication": "SHA256_SOURCE_AND_EXACT_SEVEN_PARAM_MATCH",
  "authenticated_model_parameters": {
    "phi0": 0.0,
    "e_reserve": 0.01,
    "gamma_rho": 0.3,
    "barrier_B": 0.3,
    "mu_x": 0.05,
    "gamma_x": 0.3,
    "m_x": 1.0
  }
}

## Provenance
Historical RHS: reference/test_support_q550_oriented_coupled_bulge_cone_audit_v1.py; exact SHA 10D native v18.52l in inputs_config.json.
No nonlinear ODE integration was performed by this audit.
