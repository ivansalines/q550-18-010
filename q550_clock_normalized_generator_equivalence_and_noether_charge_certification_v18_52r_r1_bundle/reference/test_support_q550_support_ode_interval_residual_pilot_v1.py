#!/usr/bin/env python3
"""
Q=550 outward-rounded support-ODE residual pilot (v1).

Goal
----
Rigorous first brick of the exact nonlinear-orbit layer.

The already-closed cocycle uses the natural cubic spline X_hat(s) built on the
support coordinate

    s = Theta/(2*pi),

so the exact nonlinear ODE in the same coordinate is

    dX/ds = H(X) = (2*pi/x) F(X).

This script interval-encloses, on selected cubic-spline segments,

    r_s(s) = X_hat'(s) - H(X_hat(s)).

It uses:
  * the exact same audit/tail/spline construction as the closed cocycle layer;
  * mpmath interval arithmetic for outward-rounded nonlinear evaluation;
  * direct interval Horner evaluation of each cubic spline piece and derivative;
  * subdivision sweeps to expose convergence of the enclosure;
  * the actual proof-required domain:
        [first dense start, last dense start + 2*P_line],
    rather than the unused end of the 70-rotation tail.

This is a rigorous residual-enclosure PILOT for the frozen spline as a reference
curve.  It does not yet prove existence of an exact nonlinear trajectory inside
a self-consistent tube.  The next step is to combine these residual bounds with
a Jacobian/log-norm (or moving-frame) bound for the error equation.

Scope
-----
Q=550, support-coordinate nonlinear orbit validation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import mpmath as mp
import numpy as np


SCRIPT_VERSION = "2026-08-11-q550-support-ode-interval-residual-pilot-v1"

DEFAULT_AUDIT = Path("test_support_mobius_monodromy_audit.py")
DEFAULT_ORBIT = Path("support_feedback_q_curtain_cache/q_550_orbit.npz")
DEFAULT_DENSE = Path("support_q550_dense_start_span_closure_raw.npz")
DEFAULT_SCOUT = Path("support_q550_nonlinear_orbit_support_defect_scout_v0_1.json")
DEFAULT_OUTPUT = Path("support_q550_support_ode_interval_residual_pilot_v1")


def load_module(path: Path, name: str) -> Any:
    path = path.resolve()
    if not path.exists():
        raise FileNotFoundError(path)
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not import {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def parse_int_list(text: str) -> list[int]:
    vals = sorted({int(x.strip()) for x in text.split(",") if x.strip()})
    if not vals or vals[0] < 1:
        raise ValueError("expected positive integers")
    return vals


def iv_float(x: float):
    """One-ULP outward enclosure of an IEEE-754 float."""
    x = float(x)
    lo = float(np.nextafter(x, -np.inf))
    hi = float(np.nextafter(x, +np.inf))
    return mp.iv.mpf([repr(lo), repr(hi)])


def iv_const(x: float):
    return iv_float(float(x))


def iv_lo(x) -> float:
    return float(x.a)


def iv_hi(x) -> float:
    return float(x.b)


def iv_absmax(x) -> float:
    return max(abs(iv_lo(x)), abs(iv_hi(x)))


def iv_width(x) -> float:
    return iv_hi(x) - iv_lo(x)


def iv_contains(x, y: float, slack_ulps: int = 4) -> bool:
    y = float(y)
    lo = iv_lo(x)
    hi = iv_hi(x)
    for _ in range(slack_ulps):
        lo = float(np.nextafter(lo, -np.inf))
        hi = float(np.nextafter(hi, +np.inf))
    return lo <= y <= hi


def spline_component_box(coeff: np.ndarray, u0: float, u1: float):
    """Interval Horner box for c0*u^3+c1*u^2+c2*u+c3."""
    u = mp.iv.mpf([repr(float(u0)), repr(float(u1))])
    c0, c1, c2, c3 = [iv_float(float(v)) for v in coeff]
    return ((c0 * u + c1) * u + c2) * u + c3


def spline_derivative_box(coeff: np.ndarray, u0: float, u1: float):
    """Interval box for 3*c0*u^2+2*c1*u+c2."""
    u = mp.iv.mpf([repr(float(u0)), repr(float(u1))])
    c0, c1, c2, _c3 = [iv_float(float(v)) for v in coeff]
    return (iv_const(3.0) * c0 * u + iv_const(2.0) * c1) * u + c2


def dU_iv(rho, e_reserve: float):
    alpha = iv_const(1.0)
    beta = iv_const(0.8)
    gamma = iv_const(0.2)
    e = iv_const(e_reserve)
    return (
        iv_const(2.0) * alpha * rho
        - iv_const(3.0) * beta * rho**2
        + iv_const(4.0) * gamma * rho**3
        - e**2 / rho**3
    )


def dV_iv(x, barrier_B: float):
    B = iv_const(barrier_B)
    return -B * (iv_const(1.0) / x - iv_const(1.0) / (iv_const(1.0) - x))


def support_rhs_iv(state: list[Any], *, kappa0: float, params: Any) -> list[Any]:
    (
        y1, v1, theta1, omega1,
        y2, v2, theta2, omega2,
        z, vz, Theta,
    ) = state

    rho1 = mp.iv.exp(y1)
    rho2 = mp.iv.exp(y2)
    x = iv_const(1.0) / (iv_const(1.0) + mp.iv.exp(-z))
    g = x * (iv_const(1.0) - x)

    if iv_lo(x) <= 0.0 or iv_hi(x) >= 1.0:
        raise FloatingPointError(f"x interval left (0,1): [{iv_lo(x)}, {iv_hi(x)}]")
    if iv_lo(rho1) <= 0.0 or iv_lo(rho2) <= 0.0:
        raise FloatingPointError("rho interval reached nonpositive values")

    phi0 = iv_const(params.phi0)
    delta = theta1 + theta2 - phi0 - Theta
    sin_delta = mp.iv.sin(delta)
    cos_delta = mp.iv.cos(delta)

    k0 = iv_const(kappa0)
    kappa_eff = k0 * x
    ratio21 = rho2 / rho1
    ratio12 = rho1 / rho2

    gamma_rho = iv_const(params.gamma_rho)
    gamma_x = iv_const(params.gamma_x)
    mu_x = iv_const(params.mu_x)
    mass_x = iv_const(params.m_x)

    F0 = v1
    F1 = (
        omega1**2
        - dU_iv(rho1, params.e_reserve) / rho1
        - kappa_eff * ratio21 * cos_delta
        - gamma_rho * v1
        - v1**2
    )
    F2 = omega1
    F3 = kappa_eff * ratio21 * sin_delta - iv_const(2.0) * v1 * omega1

    F4 = v2
    F5 = (
        omega2**2
        - dU_iv(rho2, params.e_reserve) / rho2
        - kappa_eff * ratio12 * cos_delta
        - gamma_rho * v2
        - v2**2
    )
    F6 = omega2
    F7 = kappa_eff * ratio12 * sin_delta - iv_const(2.0) * v2 * omega2

    F8 = vz
    activity = rho1 * rho2 * sin_delta
    x_dot = g * vz
    x_ddot = (
        -dV_iv(x, params.barrier_B)
        + mu_x * activity
        - gamma_x * x_dot
    ) / mass_x
    F9 = x_ddot / g - (iv_const(1.0) - iv_const(2.0) * x) * vz**2
    F10 = x

    factor = (iv_const(2.0) * mp.iv.pi) / x
    F = [F0, F1, F2, F3, F4, F5, F6, F7, F8, F9, F10]
    H = [factor * f for f in F]

    # Exact algebraic cancellation in the support coordinate:
    # (2*pi/x) * Theta_dot = (2*pi/x) * x = 2*pi.
    H[10] = iv_const(2.0) * mp.iv.pi
    return H


def model_derivative_float(state: np.ndarray, *, kappa0: float, params: Any, audit: Any):
    (
        y1, v1, theta1, omega1,
        y2, v2, theta2, omega2,
        z, vz, Theta,
    ) = np.asarray(state, dtype=float)

    rho1 = math.exp(y1)
    rho2 = math.exp(y2)
    x = float(audit.sigmoid(z))
    g = x * (1.0 - x)
    delta = theta1 + theta2 - params.phi0 - Theta
    sd, cd = math.sin(delta), math.cos(delta)
    ke = kappa0 * x

    rho1_dot = rho1 * v1
    rho2_dot = rho2 * v2
    rho1_ddot = (
        rho1 * omega1**2
        - audit.dU(rho1, params.e_reserve)
        - ke * rho2 * cd
        - params.gamma_rho * rho1_dot
    )
    rho2_ddot = (
        rho2 * omega2**2
        - audit.dU(rho2, params.e_reserve)
        - ke * rho1 * cd
        - params.gamma_rho * rho2_dot
    )
    v1_dot = rho1_ddot / rho1 - v1**2
    v2_dot = rho2_ddot / rho2 - v2**2
    omega1_dot = ke * rho2 * sd / rho1 - 2.0 * v1 * omega1
    omega2_dot = ke * rho1 * sd / rho2 - 2.0 * v2 * omega2

    activity = rho1 * rho2 * sd
    x_dot = g * vz
    x_ddot = (
        -audit.dV_barrier(x, params.barrier_B)
        + params.mu_x * activity
        - params.gamma_x * x_dot
    ) / params.m_x
    vz_dot = x_ddot / g - (1.0 - 2.0 * x) * vz**2

    return np.array(
        [v1, v1_dot, omega1, omega1_dot,
         v2, v2_dot, omega2, omega2_dot,
         vz, vz_dot, x],
        dtype=float,
    )


def support_rhs_float(state: np.ndarray, *, kappa0: float, params: Any, audit: Any):
    x = float(audit.sigmoid(float(state[8])))
    return (2.0 * math.pi / x) * model_derivative_float(
        state, kappa0=kappa0, params=params, audit=audit
    )


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fields: list[str] = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def select_segments(
    *,
    scout_json: Path,
    required_start: float,
    required_end: float,
    explicit: str | None,
    top_k: int,
) -> tuple[list[int], dict[int, float]]:
    payload = json.loads(scout_json.read_text())
    hard = payload["top_hard_support_segments"]

    scout_map: dict[int, float] = {
        int(row["segment_index"]): float(row["max_defect_norm"])
        for row in hard
    }

    if explicit:
        selected = parse_int_list(explicit)
        return selected, scout_map

    relevant = [
        row for row in hard
        if float(row["right"]) >= required_start
        and float(row["left"]) <= required_end
    ]
    selected = [int(row["segment_index"]) for row in relevant[:top_k]]
    if len(selected) < top_k:
        # The scout JSON stores only top hard segments. Fall back to what it has.
        pass
    return selected, scout_map


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--audit-script", type=Path, default=DEFAULT_AUDIT)
    p.add_argument("--orbit", type=Path, default=DEFAULT_ORBIT)
    p.add_argument("--dense-raw", type=Path, default=DEFAULT_DENSE)
    p.add_argument("--scout-json", type=Path, default=DEFAULT_SCOUT)
    p.add_argument("--q", type=float, default=550.0)
    p.add_argument("--gamma-rho", type=float, default=0.3)
    p.add_argument("--tail-periods", type=int, default=70)
    p.add_argument("--steps-per-period", type=int, default=160)
    p.add_argument("--samples-per-rotation", type=int, default=720)
    p.add_argument("--initial-x-perturbation", type=float, default=1e-6)
    p.add_argument("--segments", default=None,
                   help="explicit comma-separated spline segment indices")
    p.add_argument("--top-k", type=int, default=5)
    p.add_argument("--subdivision-sweep", default="1,2,4,8,16,32")
    p.add_argument("--sample-checks-per-subcell", type=int, default=3)
    p.add_argument("--iv-dps", type=int, default=60)
    p.add_argument("--output-prefix", type=Path, default=DEFAULT_OUTPUT)
    args = p.parse_args()

    mp.iv.dps = args.iv_dps

    audit = load_module(args.audit_script, "q550_interval_residual_audit")
    params = audit.ModelParameters(gamma_rho=args.gamma_rho)
    orbit = audit.load_orbit(args.q, args.orbit)
    tail = audit.collect_uniform_tail(
        orbit,
        tail_periods=args.tail_periods,
        steps_per_period=args.steps_per_period,
        samples_per_rotation=args.samples_per_rotation,
        initial_x_perturbation=args.initial_x_perturbation,
        params=params,
        progress_every=100,
    )
    spline = audit.build_state_spline(tail)

    with np.load(args.dense_raw, allow_pickle=True) as d:
        starts = np.asarray(d["starts"], dtype=float)
        p_line = float(d["p_line"])

    required_start = float(starts[0])
    required_end = float(starts[-1] + 2.0 * p_line)
    tail_start = float(tail.s[0])
    tail_end = float(tail.s[-1])

    if required_start < tail_start or required_end > tail_end:
        raise RuntimeError("proof-required support domain is outside the reconstructed tail")

    selected, scout_map = select_segments(
        scout_json=args.scout_json,
        required_start=required_start,
        required_end=required_end,
        explicit=args.segments,
        top_k=args.top_k,
    )
    if not selected:
        raise RuntimeError("no relevant hard segments selected")

    subdivisions = parse_int_list(args.subdivision_sweep)
    kappa0 = (55.0 / 6.0) * (10.0 / args.q)

    knots = np.asarray(spline.x, dtype=float)
    coeffs = np.asarray(spline.c, dtype=float)
    if coeffs.shape != (4, knots.size - 1, 11):
        raise RuntimeError(f"unexpected spline coefficient shape {coeffs.shape}")

    print(f"Script version        = {SCRIPT_VERSION}")
    print(f"Q                     = {args.q:g}")
    print(f"audit SHA256          = {sha256(args.audit_script)}")
    print(f"orbit SHA256          = {sha256(args.orbit)}")
    print(f"dense SHA256          = {sha256(args.dense_raw)}")
    print(f"interval dps          = {args.iv_dps}")
    print(f"tail domain           = [{tail_start:.15f}, {tail_end:.15f}]")
    print(f"required proof domain = [{required_start:.15f}, {required_end:.15f}]")
    print(f"max transport length  = {2.0*p_line:.15e}")
    print(f"selected segments     = {selected}")
    print(f"subdivision sweep     = {subdivisions}")
    print("residual              = X_hat'(s) - (2*pi/x) F(X_hat(s))")
    print()

    summary_rows: list[dict[str, Any]] = []
    subcell_rows: list[dict[str, Any]] = []
    all_containment = True

    for seg in selected:
        if not (0 <= seg < knots.size - 1):
            raise IndexError(f"invalid segment {seg}")
        seg_left = max(float(knots[seg]), required_start)
        seg_right = min(float(knots[seg + 1]), required_end)
        if not seg_left < seg_right:
            print(f"segment {seg}: outside required proof domain, skipped")
            continue

        knot_left = float(knots[seg])
        full_width = float(knots[seg + 1] - knot_left)
        scout_max = scout_map.get(seg, math.nan)

        for M in subdivisions:
            max_norm_upper = 0.0
            max_component_upper = 0.0
            worst_component = -1
            worst_subcell = -1
            max_state_width = 0.0
            max_derivative_width = 0.0
            min_x_lower = 1.0
            max_x_upper = 0.0
            sample_ok = True

            for j in range(M):
                a = seg_left + (seg_right - seg_left) * j / M
                b = seg_left + (seg_right - seg_left) * (j + 1) / M
                u0 = a - knot_left
                u1 = b - knot_left

                state_box = [
                    spline_component_box(coeffs[:, seg, k], u0, u1)
                    for k in range(11)
                ]
                deriv_box = [
                    spline_derivative_box(coeffs[:, seg, k], u0, u1)
                    for k in range(11)
                ]
                rhs_box = support_rhs_iv(
                    state_box, kappa0=kappa0, params=params
                )
                residual_box = [
                    deriv_box[k] - rhs_box[k] for k in range(11)
                ]

                comp_upper = [iv_absmax(r) for r in residual_box]
                norm_upper = math.sqrt(sum(v*v for v in comp_upper))
                local_worst_component = int(np.argmax(comp_upper))

                if norm_upper > max_norm_upper:
                    max_norm_upper = norm_upper
                    worst_component = local_worst_component
                    worst_subcell = j
                max_component_upper = max(max_component_upper, max(comp_upper))
                max_state_width = max(max_state_width, *(iv_width(x) for x in state_box))
                max_derivative_width = max(
                    max_derivative_width, *(iv_width(x) for x in deriv_box)
                )

                x_box = iv_const(1.0) / (
                    iv_const(1.0) + mp.iv.exp(-state_box[8])
                )
                min_x_lower = min(min_x_lower, iv_lo(x_box))
                max_x_upper = max(max_x_upper, iv_hi(x_box))

                # Numerical containment sanity check at a few interior points.
                if args.sample_checks_per_subcell > 0:
                    for q in np.linspace(
                        a, b, args.sample_checks_per_subcell + 2
                    )[1:-1]:
                        state = np.asarray(spline(float(q)), dtype=float)
                        deriv = np.asarray(spline(float(q), 1), dtype=float)
                        target = support_rhs_float(
                            state, kappa0=kappa0, params=params, audit=audit
                        )
                        r = deriv - target
                        for k in range(11):
                            if not iv_contains(residual_box[k], float(r[k])):
                                sample_ok = False
                                break
                        if not sample_ok:
                            break

                subcell_rows.append({
                    "segment_index": seg,
                    "subdivisions": M,
                    "subcell_index": j,
                    "left": a,
                    "right": b,
                    "residual_norm_upper": norm_upper,
                    "worst_component": local_worst_component,
                    "max_component_upper": max(comp_upper),
                    "x_lower": iv_lo(x_box),
                    "x_upper": iv_hi(x_box),
                    "sample_containment_passed": sample_ok,
                })

            all_containment = all_containment and sample_ok
            inflation = (
                max_norm_upper / scout_max
                if np.isfinite(scout_max) and scout_max > 0.0
                else math.nan
            )
            summary = {
                "segment_index": seg,
                "segment_left": seg_left,
                "segment_right": seg_right,
                "subdivisions": M,
                "interval_residual_norm_upper": max_norm_upper,
                "floating_scout_max_defect": scout_max,
                "interval_to_scout_inflation": inflation,
                "maximum_component_upper": max_component_upper,
                "worst_component": worst_component,
                "worst_subcell": worst_subcell,
                "maximum_state_box_width": max_state_width,
                "maximum_derivative_box_width": max_derivative_width,
                "minimum_x_lower": min_x_lower,
                "maximum_x_upper": max_x_upper,
                "sample_containment_passed": sample_ok,
            }
            summary_rows.append(summary)

            print(
                f"seg={seg:5d} M={M:3d} "
                f"Rupper={max_norm_upper:.6e} "
                f"scout={scout_max:.6e} "
                f"infl={inflation:.3f} "
                f"comp={worst_component:2d} "
                f"x=[{min_x_lower:.6f},{max_x_upper:.6f}] "
                f"sampleOK={sample_ok}",
                flush=True,
            )

    largest = max(subdivisions)
    largest_rows = [r for r in summary_rows if r["subdivisions"] == largest]
    worst_largest = max(largest_rows, key=lambda r: r["interval_residual_norm_upper"])
    classification = (
        "OUTWARD-ROUNDED SUPPORT-ODE RESIDUAL PILOT PASSED: "
        "finite interval residual enclosures obtained on all selected proof-domain "
        "hard segments; sampled point residuals are contained. "
        "NEXT: self-consistent exact-orbit error tube."
        if all_containment
        else
        "SUPPORT-ODE INTERVAL RESIDUAL PILOT NEEDS REVIEW: sampled containment failed."
    )

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    write_csv(Path(str(prefix) + "_summary.csv"), summary_rows)
    write_csv(Path(str(prefix) + "_subcells.csv"), subcell_rows)

    payload = {
        "script_version": SCRIPT_VERSION,
        "scope": "outward-rounded interval residual pilot for exact nonlinear orbit validation",
        "q": args.q,
        "audit_sha256": sha256(args.audit_script),
        "orbit_sha256": sha256(args.orbit),
        "dense_sha256": sha256(args.dense_raw),
        "support_coordinate": "s=Theta/(2*pi)",
        "support_ode": "dX/ds=(2*pi/x)F(X)",
        "required_proof_domain": [required_start, required_end],
        "tail_domain": [tail_start, tail_end],
        "selected_segments": selected,
        "subdivision_sweep": subdivisions,
        "interval_dps": args.iv_dps,
        "summaries": summary_rows,
        "largest_subdivision_worst": worst_largest,
        "all_sample_containment_passed": all_containment,
        "classification": classification,
    }
    Path(str(prefix) + ".json").write_text(json.dumps(payload, indent=2) + "\n")

    with Path(str(prefix) + "_report.txt").open("w") as f:
        f.write(f"Script version: {SCRIPT_VERSION}\n\n")
        f.write("Q=550 OUTWARD-ROUNDED SUPPORT-ODE RESIDUAL PILOT\n")
        f.write("================================================\n\n")
        f.write(f"required proof domain = [{required_start:.15f}, {required_end:.15f}]\n")
        f.write(f"tail domain           = [{tail_start:.15f}, {tail_end:.15f}]\n")
        f.write(f"selected segments     = {selected}\n")
        f.write(f"subdivision sweep     = {subdivisions}\n")
        f.write(f"interval dps          = {args.iv_dps}\n\n")
        f.write("Residual\n--------\n")
        f.write("r_s(s) = X_hat'(s) - (2*pi/x) F(X_hat(s))\n\n")
        f.write("Summaries\n---------\n")
        for r in summary_rows:
            f.write(
                f"seg={int(r['segment_index']):5d} "
                f"M={int(r['subdivisions']):3d} "
                f"Rupper={r['interval_residual_norm_upper']:.9e} "
                f"scout={r['floating_scout_max_defect']:.9e} "
                f"infl={r['interval_to_scout_inflation']:.6f} "
                f"component={int(r['worst_component'])} "
                f"sampleOK={bool(r['sample_containment_passed'])}\n"
            )
        f.write("\n")
        f.write(f"classification = {classification}\n\n")
        f.write("Proof boundary\n--------------\n")
        f.write(
            "This pilot rigorously encloses the residual of the frozen support spline "
            "with interval arithmetic on selected hard cubic segments. It does not yet "
            "prove that an exact nonlinear solution stays in a tube around that spline. "
            "The next step is the nonlinear error equation e'=H(X_hat+e)-H(X_hat)-r_s "
            "with a self-consistent interval/log-norm or moving-frame tube.\n"
        )

    print("=" * 78)
    print(
        f"largest M={largest}: worst selected segment={int(worst_largest['segment_index'])} "
        f"Rupper={worst_largest['interval_residual_norm_upper']:.6e} "
        f"infl={worst_largest['interval_to_scout_inflation']:.3f}"
    )
    print("classification = " + classification)
    print(f"report          = {prefix}_report.txt")
    print(f"JSON            = {prefix}.json")
    print(f"summary CSV     = {prefix}_summary.csv")
    print(f"subcells CSV    = {prefix}_subcells.csv")


if __name__ == "__main__":
    main()
