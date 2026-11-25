#!/usr/bin/env python3
"""Q=550 temporally centered non-coplanar bulge-cone audit (v1).

The previous pair-polarization audit placed the input and output SVD axes in
one active 2D chart.  That is useful for the operator, but it cannot decide
whether the two trajectory bulges occupy the same physical plane.

This audit reconstructs the two consecutive P_line arcs around the shared
center s_c = s + P_line directly from the central nonlinear orbit:

    left bulge:   q(s)          -> q(s + P_line),
    right bulge:  q(s + P_line) -> q(s + 2 P_line).

For each bulge it fits an affine rank-two orbital plane in the full reduced
9D geometry.  It then measures the two principal angles between those planes,
constructs forward cone rays through the maximum chord-deviation points, and
compares their cone axis with the true center tangent.  The relative phase is
measured both from the two bulge maxima and by a profile-correlation lag.

Only after the 9D geometry has been measured are the cone rays projected into
the active operator chart and compared with the coupled-TM SVD pair.  Thus a
2D chart can observe the cone but cannot force the two bulges to be coplanar.

This is a floating-point geometric audit.  It does not interval-certify plane
separation, cone orientation, nonlinear return, nesting, or hysteresis.
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

import matplotlib.pyplot as plt
import numpy as np


SCRIPT_VERSION = "2026-08-31-q550-oriented-coupled-bulge-cone-audit-v1"
DEFAULT_AUDIT = Path("test_support_mobius_monodromy_audit.py")
DEFAULT_ORBIT = Path("support_feedback_q_curtain_cache/q_550_orbit.npz")
DEFAULT_DENSE = Path("support_q550_dense_start_span_closure_raw.npz")
DEFAULT_STEPS = Path("support_q550_coupled_tm_distinction_lift_cell49_v1_steps.csv")
DEFAULT_OUTPUT = Path("support_q550_oriented_coupled_bulge_cone_cell49_v1")


def load_module(path: Path, name: str) -> Any:
    resolved = path.resolve()
    if not resolved.exists():
        raise FileNotFoundError(f"Required script not found: {resolved}")
    spec = importlib.util.spec_from_file_location(name, resolved)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not import {resolved}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize(vector: np.ndarray, *, label: str) -> np.ndarray:
    vector = np.asarray(vector, dtype=float)
    norm = float(np.linalg.norm(vector))
    if not math.isfinite(norm) or norm <= 1e-14:
        raise FloatingPointError(f"Degenerate {label}: norm={norm}")
    return vector / norm


def clamp(value: float) -> float:
    return max(-1.0, min(1.0, float(value)))


def acute_angle_deg(first: np.ndarray, second: np.ndarray) -> float:
    first = normalize(first, label="first acute-angle vector")
    second = normalize(second, label="second acute-angle vector")
    return math.degrees(math.acos(abs(clamp(float(first @ second)))))


def directed_angle_deg(first: np.ndarray, second: np.ndarray) -> float:
    first = normalize(first, label="first directed-angle vector")
    second = normalize(second, label="second directed-angle vector")
    if first.size != 2 or second.size != 2:
        raise ValueError("directed_angle_deg requires two-dimensional vectors")
    cross = float(first[0] * second[1] - first[1] * second[0])
    return math.degrees(math.atan2(cross, float(first @ second)))


def support_rhs_float(
    state: np.ndarray, *, kappa0: float, params: Any, audit: Any
) -> np.ndarray:
    """Central model derivative in support coordinate s=Theta/(2*pi)."""
    (
        y1,
        v1,
        theta1,
        omega1,
        y2,
        v2,
        theta2,
        omega2,
        z,
        vz,
        theta_global,
    ) = np.asarray(state, dtype=float)
    rho1 = math.exp(y1)
    rho2 = math.exp(y2)
    x = float(audit.sigmoid(z))
    g = x * (1.0 - x)
    delta = theta1 + theta2 - params.phi0 - theta_global
    sin_delta, cos_delta = math.sin(delta), math.cos(delta)
    kappa_effective = kappa0 * x
    rho1_dot = rho1 * v1
    rho2_dot = rho2 * v2
    rho1_ddot = (
        rho1 * omega1**2
        - audit.dU(rho1, params.e_reserve)
        - kappa_effective * rho2 * cos_delta
        - params.gamma_rho * rho1_dot
    )
    rho2_ddot = (
        rho2 * omega2**2
        - audit.dU(rho2, params.e_reserve)
        - kappa_effective * rho1 * cos_delta
        - params.gamma_rho * rho2_dot
    )
    v1_dot = rho1_ddot / rho1 - v1**2
    v2_dot = rho2_ddot / rho2 - v2**2
    omega1_dot = (
        kappa_effective * rho2 * sin_delta / rho1 - 2.0 * v1 * omega1
    )
    omega2_dot = (
        kappa_effective * rho1 * sin_delta / rho2 - 2.0 * v2 * omega2
    )
    activity = rho1 * rho2 * sin_delta
    x_dot = g * vz
    x_ddot = (
        -audit.dV_barrier(x, params.barrier_B)
        + params.mu_x * activity
        - params.gamma_x * x_dot
    ) / params.m_x
    vz_dot = x_ddot / g - (1.0 - 2.0 * x) * vz**2
    time_derivative = np.asarray(
        [
            v1,
            v1_dot,
            omega1,
            omega1_dot,
            v2,
            v2_dot,
            omega2,
            omega2_dot,
            vz,
            vz_dot,
            x,
        ],
        dtype=float,
    )
    return (2.0 * math.pi / x) * time_derivative


def read_selected_steps(
    path: Path, *, cell: int, requested_half_steps: int, stride: int
) -> tuple[list[dict[str, Any]], int]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    available = sorted(
        {int(row["half_steps"]) for row in rows if int(row["cell_index"]) == cell}
    )
    if not available:
        raise ValueError(f"No coupled-TM rows for cell {cell} in {path}")
    chosen = requested_half_steps if requested_half_steps > 0 else available[-1]
    if chosen not in available:
        raise ValueError(f"half_steps={chosen} unavailable; choices={available}")
    selected = [
        row
        for row in rows
        if int(row["cell_index"]) == cell and int(row["half_steps"]) == chosen
    ]
    selected.sort(
        key=lambda row: 0.5 * (float(row["physical_left"]) + float(row["physical_right"]))
    )
    expected = 2 * chosen
    if len(selected) != expected:
        raise ValueError(f"Expected {expected} rows, found {len(selected)}")
    selected = selected[::stride]
    for row in selected:
        row["sample_s"] = 0.5 * (
            float(row["physical_left"]) + float(row["physical_right"])
        )
    return selected, chosen


def fit_orbital_plane(points: np.ndarray, vertex: np.ndarray) -> dict[str, Any]:
    offsets = np.asarray(points, dtype=float) - np.asarray(vertex, dtype=float)[None, :]
    _, singular_values, vh = np.linalg.svd(offsets, full_matrices=False)
    if singular_values.size < 3 or singular_values[1] <= 1e-14:
        raise FloatingPointError("Bulge does not define a stable rank-two plane")
    plane = np.asarray(vh[:2, :], dtype=float).T
    energy = singular_values**2
    explained = float(np.sum(energy[:2]) / max(float(np.sum(energy)), 1e-300))
    residual_ratio = float(singular_values[2] / singular_values[1])
    return {
        "plane": plane,
        "singular_values": singular_values,
        "explained_fraction": explained,
        "residual_ratio": residual_ratio,
    }


def plane_principal_angles_deg(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    singular_values = np.linalg.svd(first.T @ second, compute_uv=False)
    singular_values = np.clip(singular_values, 0.0, 1.0)
    return np.sort(np.degrees(np.arccos(singular_values)))


def chord_profile(points: np.ndarray, phases: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    chord = (
        (1.0 - phases)[:, None] * points[0][None, :]
        + phases[:, None] * points[-1][None, :]
    )
    residual = points - chord
    amplitude = np.linalg.norm(residual, axis=1)
    return amplitude, residual


def best_profile_lag(
    phase: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
    maximum_lag: float,
    lag_samples: int,
) -> tuple[float, float]:
    left_scale = max(float(np.max(left)), 1e-300)
    right_scale = max(float(np.max(right)), 1e-300)
    left = np.asarray(left, dtype=float) / left_scale
    right = np.asarray(right, dtype=float) / right_scale
    best_lag = 0.0
    best_corr = -math.inf
    for lag in np.linspace(-maximum_lag, maximum_lag, lag_samples):
        shifted = phase + lag
        mask = (shifted >= 0.0) & (shifted <= 1.0)
        if np.count_nonzero(mask) < max(5, phase.size // 2):
            continue
        a = left[mask]
        b = np.interp(shifted[mask], phase, right)
        a = a - np.mean(a)
        b = b - np.mean(b)
        denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
        if denominator <= 1e-15:
            continue
        corr = float(a @ b) / denominator
        if corr > best_corr:
            best_corr = corr
            best_lag = float(lag)
    if not math.isfinite(best_corr):
        raise FloatingPointError("Could not determine bulge-profile phase lag")
    return best_lag, best_corr


def summarize(values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    return {
        "minimum": float(np.min(values)),
        "median": float(np.median(values)),
        "mean": float(np.mean(values)),
        "maximum": float(np.max(values)),
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit-script", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--orbit", type=Path, default=DEFAULT_ORBIT)
    parser.add_argument("--dense-raw", type=Path, default=DEFAULT_DENSE)
    parser.add_argument("--steps-csv", type=Path, default=DEFAULT_STEPS)
    parser.add_argument("--q", type=float, default=550.0)
    parser.add_argument("--cell", type=int, default=49)
    parser.add_argument("--half-steps", type=int, default=0)
    parser.add_argument("--sample-stride", type=int, default=4)
    parser.add_argument("--arc-samples", type=int, default=65)
    parser.add_argument("--maximum-profile-lag", type=float, default=0.25)
    parser.add_argument("--profile-lag-samples", type=int, default=201)
    parser.add_argument("--gamma-rho", type=float, default=0.3)
    parser.add_argument("--tail-periods", type=int, default=70)
    parser.add_argument("--steps-per-period", type=int, default=160)
    parser.add_argument("--samples-per-rotation", type=int, default=720)
    parser.add_argument("--initial-x-perturbation", type=float, default=1e-6)
    parser.add_argument(
        "--minimum-plane-separation-deg",
        type=float,
        default=1e-6,
        help="absolute numerical floor for resolved non-coplanarity",
    )
    parser.add_argument(
        "--plane-resolution-safety-factor",
        type=float,
        default=8.0,
        help="required separation/error ratio under arc-grid coarsening",
    )
    parser.add_argument("--maximum-axis-angle-deg", type=float, default=15.0)
    parser.add_argument("--output-prefix", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    if args.sample_stride < 1:
        parser.error("--sample-stride must be positive")
    if args.arc_samples < 9 or args.arc_samples % 2 == 0:
        parser.error("--arc-samples must be an odd integer >= 9")
    if not 0.0 <= args.maximum_profile_lag < 0.5:
        parser.error("--maximum-profile-lag must lie in [0,0.5)")
    if args.profile_lag_samples < 3:
        parser.error("--profile-lag-samples must be >= 3")
    if args.minimum_plane_separation_deg < 0.0:
        parser.error("--minimum-plane-separation-deg must be nonnegative")
    if args.plane_resolution_safety_factor <= 1.0:
        parser.error("--plane-resolution-safety-factor must exceed one")

    audit = load_module(args.audit_script, "bulge_cone_audit")
    selected, chosen_half_steps = read_selected_steps(
        args.steps_csv,
        cell=args.cell,
        requested_half_steps=args.half_steps,
        stride=args.sample_stride,
    )

    with np.load(args.dense_raw, allow_pickle=True) as dense:
        starts = np.asarray(dense["starts"], dtype=float)
        basis = np.asarray(dense["basis"], dtype=float)
        labels = [str(item) for item in dense["basis_labels"]]
        p_line = float(dense["p_line"])
    if basis.shape != (9, 9):
        raise ValueError(f"Expected 9x9 dense basis, got {basis.shape}")
    if labels[:2] != ["lower angular bascule", "upper counterbalance"]:
        raise ValueError(f"Unexpected active labels: {labels[:2]}")
    if args.cell < 0 or args.cell + 1 >= starts.size:
        raise ValueError(f"cell {args.cell} outside dense grid")

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
    extraction = np.asarray(audit.reduced_extraction(), dtype=float)
    basis_inverse = np.linalg.inv(basis)
    basis_active = basis[:, :2]
    kappa0 = (55.0 / 6.0) * (10.0 / args.q)
    phase_grid = np.linspace(0.0, 1.0, args.arc_samples)

    def reduced_state(s_values: np.ndarray) -> np.ndarray:
        state = np.asarray(spline(np.asarray(s_values, dtype=float)), dtype=float)
        if state.ndim == 1:
            state = state[:, None]
        if state.shape[-1] == 11:
            return state @ extraction.T
        if state.shape[0] == 11:
            return (extraction @ state).T
        else:
            raise ValueError(f"Unexpected spline state shape: {state.shape}")

    def reduced_tangent(s_value: float) -> np.ndarray:
        state = np.asarray(spline(float(s_value)), dtype=float)
        derivative = np.asarray(
            support_rhs_float(state, kappa0=kappa0, params=params, audit=audit),
            dtype=float,
        )
        return extraction @ derivative

    output_rows: list[dict[str, Any]] = []
    previous_pair_bisector: np.ndarray | None = None

    for row in selected:
        s_input = float(row["sample_s"])
        s_center = s_input + p_line
        s_output = s_input + 2.0 * p_line
        left_s = s_input + phase_grid * p_line
        right_s = s_center + phase_grid * p_line
        left_points = reduced_state(left_s)
        right_points = reduced_state(right_s)
        center_state = 0.5 * (left_points[-1] + right_points[0])

        left_fit = fit_orbital_plane(left_points, center_state)
        right_fit = fit_orbital_plane(right_points, center_state)
        plane_angles = plane_principal_angles_deg(
            left_fit["plane"], right_fit["plane"]
        )
        coarse_left_fit = fit_orbital_plane(left_points[::2], center_state)
        coarse_right_fit = fit_orbital_plane(right_points[::2], center_state)
        coarse_plane_angles = plane_principal_angles_deg(
            coarse_left_fit["plane"], coarse_right_fit["plane"]
        )
        plane_resolution_delta = np.abs(plane_angles - coarse_plane_angles)

        left_profile, _ = chord_profile(left_points, phase_grid)
        right_profile, _ = chord_profile(right_points, phase_grid)
        left_peak_index = int(np.argmax(left_profile))
        right_peak_index = int(np.argmax(right_profile))

        # Both rays are oriented with the temporal advance: the left ray points
        # from its bulge maximum to the common center, the right ray from the
        # common center to its bulge maximum.
        left_ray = normalize(
            center_state - left_points[left_peak_index], label="left cone ray"
        )
        right_ray = normalize(
            right_points[right_peak_index] - center_state, label="right cone ray"
        )
        tangent = normalize(reduced_tangent(s_center), label="center tangent")
        if float(left_ray @ tangent) < 0.0:
            left_ray = -left_ray
        if float(right_ray @ tangent) < 0.0:
            right_ray = -right_ray
        if float(left_ray @ right_ray) < 0.0:
            right_ray = -right_ray
        cone_axis = normalize(left_ray + right_ray, label="oriented cone axis")
        if float(cone_axis @ tangent) < 0.0:
            cone_axis = -cone_axis

        cone_opening = acute_angle_deg(left_ray, right_ray)
        left_half_angle = acute_angle_deg(cone_axis, left_ray)
        right_half_angle = acute_angle_deg(cone_axis, right_ray)
        axis_tangent_angle = acute_angle_deg(cone_axis, tangent)

        left_coeff = basis_inverse @ left_ray
        right_coeff = basis_inverse @ right_ray
        axis_coeff = basis_inverse @ cone_axis
        tangent_coeff = basis_inverse @ tangent
        left_active = normalize(left_coeff[:2], label="left active cone ray")
        right_active = normalize(right_coeff[:2], label="right active cone ray")
        axis_active = normalize(axis_coeff[:2], label="active cone axis")
        tangent_active = normalize(tangent_coeff[:2], label="active center tangent")
        signed_opening = directed_angle_deg(left_active, right_active)
        active_axis_tangent_angle = acute_angle_deg(axis_active, tangent_active)

        u = normalize(
            np.asarray([float(row["witness_u0"]), float(row["witness_u1"])]),
            label="output SVD axis",
        )
        v = normalize(
            np.asarray([float(row["witness_v0"]), float(row["witness_v1"])]),
            label="input SVD axis",
        )
        if float(u @ v) < 0.0:
            u = -u
        pair_bisector = normalize(u + v, label="SVD pair bisector")
        if previous_pair_bisector is not None and float(
            pair_bisector @ previous_pair_bisector
        ) < 0.0:
            pair_bisector = -pair_bisector
        if previous_pair_bisector is None and float(pair_bisector @ axis_active) < 0.0:
            pair_bisector = -pair_bisector
        previous_pair_bisector = pair_bisector.copy()
        pair_cone_angle = acute_angle_deg(pair_bisector, axis_active)

        best_lag_fraction, best_profile_correlation = best_profile_lag(
            phase_grid,
            left_profile,
            right_profile,
            args.maximum_profile_lag,
            args.profile_lag_samples,
        )
        peak_phase_lag = float(
            phase_grid[right_peak_index] - phase_grid[left_peak_index]
        )

        output_rows.append(
            {
                "cell_index": args.cell,
                "half_steps": chosen_half_steps,
                "side": row["side"],
                "step_index": int(row["step_index"]),
                "s_input": s_input,
                "s_center": s_center,
                "s_output": s_output,
                "left_plane_principal_angle_deg": float(plane_angles[0]),
                "right_plane_principal_angle_deg": float(plane_angles[1]),
                "plane_angle1_resolution_delta_deg": float(
                    plane_resolution_delta[0]
                ),
                "plane_angle2_resolution_delta_deg": float(
                    plane_resolution_delta[1]
                ),
                "plane_twist_rate_deg_per_support": float(
                    plane_angles[1] / p_line
                ),
                "plane_chordal_distance": float(
                    np.linalg.norm(np.sin(np.radians(plane_angles)))
                ),
                "left_plane_explained_fraction": float(
                    left_fit["explained_fraction"]
                ),
                "right_plane_explained_fraction": float(
                    right_fit["explained_fraction"]
                ),
                "left_plane_rank3_over_rank2": float(left_fit["residual_ratio"]),
                "right_plane_rank3_over_rank2": float(right_fit["residual_ratio"]),
                "cone_opening_deg": cone_opening,
                "cone_left_half_angle_deg": left_half_angle,
                "cone_right_half_angle_deg": right_half_angle,
                "cone_axis_tangent_angle_deg": axis_tangent_angle,
                "active_cone_axis_tangent_angle_deg": active_axis_tangent_angle,
                "active_signed_cone_opening_deg": signed_opening,
                "left_bulge_peak_phase": float(phase_grid[left_peak_index]),
                "right_bulge_peak_phase": float(phase_grid[right_peak_index]),
                "peak_phase_lag_fraction": peak_phase_lag,
                "peak_phase_lag_support": peak_phase_lag * p_line,
                "profile_best_lag_fraction": best_lag_fraction,
                "profile_best_lag_support": best_lag_fraction * p_line,
                "profile_best_correlation": best_profile_correlation,
                "left_bulge_amplitude": float(left_profile[left_peak_index]),
                "right_bulge_amplitude": float(right_profile[right_peak_index]),
                "bulge_amplitude_ratio_right_over_left": float(
                    right_profile[right_peak_index]
                    / max(left_profile[left_peak_index], 1e-300)
                ),
                "cone_axis_active_fraction": float(
                    np.linalg.norm(axis_coeff[:2])
                    / max(float(np.linalg.norm(axis_coeff)), 1e-300)
                ),
                "center_tangent_active_fraction": float(
                    np.linalg.norm(tangent_coeff[:2])
                    / max(float(np.linalg.norm(tangent_coeff)), 1e-300)
                ),
                "svd_pair_to_geometric_cone_axis_angle_deg": pair_cone_angle,
                "witness_singular_gap": float(row["witness_singular_gap"]),
            }
        )

    arrays = {
        key: np.asarray([float(row[key]) for row in output_rows], dtype=float)
        for key in (
            "left_plane_principal_angle_deg",
            "right_plane_principal_angle_deg",
            "plane_angle1_resolution_delta_deg",
            "plane_angle2_resolution_delta_deg",
            "plane_twist_rate_deg_per_support",
            "plane_chordal_distance",
            "left_plane_explained_fraction",
            "right_plane_explained_fraction",
            "cone_opening_deg",
            "cone_axis_tangent_angle_deg",
            "active_cone_axis_tangent_angle_deg",
            "active_signed_cone_opening_deg",
            "peak_phase_lag_fraction",
            "profile_best_lag_fraction",
            "profile_best_correlation",
            "bulge_amplitude_ratio_right_over_left",
            "svd_pair_to_geometric_cone_axis_angle_deg",
        )
    }

    resolved_plane_threshold = max(
        args.minimum_plane_separation_deg,
        args.plane_resolution_safety_factor
        * float(np.max(arrays["plane_angle2_resolution_delta_deg"])),
    )
    noncoplanar = bool(
        float(np.min(arrays["right_plane_principal_angle_deg"]))
        > resolved_plane_threshold
    )
    oriented = bool(
        float(np.max(arrays["active_cone_axis_tangent_angle_deg"]))
        <= args.maximum_axis_angle_deg
    )
    phase_signed = bool(
        np.all(arrays["profile_best_lag_fraction"] > 0.0)
        or np.all(arrays["profile_best_lag_fraction"] < 0.0)
    )
    if noncoplanar and oriented:
        classification = (
            "TWO NON-COPLANAR COUPLED BULGES OPEN AN ORIENTED ADVANCING CONE: "
            "NUMERICALLY SUPPORTED"
        )
    else:
        classification = "NON-COPLANAR ORIENTED BULGE CONE: NOT ESTABLISHED"

    summary = {
        "script_version": SCRIPT_VERSION,
        "classification": classification,
        "status": "floating-point geometric diagnostic",
        "q": args.q,
        "cell_index": args.cell,
        "cell_interval": [float(starts[args.cell]), float(starts[args.cell + 1])],
        "p_line": p_line,
        "sample_count": len(output_rows),
        "selected_half_steps": chosen_half_steps,
        "sample_stride": args.sample_stride,
        "arc_samples": args.arc_samples,
        "active_labels": labels[:2],
        "geometry": {
            "left_bulge": "orbit arc [s,s+P_line] fitted in reduced 9D",
            "right_bulge": "orbit arc [s+P_line,s+2P_line] fitted in reduced 9D",
            "shared_vertex": "q(s+P_line)",
            "cone_rays": "temporally forward rays through maximum chord-deviation points",
        },
        "minimum_plane_separation_deg": args.minimum_plane_separation_deg,
        "plane_resolution_safety_factor": args.plane_resolution_safety_factor,
        "resolved_plane_threshold_deg": resolved_plane_threshold,
        "maximum_axis_angle_deg": args.maximum_axis_angle_deg,
        "noncoplanar_gate": noncoplanar,
        "advance_orientation_gate": oriented,
        "profile_lag_has_uniform_sign": phase_signed,
        "metrics": {key: summarize(value) for key, value in arrays.items()},
        "sources": {
            "orbit": str(args.orbit),
            "orbit_sha256": sha256(args.orbit),
            "dense_raw": str(args.dense_raw),
            "dense_raw_sha256": sha256(args.dense_raw),
            "coupled_steps": str(args.steps_csv),
            "coupled_steps_sha256": sha256(args.steps_csv),
            "audit_script": str(args.audit_script),
            "audit_script_sha256": sha256(args.audit_script),
        },
        "proof_boundary": (
            "central floating-point orbit only; no interval plane, cone, phase, "
            "return, nesting, or hysteresis certificate"
        ),
    }

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = Path(str(prefix) + "_samples.csv")
    json_path = Path(str(prefix) + ".json")
    report_path = Path(str(prefix) + "_report.txt")
    png_path = Path(str(prefix) + ".png")
    write_csv(csv_path, output_rows)
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    lines = [
        f"Script version                         = {SCRIPT_VERSION}",
        f"classification                         = {classification}",
        f"cell                                   = {args.cell}",
        f"samples                                = {len(output_rows)}",
        f"P_line                                 = {p_line:.15f}",
        "",
        "Two orbital planes in reduced 9D",
        f"  principal angle 1 range              = [{np.min(arrays['left_plane_principal_angle_deg']):.9f}, {np.max(arrays['left_plane_principal_angle_deg']):.9f}] deg",
        f"  principal angle 2 range              = [{np.min(arrays['right_plane_principal_angle_deg']):.9f}, {np.max(arrays['right_plane_principal_angle_deg']):.9f}] deg",
        f"  angle 2 resolution delta max         = {np.max(arrays['plane_angle2_resolution_delta_deg']):.12e} deg",
        f"  resolved separation threshold        = {resolved_plane_threshold:.12e} deg",
        f"  plane twist rate range               = [{np.min(arrays['plane_twist_rate_deg_per_support']):.9f}, {np.max(arrays['plane_twist_rate_deg_per_support']):.9f}] deg/support",
        f"  left plane explained fraction min    = {np.min(arrays['left_plane_explained_fraction']):.12e}",
        f"  right plane explained fraction min   = {np.min(arrays['right_plane_explained_fraction']):.12e}",
        "",
        "Oriented cone",
        f"  cone opening range                   = [{np.min(arrays['cone_opening_deg']):.9f}, {np.max(arrays['cone_opening_deg']):.9f}] deg",
        f"  active axis / advance angle max      = {np.max(arrays['active_cone_axis_tangent_angle_deg']):.9f} deg",
        f"  full 9D axis / tangent angle range   = [{np.min(arrays['cone_axis_tangent_angle_deg']):.9f}, {np.max(arrays['cone_axis_tangent_angle_deg']):.9f}] deg",
        f"  SVD pair / geometric cone angle max  = {np.max(arrays['svd_pair_to_geometric_cone_axis_angle_deg']):.9f} deg",
        "",
        "Relative phase",
        f"  peak lag fraction range              = [{np.min(arrays['peak_phase_lag_fraction']):+.9f}, {np.max(arrays['peak_phase_lag_fraction']):+.9f}]",
        f"  profile lag fraction range           = [{np.min(arrays['profile_best_lag_fraction']):+.9f}, {np.max(arrays['profile_best_lag_fraction']):+.9f}]",
        f"  profile correlation min              = {np.min(arrays['profile_best_correlation']):.12e}",
        f"  right/left amplitude ratio range     = [{np.min(arrays['bulge_amplitude_ratio_right_over_left']):.9f}, {np.max(arrays['bulge_amplitude_ratio_right_over_left']):.9f}]",
        "",
        f"noncoplanar gate                       = {noncoplanar}",
        f"advance orientation gate               = {oriented}",
        f"profile lag has uniform sign           = {phase_signed}",
        "",
        "Boundary: floating-point orbit geometry; no interval cone or hysteresis theorem.",
    ]
    report_path.write_text("\n".join(lines) + "\n")

    centers = np.asarray([float(row["s_center"]) for row in output_rows])
    fig, axes = plt.subplots(4, 1, figsize=(11, 12), sharex=True, constrained_layout=True)
    axes[0].plot(centers, arrays["left_plane_principal_angle_deg"], label="plane principal angle 1")
    axes[0].plot(centers, arrays["right_plane_principal_angle_deg"], label="plane principal angle 2")
    axes[0].axhline(resolved_plane_threshold, color="black", linestyle="--", linewidth=1, label="resolved non-coplanar gate")
    axes[0].set_ylabel("plane angle [deg]")
    axes[0].legend(loc="best")

    axes[1].plot(centers, arrays["cone_opening_deg"], label="cone opening")
    axes[1].plot(centers, arrays["active_cone_axis_tangent_angle_deg"], label="axis / active advance")
    axes[1].plot(centers, arrays["svd_pair_to_geometric_cone_axis_angle_deg"], label="SVD pair / cone axis")
    axes[1].axhline(args.maximum_axis_angle_deg, color="black", linestyle="--", linewidth=1, label="orientation gate")
    axes[1].set_ylabel("angle [deg]")
    axes[1].legend(loc="best")

    axes[2].plot(centers, arrays["peak_phase_lag_fraction"], label="peak lag")
    axes[2].plot(centers, arrays["profile_best_lag_fraction"], label="profile-correlation lag")
    axes[2].axhline(0.0, color="black", linewidth=1)
    axes[2].set_ylabel("lag / P_line")
    axes[2].legend(loc="best")

    axes[3].plot(centers, arrays["bulge_amplitude_ratio_right_over_left"], label="right / left bulge amplitude")
    axes[3].plot(centers, arrays["profile_best_correlation"], label="profile correlation")
    axes[3].axhline(1.0, color="black", linewidth=1)
    axes[3].set_ylabel("ratio / correlation")
    axes[3].set_xlabel("shared-center support coordinate s + P_line")
    axes[3].legend(loc="best")
    fig.suptitle("Q=550 cell 49: two non-coplanar coupled bulges and oriented cone")
    fig.savefig(png_path, dpi=180)
    plt.close(fig)

    print("\n".join(lines))
    print(f"samples CSV                            = {csv_path}")
    print(f"JSON ledger                            = {json_path}")
    print(f"plot                                   = {png_path}")
    print(f"report                                 = {report_path}")


if __name__ == "__main__":
    main()
