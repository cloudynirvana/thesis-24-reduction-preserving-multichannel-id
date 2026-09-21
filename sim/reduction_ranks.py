#!/usr/bin/env python3
"""Full versus reduced Fisher ranks for one metabolic ODE.

The generator extends the four-state steady map used in Thesis #9 by one
fast modifier, then removes that modifier by the quasi-steady product used
in Thesis #12. Channel schedules are applied to both parameter vectors.
Ranks that clear the practical cut on the full vector and miss it on the
reduced vector are flagged as artefacts of the unreduced coordinates.

Noise-free Gaussian information of the mean. One seeded draw is stored only
as a noisy check. Synthetic inputs. Not a CCLE file. Not a medical device.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import cholesky, solve_triangular

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

SEED = 20260921
RNG = np.random.default_rng(SEED)

# Log-parameter order for the unreduced field.
FULL_NAMES = ["k_gp", "k_pl", "k_pq", "d_l", "d_q", "a", "b", "h"]
RED_NAMES = ["k_gp", "k_pl", "k_pq", "d_l", "d_q", "kappa"]
THETA0 = np.array([0.80, 0.50, 0.30, 0.40, 0.60, 1.60, 8.00, 0.45], dtype=float)

# Synthetic input labels. The numbers are not CCLE or DepMap measurements.
LINEAGES = ["BRCA", "LUAD", "COAD", "AML"]
U_G = {"BRCA": 1.00, "LUAD": 1.40, "COAD": 0.80, "AML": 1.20}
U_Q = {"BRCA": 0.40, "LUAD": 0.55, "COAD": 0.50, "AML": 0.30}

N_REP = {"S": 6, "M": 6, "Z": 6, "D": 4}
SIGMA = {"S": 0.08, "M": 0.15, "Z": 0.15, "D": 0.10}
# Declared before the spectra were read. Same fractions as Thesis #9 and #12.
NUM_CUT = 1e-8
PRACTICAL_CUT = 1e-3
SUPPORT_MASS = 0.50
FD_STEADY = 1e-6
FD_TRANSIENT = 1e-4

# Correlated log-concentration noise on (G, P, L, Q). Not estimated from a panel.
CORR = np.array(
    [
        [1.00, 0.55, 0.45, 0.15],
        [0.55, 1.00, 0.60, 0.20],
        [0.45, 0.60, 1.00, 0.10],
        [0.15, 0.20, 0.10, 1.00],
    ],
    dtype=float,
)

Y_GPLQ = np.array([1.20, 0.40, 0.30, 0.50], dtype=float)
Z_OFF = 0.02
T_EARLY = np.array([0.05, 0.15, 0.40, 1.0, 3.0, 8.0, 15.0])
T_LATE = np.array([1.0, 2.0, 4.0, 8.0, 15.0])
# Spends the steady-null ratio rho = a/b as an extra clearance of Q.
BAD_GAMMA = 0.50
CHANNELS4 = ["G", "P", "L", "Q"]


def unpack_full(th: np.ndarray) -> tuple[float, ...]:
    return tuple(float(v) for v in th)


def kappa_of(th: np.ndarray) -> float:
    _k_gp, _k_pl, _k_pq, _d_l, _d_q, a, b, h = unpack_full(th)
    return h * a / b


def rho_of(th: np.ndarray) -> float:
    _k_gp, _k_pl, _k_pq, _d_l, _d_q, a, b, _h = unpack_full(th)
    return a / b


def phi_of(th: np.ndarray) -> np.ndarray:
    k_gp, k_pl, k_pq, d_l, d_q, a, b, h = unpack_full(th)
    return np.array([k_gp, k_pl, k_pq, d_l, d_q, h * a / b], dtype=float)


def steady_reduced(phi: np.ndarray, u_g: float, u_q: float) -> np.ndarray:
    """Algebraic steady state of (G, P, L, Q). Positive root of the quadratic in P."""
    k_gp, k_pl, k_pq, d_l, d_q, kappa = (float(v) for v in phi)
    if min(k_gp, k_pl, k_pq, d_l, d_q) <= 0.0 or kappa <= 0.0:
        return np.full(4, np.nan)
    a_coef = k_pl * kappa
    b_coef = k_pl + k_pq
    disc = b_coef * b_coef + 4.0 * a_coef * u_g
    if disc <= 0.0:
        return np.full(4, np.nan)
    p = (-b_coef + np.sqrt(disc)) / (2.0 * a_coef)
    g = u_g / k_gp
    lac = k_pl * (1.0 + kappa * p) * p / d_l
    q = (u_q + k_pq * p) / d_q
    out = np.array([g, p, lac, q], dtype=float)
    if np.any(out <= 0.0) or not np.all(np.isfinite(out)):
        return np.full(4, np.nan)
    return out


def steady_full(th: np.ndarray, u_g: float, u_q: float) -> np.ndarray:
    """Steady (G, P, L, Q, Z). Z is slaved, so the first four match the reduced map."""
    base = steady_reduced(phi_of(th), u_g, u_q)
    if not np.all(np.isfinite(base)):
        return np.full(5, np.nan)
    z = rho_of(th) * base[1]
    return np.array([base[0], base[1], base[2], base[3], z], dtype=float)


def steady_bad(th: np.ndarray, u_g: float, u_q: float) -> np.ndarray:
    """Same steady map, except rho is spent as extra clearance of Q."""
    base = steady_full(th, u_g, u_q)
    if not np.all(np.isfinite(base)):
        return base
    _k_gp, _k_pl, k_pq, _d_l, d_q, _a, _b, _h = unpack_full(th)
    p = base[1]
    q = (u_q + k_pq * p) / (d_q + BAD_GAMMA * rho_of(th))
    out = base.copy()
    out[3] = q
    return out


def rhs_full(t: float, y: np.ndarray, th: np.ndarray, u_g: float, u_q: float) -> np.ndarray:
    g, p, lac, q, z = y
    k_gp, k_pl, k_pq, d_l, d_q, a, b, h = unpack_full(th)
    branch = k_pl * (1.0 + h * z)
    return np.array(
        [
            u_g - k_gp * g,
            k_gp * g - (branch + k_pq) * p,
            branch * p - d_l * lac,
            u_q + k_pq * p - d_q * q,
            a * p - b * z,
        ],
        dtype=float,
    )


def rhs_reduced(t: float, y: np.ndarray, phi: np.ndarray, u_g: float, u_q: float) -> np.ndarray:
    g, p, lac, q = y
    k_gp, k_pl, k_pq, d_l, d_q, kappa = (float(v) for v in phi)
    flux = k_pl * (1.0 + kappa * p) * p
    return np.array(
        [
            u_g - k_gp * g,
            k_gp * g - flux - k_pq * p,
            flux - d_l * lac,
            u_q + k_pq * p - d_q * q,
        ],
        dtype=float,
    )


def _integrate(fun, y0: np.ndarray, arg: np.ndarray, u_g: float, u_q: float, t_eval: np.ndarray) -> np.ndarray:
    sol = solve_ivp(
        fun,
        (0.0, float(t_eval[-1])),
        y0,
        t_eval=t_eval,
        args=(arg, u_g, u_q),
        method="LSODA",
        rtol=1e-8,
        atol=1e-9,
        max_step=0.05,
    )
    if (not sol.success) or sol.y.shape[1] != t_eval.size or not np.all(np.isfinite(sol.y)):
        return np.full((y0.size, t_eval.size), np.nan)
    return sol.y


def whiten4() -> np.ndarray:
    cov = (SIGMA["M"] ** 2) * CORR
    chol = cholesky(cov, lower=True)
    return solve_triangular(chol, np.eye(4), lower=True)


W4 = whiten4()


def _log_safe(x: np.ndarray) -> np.ndarray:
    return np.log(np.maximum(x, 1e-12))


def predict(th_or_phi: np.ndarray, schedule: str, kind: str) -> np.ndarray:
    """Noise-free mean in the coordinates the Fisher sees.

    kind is 'full', 'reduced', or 'bad'. Schedule Z is refused on the reduced field.
    """
    vec = np.asarray(th_or_phi, dtype=float)
    blocks = []
    for name in LINEAGES:
        u_g, u_q = U_G[name], U_Q[name]
        if schedule in ("S", "M", "Z", "Mbeta"):
            if kind == "full":
                st = steady_full(vec, u_g, u_q)
            elif kind == "reduced":
                st = steady_reduced(vec, u_g, u_q)
            elif kind == "bad":
                st = steady_bad(vec, u_g, u_q)
            else:
                raise ValueError(kind)
            if schedule == "S":
                blocks.append([np.log(st[2] / st[0])])
            elif schedule == "M":
                blocks.append(_log_safe(st[:4]))
            elif schedule == "Mbeta":
                blocks.append(_log_safe(st[:4]))
            elif schedule == "Z":
                if kind != "full":
                    raise RuntimeError("Z is not a state of the reduced field")
                blocks.append(_log_safe(st))
            continue
        if schedule in ("Dearly", "Dlate"):
            times = T_EARLY if schedule == "Dearly" else T_LATE
            if kind == "full":
                z0 = Z_OFF if schedule == "Dearly" else rho_of(THETA0) * Y_GPLQ[1]
                # Dlate always starts off the manifold; the layer has died by t = 1.
                if schedule == "Dlate":
                    z0 = Z_OFF
                y0 = np.concatenate([Y_GPLQ, [z0]])
                traj = _integrate(rhs_full, y0, vec, u_g, u_q, times)
                obs = _log_safe(traj[:4, :]).T.reshape(-1)
            elif kind == "reduced":
                traj = _integrate(rhs_reduced, Y_GPLQ, vec, u_g, u_q, times)
                obs = _log_safe(traj).T.reshape(-1)
            else:
                raise ValueError(kind)
            blocks.append(obs)
            continue
        raise ValueError(schedule)
    return np.concatenate(blocks)


def predict_early(th: np.ndarray, z0: float) -> np.ndarray:
    """Full-model early time course at a declared Z(0)."""
    blocks = []
    for name in LINEAGES:
        y0 = np.concatenate([Y_GPLQ, [z0]])
        traj = _integrate(rhs_full, y0, th, U_G[name], U_Q[name], T_EARLY)
        blocks.append(_log_safe(traj[:4, :]).T.reshape(-1))
    return np.concatenate(blocks)


def _weight_matrix(schedule: str, n_obs_one: int) -> np.ndarray:
    """Return the matrix R such that the Fisher is (R J)^T (R J), including replicates."""
    n_lin = len(LINEAGES)
    if schedule == "S":
        scale = np.sqrt(N_REP["S"]) / SIGMA["S"]
        return scale * np.eye(n_lin)
    if schedule == "M":
        # Block-diagonal whitening, then replicate weight.
        blocks = [np.sqrt(N_REP["M"]) * W4 for _ in range(n_lin)]
        return _blkdiag(blocks)
    if schedule == "Mbeta":
        blocks = [np.sqrt(N_REP["M"]) * W4 for _ in range(n_lin)]
        return _blkdiag(blocks)
    if schedule == "Z":
        # Independent noise on five logs. The four-channel correlation is not copied onto Z.
        scale = np.sqrt(N_REP["Z"]) / SIGMA["Z"]
        return scale * np.eye(n_lin * 5)
    if schedule in ("Dearly", "Dlate"):
        scale = np.sqrt(N_REP["D"]) / SIGMA["D"]
        return scale * np.eye(n_obs_one * n_lin)
    raise ValueError(schedule)


def _blkdiag(blocks: list[np.ndarray]) -> np.ndarray:
    n = sum(b.shape[0] for b in blocks)
    out = np.zeros((n, n), dtype=float)
    i = 0
    for b in blocks:
        k = b.shape[0]
        out[i : i + k, i : i + k] = b
        i += k
    return out


def _jacobian(predict_fun, x: np.ndarray, step: float) -> tuple[np.ndarray, np.ndarray]:
    mu = predict_fun(x)
    if not np.all(np.isfinite(mu)):
        raise RuntimeError("non-finite prediction at the evaluation point")
    sens = np.zeros((mu.size, x.size), dtype=float)
    for j in range(x.size):
        up = x.copy()
        dn = x.copy()
        up[j] = x[j] * np.exp(step)
        dn[j] = x[j] * np.exp(-step)
        sp = predict_fun(up)
        sm = predict_fun(dn)
        if not (np.all(np.isfinite(sp)) and np.all(np.isfinite(sm))):
            raise RuntimeError(f"non-finite perturbation of coordinate {j}")
        sens[:, j] = (sp - sm) / (2.0 * step)
    return mu, sens


def fisher_from_sensitivity(sens: np.ndarray, weight: np.ndarray) -> np.ndarray:
    js = weight @ sens
    fim = js.T @ js
    return 0.5 * (fim + fim.T)


def analytic_axes() -> dict[str, np.ndarray]:
    """Orthonormal directions in the (a, b, h) block. Kinetic entries are zero."""
    e_kappa = np.zeros(8)
    e_speed = np.zeros(8)
    e_rho = np.zeros(8)
    e_kappa[5:] = np.array([1.0, -1.0, 1.0]) / np.sqrt(3.0)
    e_speed[5:] = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
    e_rho[5:] = np.array([1.0, -1.0, -2.0]) / np.sqrt(6.0)
    return {"kappa": e_kappa, "speed": e_speed, "rho": e_rho}


def reduced_support_mass(evec: np.ndarray) -> float:
    """Squared mass on the five rates plus the kappa direction in (a, b, h)."""
    axes = analytic_axes()
    kinetic = float(np.sum(evec[:5] ** 2))
    along_kappa = float(np.dot(evec, axes["kappa"]) ** 2)
    return kinetic + along_kappa


def spectrum_report(fim: np.ndarray, names: list[str]) -> dict:
    evals, evecs = np.linalg.eigh(0.5 * (fim + fim.T))
    order = np.argsort(evals)[::-1]
    evals = np.clip(evals[order], 0.0, None)
    evecs = evecs[:, order]
    lam_max = float(evals[0]) if evals.size else 0.0
    ratios = evals / lam_max if lam_max > 0 else np.zeros_like(evals)
    numerical_rank = int(np.sum(ratios > NUM_CUT))
    practical_rank = int(np.sum(ratios > PRACTICAL_CUT))
    stiff = ratios > PRACTICAL_CUT
    supported = []
    masses = []
    for j in range(evecs.shape[1]):
        mass = reduced_support_mass(evecs[:, j]) if len(names) == 8 else 1.0
        masses.append(mass)
        if stiff[j] and mass >= SUPPORT_MASS:
            supported.append(j)
    retained = evals[stiff]
    if retained.size >= 2 and retained[-1] > 0:
        cond = float(retained[0] / retained[-1])
    elif retained.size == 1:
        cond = 1.0
    else:
        cond = None
    axes = analytic_axes()
    alignment = {}
    if len(names) == 8:
        for key, vec in axes.items():
            alignment[key] = [float(np.dot(evecs[:, j], vec) ** 2) for j in range(evecs.shape[1])]
    return {
        "names": names,
        "n_param": len(names),
        "eigenvalues_desc": [float(v) for v in evals],
        "log10_eigenvalues": [float(np.log10(v + 1e-30)) for v in evals],
        "ratios_to_max": [float(v) for v in ratios],
        "numerical_rank": numerical_rank,
        "practical_rank": practical_rank,
        "lambda_max": lam_max,
        "condition_of_stiff_block": cond,
        "support_mass": masses,
        "n_practical_supported_on_reduced_coordinates": len(supported) if len(names) == 8 else practical_rank,
        "alignment_sq_with_analytic_axes": alignment,
        "eigenvectors_desc_columns": evecs.tolist(),
    }


def chi_square(mu: np.ndarray, mu_ref: np.ndarray, schedule: str) -> float:
    weight = _weight_matrix(schedule, mu_ref.size // len(LINEAGES))
    r = weight @ (mu - mu_ref)
    return float(r @ r)


def fisher_pair(schedule: str, step: float) -> dict:
    if schedule == "Z":
        mu, sens = _jacobian(lambda th: predict(th, "Z", "full"), THETA0, step)
        weight = _weight_matrix("Z", mu.size // len(LINEAGES))
        fim = fisher_from_sensitivity(sens, weight)
        spec = spectrum_report(fim, FULL_NAMES)
        return {"full": spec, "reduced": None, "schedule": schedule}
    mu_f, sens_f = _jacobian(lambda th: predict(th, schedule, "full"), THETA0, step)
    phi = phi_of(THETA0)
    mu_r, sens_r = _jacobian(lambda p: predict(p, schedule, "reduced"), phi, step)
    w = _weight_matrix(schedule, mu_f.size // len(LINEAGES))
    spec_f = spectrum_report(fisher_from_sensitivity(sens_f, w), FULL_NAMES)
    spec_r = spectrum_report(fisher_from_sensitivity(sens_r, w), RED_NAMES)
    # Drop eigenvectors from the pair summary later; keep them for the flow only on M.
    return {
        "full": spec_f,
        "reduced": spec_r,
        "schedule": schedule,
        "steady_log_rmse": float(np.sqrt(np.mean((mu_f - mu_r) ** 2))) if schedule in ("S", "M", "Mbeta") else None,
    }


def profiled_kinetic_spectrum(kind: str) -> dict:
    """Four-channel snapshot plus one loading per lineage, then the Schur complement on the rates."""
    if kind == "full":
        x0 = THETA0
        n_kin = 8
        names = FULL_NAMES

        def pred(x: np.ndarray) -> np.ndarray:
            return predict(x, "M", "full")

    else:
        x0 = phi_of(THETA0)
        n_kin = 6
        names = RED_NAMES

        def pred(x: np.ndarray) -> np.ndarray:
            return predict(x, "M", "reduced")

    _mu, sens_k = _jacobian(pred, x0, FD_STEADY)
    n_beta = len(LINEAGES)
    n_obs = sens_k.shape[0]
    sens = np.zeros((n_obs, n_kin + n_beta), dtype=float)
    sens[:, :n_kin] = sens_k
    block = 4
    for i in range(n_beta):
        sens[i * block : (i + 1) * block, n_kin + i] = 1.0
    weight = _weight_matrix("M", block)
    fim = fisher_from_sensitivity(sens, weight)
    joint_names = names + [f"beta_{lin}" for lin in LINEAGES]
    joint = spectrum_report(fim, joint_names if kind == "full" else joint_names)
    # spectrum_report's support mass assumes length 8. Override for the joint matrix.
    joint["n_practical_supported_on_reduced_coordinates"] = None
    joint["support_mass"] = None
    joint["alignment_sq_with_analytic_axes"] = {}
    f_tt = fim[:n_kin, :n_kin]
    f_tb = fim[:n_kin, n_kin:]
    f_bb = fim[n_kin:, n_kin:]
    schur = f_tt - f_tb @ np.linalg.solve(f_bb, f_tb.T)
    profiled = spectrum_report(schur, names)
    return {"joint": _strip_evecs(joint), "profiled_kinetics": _strip_evecs(profiled)}


def _strip_evecs(spec: dict) -> dict:
    out = dict(spec)
    out.pop("eigenvectors_desc_columns", None)
    return out


def principal_angles(spec: dict) -> dict:
    """Angles between the practical sloppy subspace and the analytic (speed, rho) plane."""
    evals = np.array(spec["ratios_to_max"])
    evecs = np.array(spec["eigenvectors_desc_columns"])
    sloppy = evecs[:, evals <= PRACTICAL_CUT]
    axes = analytic_axes()
    plane = np.column_stack([axes["speed"], axes["rho"]])
    if sloppy.size == 0:
        return {"n_sloppy": 0}
    # Singular values of plane^T sloppy are cosines of principal angles.
    sv = np.linalg.svd(plane.T @ sloppy, compute_uv=False)
    sv = np.clip(sv, 0.0, 1.0)
    angles = np.degrees(np.arccos(sv))
    return {
        "n_sloppy": int(sloppy.shape[1]),
        "cosines": [float(v) for v in sv],
        "angles_degrees": [float(v) for v in angles],
    }


def null_flow(spec_m: dict) -> dict:
    """Walk the exact rho direction, which keeps kappa and the steady (G, P, L, Q) map."""
    axes = analytic_axes()
    direction = axes["rho"]
    mu_m = predict(THETA0, "M", "full")
    mu_z = predict(THETA0, "Z", "full")
    mu_s = predict(THETA0, "S", "full")
    rows = []
    for step in range(0, 13):
        # direction is a unit vector. Step length 0.12 in log-parameter space.
        th = THETA0 * np.exp((0.12 * step) * direction)
        row = {
            "step": step,
            "log_step": 0.12 * step,
            "kappa": kappa_of(th),
            "rho": rho_of(th),
            "b": float(th[6]),
            "chi_M": chi_square(predict(th, "M", "full"), mu_m, "M"),
            "chi_S": chi_square(predict(th, "S", "full"), mu_s, "S"),
            "chi_Z": chi_square(predict(th, "Z", "full"), mu_z, "Z"),
            "parameter_ratios": (th / THETA0).tolist(),
        }
        rows.append(row)
    # Numerical trailing eigenvector, one step, to compare with the exact plane.
    evecs = np.array(spec_m["eigenvectors_desc_columns"])
    trailing = evecs[:, -1]
    # Flip so the rho-component is positive.
    if np.dot(trailing, direction) < 0:
        trailing = -trailing
    return {
        "exact_rho_walk": rows,
        "trailing_squared_overlap_with_rho": float(np.dot(trailing, direction) ** 2),
        "trailing_squared_overlap_with_speed": float(np.dot(trailing, axes["speed"]) ** 2),
        "trailing_squared_overlap_with_kappa": float(np.dot(trailing, axes["kappa"]) ** 2),
        "trailing_dlog_kappa": float(np.dot(trailing, axes["kappa"]) * np.sqrt(3.0)),
    }


def bad_reduction_check() -> dict:
    mu_m = predict(THETA0, "M", "full")
    mu_s = predict(THETA0, "S", "full")
    mu_bad_m = predict(THETA0, "M", "bad")
    mu_bad_s = predict(THETA0, "S", "bad")
    q_rel = []
    for name in LINEAGES:
        good = steady_full(THETA0, U_G[name], U_Q[name])[3]
        bad = steady_bad(THETA0, U_G[name], U_Q[name])[3]
        q_rel.append(float((bad - good) / good))
    return {
        "gamma": BAD_GAMMA,
        "rho": rho_of(THETA0),
        "extra_clearance": BAD_GAMMA * rho_of(THETA0),
        "chi_S": chi_square(mu_bad_s, mu_s, "S"),
        "chi_M": chi_square(mu_bad_m, mu_m, "M"),
        "relative_Q_shift": {name: q_rel[i] for i, name in enumerate(LINEAGES)},
        "max_abs_relative_Q_shift": float(np.max(np.abs(q_rel))),
    }


def early_initial_layer() -> dict:
    """Same early observation times. Two hidden initial conditions on the full field."""
    z_on = rho_of(THETA0) * float(Y_GPLQ[1])
    mu_off, sens_off = _jacobian(lambda th: predict_early(th, Z_OFF), THETA0, FD_TRANSIENT)
    mu_on, sens_on = _jacobian(lambda th: predict_early(th, z_on), THETA0, FD_TRANSIENT)
    phi = phi_of(THETA0)
    _mu_r, sens_r = _jacobian(lambda p: predict(p, "Dearly", "reduced"), phi, FD_TRANSIENT)
    w = _weight_matrix("Dearly", mu_off.size // len(LINEAGES))
    spec_off = spectrum_report(fisher_from_sensitivity(sens_off, w), FULL_NAMES)
    spec_on = spectrum_report(fisher_from_sensitivity(sens_on, w), FULL_NAMES)
    spec_r = spectrum_report(fisher_from_sensitivity(sens_r, w), RED_NAMES)
    rmse_off = float(np.sqrt(np.mean((mu_off - predict(phi, "Dearly", "reduced")) ** 2)))
    rmse_on = float(np.sqrt(np.mean((mu_on - predict(phi, "Dearly", "reduced")) ** 2)))
    return {
        "z0_off": Z_OFF,
        "z0_on": z_on,
        "log_rmse_off_vs_reduced": rmse_off,
        "log_rmse_on_vs_reduced": rmse_on,
        "full_off": spec_off,
        "full_on": spec_on,
        "reduced": spec_r,
    }


def flag_rank(r_full: int, r_red: int | None) -> str:
    if r_red is None:
        return "vanishes"
    if r_red == r_full:
        return "survives"
    if r_red < r_full:
        return "vanishes"
    return "emerges_after_reduction"


def participation_table(spec: dict) -> list[dict]:
    ev = np.array(spec["eigenvectors_desc_columns"])
    rows = []
    for j, ratio in enumerate(spec["ratios_to_max"]):
        comp = ev[:, j] ** 2
        order = np.argsort(comp)[::-1]
        rows.append(
            {
                "index": j + 1,
                "ratio_to_max": float(ratio),
                "squared_components": {spec["names"][i]: float(comp[i]) for i in order},
            }
        )
    return rows


def schedule_record(name: str, pair: dict, observed_on_reduced: bool) -> dict:
    full = pair["full"]
    red = pair["reduced"]
    r_full = int(full["practical_rank"])
    n_full = int(full["n_param"])
    supported = full.get("n_practical_supported_on_reduced_coordinates")
    if not observed_on_reduced or red is None:
        vanished = r_full  # the whole schedule is an unreduced coordinate
        # More precise: eigenvalues supported only on unreduced axes.
        artefact = None if supported is None else int(r_full - supported)
        return {
            "schedule": name,
            "observed_on_reduced_model": False,
            "practical_rank_full": r_full,
            "numerical_rank_full": int(full["numerical_rank"]),
            "n_param_full": n_full,
            "practical_rank_reduced": None,
            "numerical_rank_reduced": None,
            "n_param_reduced": None,
            "n_supported_on_reduced_coordinates": supported,
            "flag": "vanishes",
            "artefact_count": artefact if artefact is not None else r_full,
            "note": "Schedule uses a state the reduction deletes.",
        }
    r_red = int(red["practical_rank"])
    artefact = int(r_full - supported) if supported is not None else max(r_full - r_red, 0)
    return {
        "schedule": name,
        "observed_on_reduced_model": True,
        "practical_rank_full": r_full,
        "numerical_rank_full": int(full["numerical_rank"]),
        "n_param_full": n_full,
        "practical_rank_reduced": r_red,
        "numerical_rank_reduced": int(red["numerical_rank"]),
        "n_param_reduced": int(red["n_param"]),
        "n_supported_on_reduced_coordinates": supported,
        "flag": flag_rank(r_full, r_red),
        "artefact_count": artefact,
        "lambda_max_full": full["lambda_max"],
        "lambda_max_reduced": red["lambda_max"],
        "condition_full": full["condition_of_stiff_block"],
        "condition_reduced": red["condition_of_stiff_block"],
        "log10_eigenvalues_full": full["log10_eigenvalues"],
        "log10_eigenvalues_reduced": red["log10_eigenvalues"],
        "ratios_full": full["ratios_to_max"],
        "ratios_reduced": red["ratios_to_max"],
    }


def steady_table() -> dict:
    rows = {}
    for name in LINEAGES:
        st = steady_full(THETA0, U_G[name], U_Q[name])
        red = steady_reduced(phi_of(THETA0), U_G[name], U_Q[name])
        rows[name] = {
            "u_G": U_G[name],
            "u_Q": U_Q[name],
            "G": float(st[0]),
            "P": float(st[1]),
            "L": float(st[2]),
            "Q": float(st[3]),
            "Z": float(st[4]),
            "L_over_G": float(st[2] / st[0]),
            "max_abs_gap_vs_reduced": float(np.max(np.abs(st[:4] - red))),
        }
    ratios = [rows[n]["L_over_G"] for n in LINEAGES]
    return {
        "by_lineage": rows,
        "ratio_span": float(max(ratios) - min(ratios)),
        "ratio_identical": bool(max(ratios) - min(ratios) < 1e-8),
    }


def _style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelsize": 10,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "figure.dpi": 140,
        }
    )


def make_figures(records: dict, flow: dict, early: dict, loaded: dict) -> None:
    _style()
    # Figure 1. Spectra.
    order = ["S", "M", "Dearly_on", "Dearly_off", "Dlate"]
    fig, axes = plt.subplots(1, len(order), figsize=(11.2, 3.3), sharey=True)
    for ax, key in zip(axes, order):
        rec = records[key]
        yf = rec["log10_eigenvalues_full"]
        yr = rec["log10_eigenvalues_reduced"]
        ax.plot(np.arange(1, len(yf) + 1), yf, "o", color="#1f4e79", ms=4.5, label="full")
        ax.plot(np.arange(1, len(yr) + 1), yr, "s", color="#c47b2b", ms=4.0, label="reduced")
        ax.axhline(np.log10(PRACTICAL_CUT) + yf[0], color="#8c2f39", lw=0.7, ls="--")
        ax.set_title(key.replace("_", "\n"), fontsize=9)
        ax.set_xlabel("index")
        ax.set_xticks(range(1, 9))
    axes[0].set_ylabel("log10 eigenvalue")
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fim_spectra_full_vs_reduced.png", bbox_inches="tight")
    plt.close(fig)

    # Figure 2. Rank bars.
    labels = ["S", "M", "M loaded\nprofiled", "D early\non manifold", "D early\noff manifold", "D late", "Z channel"]
    keys_full = [
        records["S"]["practical_rank_full"],
        records["M"]["practical_rank_full"],
        loaded["full"]["profiled_kinetics"]["practical_rank"],
        records["Dearly_on"]["practical_rank_full"],
        records["Dearly_off"]["practical_rank_full"],
        records["Dlate"]["practical_rank_full"],
        records["Z"]["practical_rank_full"],
    ]
    keys_red = [
        records["S"]["practical_rank_reduced"],
        records["M"]["practical_rank_reduced"],
        loaded["reduced"]["profiled_kinetics"]["practical_rank"],
        records["Dearly_on"]["practical_rank_reduced"],
        records["Dearly_off"]["practical_rank_reduced"],
        records["Dlate"]["practical_rank_reduced"],
        np.nan,
    ]
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(8.4, 3.8))
    ax.bar(x - 0.18, keys_full, width=0.36, color="#1f4e79", label="full, 8 rates")
    ax.bar(x + 0.18, keys_red, width=0.36, color="#c47b2b", label="reduced, 6 rates")
    ax.text(x[-1] + 0.18, 0.15, "deleted", ha="center", va="bottom", fontsize=8, color="#8c2f39")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("practical Fisher rank")
    ax.set_ylim(0, 9)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "rank_bars.png", bbox_inches="tight")
    plt.close(fig)

    # Figure 3. Alignment of full-model eigenvectors for M and Dearly_off.
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.4), sharey=True)
    for ax, key, title in (
        (axes[0], "M", "steady four channels"),
        (axes[1], "Dearly_off", "early time course, Z(0) off"),
        (axes[2], "Z", "steady four channels plus Z"),
    ):
        al = records[key]["alignment_full"]
        idx = np.arange(1, 9)
        ax.plot(idx, al["kappa"], "o-", color="#1f4e79", ms=4, label="κ direction")
        ax.plot(idx, al["rho"], "s-", color="#8c2f39", ms=4, label="ρ direction")
        ax.plot(idx, al["speed"], "^-", color="#2f6f4e", ms=4, label="speed direction")
        ax.set_xlabel("eigenvector index (stiff to sloppy)")
        ax.set_title(title, fontsize=10)
        ax.set_ylim(-0.05, 1.05)
    axes[0].set_ylabel("squared overlap")
    axes[2].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "sloppy_alignment.png", bbox_inches="tight")
    plt.close(fig)

    # Figure 4. Null walk.
    rows = flow["exact_rho_walk"]
    steps = [r["log_step"] for r in rows]
    fig, ax1 = plt.subplots(figsize=(6.2, 3.6))
    ax1.plot(steps, [r["chi_M"] for r in rows], color="#1f4e79", label="χ² on M")
    ax1.plot(steps, [r["chi_Z"] for r in rows], color="#8c2f39", label="χ² on Z")
    ax1.set_xlabel("arc length along the exact ρ direction")
    ax1.set_ylabel("noise-weighted χ²")
    ax2 = ax1.twinx()
    ax2.plot(steps, [r["rho"] / rows[0]["rho"] for r in rows], color="#c47b2b", ls="--", label="ρ / ρ₀")
    ax2.plot(steps, [r["kappa"] / rows[0]["kappa"] for r in rows], color="#2f6f4e", ls=":", label="κ / κ₀")
    ax2.set_ylabel("ratio to the generating value")
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "rho_null_walk.png", bbox_inches="tight")
    plt.close(fig)

    # Figure 5. One lineage, lactate, full versus reduced.
    times = np.linspace(0.0, 8.0, 240)
    name = "BRCA"
    y_off = np.concatenate([Y_GPLQ, [Z_OFF]])
    y_on = np.concatenate([Y_GPLQ, [rho_of(THETA0) * Y_GPLQ[1]]])
    tr_off = _integrate(rhs_full, y_off, THETA0, U_G[name], U_Q[name], times)
    tr_on = _integrate(rhs_full, y_on, THETA0, U_G[name], U_Q[name], times)
    tr_red = _integrate(rhs_reduced, Y_GPLQ, phi_of(THETA0), U_G[name], U_Q[name], times)
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.plot(times, tr_off[2], color="#8c2f39", label="full, Z(0) off")
    ax.plot(times, tr_on[2], color="#1f4e79", label="full, Z(0) on manifold")
    ax.plot(times, tr_red[2], color="#c47b2b", ls="--", label="reduced")
    for t in T_EARLY:
        ax.axvline(t, color="#bbbbbb", lw=0.5)
    ax.set_xlabel("time")
    ax.set_ylabel("lactate, BRCA input")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "lactate_initial_layer.png", bbox_inches="tight")
    plt.close(fig)
    del early


def noisy_check() -> dict:
    """One draw on schedule M. Not used for the rank table."""
    mu = predict(THETA0, "M", "full")
    # Whiten, add standard normal, unwhiten. Replicates are not stacked; one panel mean.
    cov = (SIGMA["M"] ** 2) * CORR
    chol = cholesky(cov, lower=True)
    draws = []
    for _name in LINEAGES:
        z = chol @ RNG.standard_normal(4)
        draws.append(z)
    noise = np.concatenate(draws)
    # Compare the score of the true point to a point moved along rho.
    direction = analytic_axes()["rho"]
    th = THETA0 * np.exp(0.4 * direction)
    mu_moved = predict(th, "M", "full")
    y = mu + noise
    return {
        "seed": SEED,
        "chi_true_point": chi_square(mu, y, "M"),
        "chi_rho_moved": chi_square(mu_moved, y, "M"),
        "chi_gap": chi_square(mu_moved, y, "M") - chi_square(mu, y, "M"),
    }


def main() -> None:
    steady = steady_table()
    gap = max(row["max_abs_gap_vs_reduced"] for row in steady["by_lineage"].values())
    if gap > 1e-9:
        raise RuntimeError(f"steady reduction is not exact: gap {gap}")

    print("steady schedules")
    pair_s = fisher_pair("S", FD_STEADY)
    pair_m = fisher_pair("M", FD_STEADY)
    print("Z schedule")
    pair_z = fisher_pair("Z", FD_STEADY)
    print("late time course")
    pair_late = fisher_pair("Dlate", FD_TRANSIENT)
    print("early layer")
    early = early_initial_layer()
    print("loadings")
    loaded = {"full": profiled_kinetic_spectrum("full"), "reduced": profiled_kinetic_spectrum("reduced")}
    print("flow and bad edit")
    flow = null_flow(pair_m["full"])
    bad = bad_reduction_check()
    noisy = noisy_check()
    angles = principal_angles(pair_m["full"])

    def pack_early(spec_full: dict, spec_red: dict) -> dict:
        return {"full": spec_full, "reduced": spec_red, "schedule": "Dearly"}

    records_raw = {
        "S": schedule_record("S", pair_s, True),
        "M": schedule_record("M", pair_m, True),
        "Dlate": schedule_record("Dlate", pair_late, True),
        "Dearly_off": schedule_record("Dearly_off", pack_early(early["full_off"], early["reduced"]), True),
        "Dearly_on": schedule_record("Dearly_on", pack_early(early["full_on"], early["reduced"]), True),
        "Z": schedule_record("Z", pair_z, False),
    }
    # Attach alignment for the two figures that need it, then strip bulky eigenvectors from the archive.
    for key, spec in (
        ("M", pair_m["full"]),
        ("Dearly_off", early["full_off"]),
        ("Dearly_on", early["full_on"]),
        ("S", pair_s["full"]),
        ("Dlate", pair_late["full"]),
        ("Z", pair_z["full"]),
    ):
        records_raw[key]["alignment_full"] = spec["alignment_sq_with_analytic_axes"]

    r_m = records_raw["M"]["practical_rank_full"]
    r_s = records_raw["S"]["practical_rank_full"]
    r_m_red = records_raw["M"]["practical_rank_reduced"]
    r_s_red = records_raw["S"]["practical_rank_reduced"]
    delta_full = r_m - r_s
    delta_red = r_m_red - r_s_red
    if delta_full > 0 and delta_red == delta_full:
        advantage = "survives"
    elif delta_full > 0 and delta_red == 0:
        advantage = "vanishes"
    elif delta_full > 0 and delta_red < delta_full:
        advantage = "shrinks"
    elif delta_red > delta_full:
        advantage = "grows"
    else:
        advantage = "absent"

    loaded_flag = flag_rank(
        loaded["full"]["profiled_kinetics"]["practical_rank"],
        loaded["reduced"]["profiled_kinetics"]["practical_rank"],
    )

    summary_rows = []
    for key in ["S", "M", "Dlate", "Dearly_on", "Dearly_off", "Z"]:
        rec = records_raw[key]
        summary_rows.append(
            {
                "schedule": key,
                "practical_rank_full": rec["practical_rank_full"],
                "n_full": rec["n_param_full"],
                "practical_rank_reduced": rec["practical_rank_reduced"],
                "n_reduced": rec["n_param_reduced"],
                "numerical_rank_full": rec["numerical_rank_full"],
                "numerical_rank_reduced": rec["numerical_rank_reduced"],
                "supported": rec["n_supported_on_reduced_coordinates"],
                "practical_flag": rec["flag"],
                "numerical_flag": flag_rank(rec["numerical_rank_full"], rec["numerical_rank_reduced"]),
                "artefact_count": rec["artefact_count"],
            }
        )

    archive = {
        "seed": SEED,
        "data": "synthetic multi-channel surrogate, not a CCLE download",
        "depends_on": ["T09", "T12"],
        "theta_full": THETA0.tolist(),
        "theta_names_full": FULL_NAMES,
        "phi_reduced": phi_of(THETA0).tolist(),
        "phi_names_reduced": RED_NAMES,
        "kappa": kappa_of(THETA0),
        "rho": rho_of(THETA0),
        "cuts": {"numerical": NUM_CUT, "practical": PRACTICAL_CUT, "support_mass": SUPPORT_MASS},
        "noise": {"sigma": SIGMA, "n_rep": N_REP, "correlation": CORR.tolist()},
        "times": {"early": T_EARLY.tolist(), "late": T_LATE.tolist()},
        "steady": steady,
        "reduction_steady_gap": gap,
        "ranks": summary_rows,
        "channel_advantage_M_minus_S": {
            "full": delta_full,
            "reduced": delta_red,
            "flag": advantage,
        },
        "loaded_profiled_flag": loaded_flag,
        "loaded": {
            "full_joint_practical_rank": loaded["full"]["joint"]["practical_rank"],
            "full_joint_n": loaded["full"]["joint"]["n_param"],
            "reduced_joint_practical_rank": loaded["reduced"]["joint"]["practical_rank"],
            "reduced_joint_n": loaded["reduced"]["joint"]["n_param"],
            "full_profiled_practical_rank": loaded["full"]["profiled_kinetics"]["practical_rank"],
            "reduced_profiled_practical_rank": loaded["reduced"]["profiled_kinetics"]["practical_rank"],
            "full_profiled_log10": loaded["full"]["profiled_kinetics"]["log10_eigenvalues"],
            "reduced_profiled_log10": loaded["reduced"]["profiled_kinetics"]["log10_eigenvalues"],
            "full_profiled_ratios": loaded["full"]["profiled_kinetics"]["ratios_to_max"],
            "reduced_profiled_ratios": loaded["reduced"]["profiled_kinetics"]["ratios_to_max"],
            "full_profiled_numerical_rank": loaded["full"]["profiled_kinetics"]["numerical_rank"],
            "reduced_profiled_numerical_rank": loaded["reduced"]["profiled_kinetics"]["numerical_rank"],
            "full_joint_numerical_rank": loaded["full"]["joint"]["numerical_rank"],
            "reduced_joint_numerical_rank": loaded["reduced"]["joint"]["numerical_rank"],
        },
        "participation": {
            "S_full": participation_table(pair_s["full"]),
            "S_reduced": participation_table(pair_s["reduced"]),
            "M_full": participation_table(pair_m["full"]),
            "M_reduced": participation_table(pair_m["reduced"]),
            "Z_full": participation_table(pair_z["full"]),
            "Dearly_off_full": participation_table(early["full_off"]),
            "Dlate_full": participation_table(pair_late["full"]),
        },
        "principal_angles_M_sloppy_vs_analytic_plane": angles,
        "early_layer": {
            "z0_off": early["z0_off"],
            "z0_on": early["z0_on"],
            "log_rmse_off_vs_reduced": early["log_rmse_off_vs_reduced"],
            "log_rmse_on_vs_reduced": early["log_rmse_on_vs_reduced"],
        },
        "flow": {
            "trailing_squared_overlap_with_rho": flow["trailing_squared_overlap_with_rho"],
            "trailing_squared_overlap_with_speed": flow["trailing_squared_overlap_with_speed"],
            "trailing_squared_overlap_with_kappa": flow["trailing_squared_overlap_with_kappa"],
            "trailing_dlog_kappa": flow["trailing_dlog_kappa"],
            "exact_rho_walk": flow["exact_rho_walk"],
        },
        "bad_reduction": bad,
        "noisy_check": noisy,
        "spectra": {
            key: {
                "log10_full": records_raw[key]["log10_eigenvalues_full"] if "log10_eigenvalues_full" in records_raw[key] else records_raw[key].get("log10_eigenvalues_full"),
                "log10_reduced": records_raw[key].get("log10_eigenvalues_reduced"),
                "ratios_full": records_raw[key].get("ratios_full"),
                "ratios_reduced": records_raw[key].get("ratios_reduced"),
                "condition_full": records_raw[key].get("condition_full"),
                "condition_reduced": records_raw[key].get("condition_reduced"),
                "lambda_max_full": records_raw[key].get("lambda_max_full"),
                "lambda_max_reduced": records_raw[key].get("lambda_max_reduced"),
                "alignment_full": records_raw[key].get("alignment_full"),
            }
            for key in records_raw
        },
        "disclaimer": "Synthetic generator. Not a medical device, dose, or cell-line measurement.",
    }
    # Z record has no reduced eigenvalues. Fill spectra from the full spec only.
    archive["spectra"]["Z"]["log10_full"] = pair_z["full"]["log10_eigenvalues"]
    archive["spectra"]["Z"]["ratios_full"] = pair_z["full"]["ratios_to_max"]
    archive["spectra"]["Z"]["condition_full"] = pair_z["full"]["condition_of_stiff_block"]
    archive["spectra"]["Z"]["lambda_max_full"] = pair_z["full"]["lambda_max"]

    make_figures(records_raw | {
        "M": records_raw["M"] | {"log10_eigenvalues_full": pair_m["full"]["log10_eigenvalues"], "log10_eigenvalues_reduced": pair_m["reduced"]["log10_eigenvalues"], "alignment_full": pair_m["full"]["alignment_sq_with_analytic_axes"]},
        "S": records_raw["S"] | {"log10_eigenvalues_full": pair_s["full"]["log10_eigenvalues"], "log10_eigenvalues_reduced": pair_s["reduced"]["log10_eigenvalues"]},
        "Dlate": records_raw["Dlate"] | {"log10_eigenvalues_full": pair_late["full"]["log10_eigenvalues"], "log10_eigenvalues_reduced": pair_late["reduced"]["log10_eigenvalues"]},
        "Dearly_off": records_raw["Dearly_off"] | {"log10_eigenvalues_full": early["full_off"]["log10_eigenvalues"], "log10_eigenvalues_reduced": early["reduced"]["log10_eigenvalues"], "alignment_full": early["full_off"]["alignment_sq_with_analytic_axes"]},
        "Dearly_on": records_raw["Dearly_on"] | {"log10_eigenvalues_full": early["full_on"]["log10_eigenvalues"], "log10_eigenvalues_reduced": early["reduced"]["log10_eigenvalues"]},
    }, flow, early, loaded)

    out = ROOT / "results.json"
    # Drop eigenvectors if any leaked.
    text = json.dumps(archive, indent=2)
    out.write_text(text + "\n", encoding="utf-8")
    print(json.dumps({"ranks": summary_rows, "advantage": archive["channel_advantage_M_minus_S"], "loaded_flag": loaded_flag, "angles": angles, "bad": bad, "early_rmse": archive["early_layer"], "kappa": archive["kappa"], "steady_ratios": {n: steady["by_lineage"][n]["L_over_G"] for n in LINEAGES}}, indent=2))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
