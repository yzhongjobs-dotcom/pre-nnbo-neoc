"""Derivative consistency checks for the DS-ANEOC-NNBO HC model.

This script independently compares the closed-form derivatives used in the
method with a small forward-mode AD implementation and central differences.
It covers F_theta, G_theta, f_s,theta, h_theta, A_theta, and the reduced-HJB
direct parameter derivative for the shifted contact equilibrium. Galerkin-root and policy-tangent checks require an
actual basis/collocation anchor and are specified as D2-D3 in the method.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np


P = 3


@dataclass
class Dual:
    val: float
    der: np.ndarray

    def __add__(self, other):
        other = as_dual(other)
        return Dual(self.val + other.val, self.der + other.der)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.val, -self.der)

    def __sub__(self, other):
        return self + (-as_dual(other))

    def __rsub__(self, other):
        return as_dual(other) + (-self)

    def __mul__(self, other):
        other = as_dual(other)
        return Dual(self.val * other.val,
                    self.der * other.val + other.der * self.val)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = as_dual(other)
        return Dual(self.val / other.val,
                    (self.der * other.val - other.der * self.val)
                    / (other.val ** 2))

    def __rtruediv__(self, other):
        return as_dual(other) / self


def as_dual(value) -> Dual:
    if isinstance(value, Dual):
        return value
    return Dual(float(value), np.zeros(P))


def dexp(x: Dual) -> Dual:
    value = math.exp(x.val)
    return Dual(value, value * x.der)


def independent(theta: np.ndarray) -> list[Dual]:
    out = []
    for j, value in enumerate(theta):
        derivative = np.zeros(P)
        derivative[j] = 1.0
        out.append(Dual(float(value), derivative))
    return out


def model_values(theta, e_tilde, r_tilde, psi, gamma,
                 delta_0, delta_s, x_dot):
    kappa, b, n = theta
    bar_delta_s = delta_s / delta_0
    bar_delta = bar_delta_s + (e_tilde + psi * r_tilde) / delta_0
    rho = bar_delta_s / bar_delta
    K = np.exp(kappa)
    B = np.exp(b)
    a = K / B
    F = np.array([a * (rho ** n - 1.0) - psi * gamma * r_tilde,
                  gamma * r_tilde])
    G = np.array([np.exp(-b) * bar_delta ** (-n), 0.0])
    h = bar_delta ** n * (K + B * x_dot)
    f_s = K * bar_delta_s ** n
    return F, G, h, f_s


def model_ad(theta, e_tilde, r_tilde, psi, gamma,
             delta_0, delta_s, x_dot):
    kappa, b, n = independent(theta)
    bar_delta_s = delta_s / delta_0
    bar_delta = bar_delta_s + (e_tilde + psi * r_tilde) / delta_0
    rho = bar_delta_s / bar_delta
    K = dexp(kappa)
    B = dexp(b)
    a = K / B
    bar_pow_n = dexp(n * math.log(bar_delta))
    rho_pow_n = dexp(n * math.log(rho))
    F0 = a * (rho_pow_n - 1.0) - psi * gamma * r_tilde
    F1 = as_dual(gamma * r_tilde)
    G0 = dexp(-b) / bar_pow_n
    G1 = as_dual(0.0)
    h = bar_pow_n * (K + B * x_dot)
    f_s = K * dexp(n * math.log(bar_delta_s))
    F = np.array([F0.val, F1.val])
    G = np.array([G0.val, G1.val])
    JF = np.vstack([F0.der, F1.der])
    JG = np.vstack([G0.der, G1.der])
    return F, G, h.val, f_s.val, JF, JG, h.der, f_s.der


def analytic_jacobians(theta, e_tilde, r_tilde, psi,
                       delta_0, delta_s, x_dot):
    kappa, b, n = theta
    bar_delta_s = delta_s / delta_0
    bar_delta = bar_delta_s + (e_tilde + psi * r_tilde) / delta_0
    rho = bar_delta_s / bar_delta
    K = math.exp(kappa)
    B = math.exp(b)
    a = math.exp(kappa - b)
    c = rho ** n
    G0 = math.exp(-b) * bar_delta ** (-n)
    h = bar_delta ** n * (K + B * x_dot)
    f_s = K * bar_delta_s ** n
    JF = np.array([[a * (c - 1.0), -a * (c - 1.0), a * c * math.log(rho)],
                   [0.0, 0.0, 0.0]])
    JG = np.array([[0.0, -G0, -math.log(bar_delta) * G0],
                   [0.0, 0.0, 0.0]])
    Jh = np.array([bar_delta ** n * K,
                   bar_delta ** n * B * x_dot,
                   h * math.log(bar_delta)])
    Jfs = np.array([f_s, 0.0, f_s * math.log(bar_delta_s)])
    return JF, JG, Jh, Jfs


def central_jacobian(fun, theta, step):
    value = np.atleast_1d(fun(theta))
    jacobian = np.zeros((value.size, theta.size))
    for j in range(theta.size):
        direction = np.zeros_like(theta)
        direction[j] = step
        jacobian[:, j] = (
            np.atleast_1d(fun(theta + direction))
            - np.atleast_1d(fun(theta - direction))
        ) / (2.0 * step)
    return jacobian


def relative_error(actual, reference):
    return np.linalg.norm(actual - reference) / max(1.0, np.linalg.norm(reference))


def run(seed=20260909, trials=64):
    rng = np.random.default_rng(seed)
    max_ad = {"F": 0.0, "G": 0.0, "h": 0.0, "f_s": 0.0,
              "A": 0.0, "HJB": 0.0}
    max_fd = {key: 0.0 for key in max_ad}
    steps = (3e-4, 1e-4, 3e-5)
    convergence_samples = []

    for _ in range(trials):
        theta = np.array([rng.uniform(-0.5, 1.0),
                          rng.uniform(-1.0, 0.5),
                          rng.uniform(0.8, 1.8)])
        e_tilde = rng.uniform(-0.002, 0.002)
        r_tilde = rng.uniform(-0.002, 0.002)
        psi = 1.0
        gamma = rng.uniform(-1.5, -0.2)
        delta_0 = 0.01
        delta_s = rng.uniform(0.008, 0.015)
        # Keep the normalized indentation strictly positive.
        if delta_s + e_tilde + psi * r_tilde <= 0.002:
            e_tilde += 0.004
        x_dot = rng.uniform(-0.03, 0.03)

        args = (e_tilde, r_tilde, psi, gamma, delta_0, delta_s, x_dot)
        F, G, h, f_s, JF_ad, JG_ad, Jh_ad, Jfs_ad = model_ad(theta, *args)
        JF, JG, Jh, Jfs = analytic_jacobians(
            theta, e_tilde, r_tilde, psi, delta_0, delta_s, x_dot
        )
        max_ad["F"] = max(max_ad["F"], relative_error(JF_ad, JF))
        max_ad["G"] = max(max_ad["G"], relative_error(JG_ad, JG))
        max_ad["h"] = max(max_ad["h"], relative_error(Jh_ad, Jh))
        max_ad["f_s"] = max(max_ad["f_s"], relative_error(Jfs_ad, Jfs))

        value_fun = lambda th: model_values(th, *args)
        for step in steps:
            JF_fd = central_jacobian(lambda th: value_fun(th)[0], theta, step)
            JG_fd = central_jacobian(lambda th: value_fun(th)[1], theta, step)
            Jh_fd = central_jacobian(lambda th: value_fun(th)[2], theta, step).ravel()
            Jfs_fd = central_jacobian(lambda th: value_fun(th)[3], theta, step).ravel()
            max_fd["F"] = max(max_fd["F"], relative_error(JF_fd, JF))
            max_fd["G"] = max(max_fd["G"], relative_error(JG_fd, JG))
            max_fd["h"] = max(max_fd["h"], relative_error(Jh_fd, Jh))
            max_fd["f_s"] = max(max_fd["f_s"], relative_error(Jfs_fd, Jfs))

        # A=G R^{-1}G^T and its direct theta derivative.
        R = rng.uniform(0.5, 2.0)
        A = np.outer(G, G) / R
        JA = np.stack([
            (np.outer(JG[:, j], G) + np.outer(G, JG[:, j])) / R
            for j in range(P)
        ], axis=2)
        JG_dual = JG_ad
        JA_ad = np.stack([
            (np.outer(JG_dual[:, j], G) + np.outer(G, JG_dual[:, j])) / R
            for j in range(P)
        ], axis=2)
        max_ad["A"] = max(max_ad["A"], relative_error(JA_ad, JA))
        JA_fd = np.stack([
            central_jacobian(
                lambda th: (np.outer(value_fun(th)[1], value_fun(th)[1]) / R).ravel(),
                theta,
                step,
            ).reshape(2, 2, P)
            for step in (steps[-1],)
        ])[0]
        max_fd["A"] = max(max_fd["A"], relative_error(JA_fd, JA))

        # Direct parameter derivative of p^T F - 1/4 p^T A p, with p fixed.
        p_value = rng.normal(size=2)
        hjb_jac = p_value @ JF - 0.25 * np.einsum("i,ijk,j->k", p_value, JA, p_value)
        hjb_ad = p_value @ JF_ad - 0.25 * np.einsum(
            "i,ijk,j->k", p_value, JA_ad, p_value
        )
        max_ad["HJB"] = max(max_ad["HJB"], relative_error(hjb_ad, hjb_jac))

        hjb_fun = lambda th: (
            p_value @ value_fun(th)[0]
            - 0.25 * p_value @ (
                np.outer(value_fun(th)[1], value_fun(th)[1]) / R
            ) @ p_value
        )
        errors = []
        for step in steps:
            hjb_fd = central_jacobian(hjb_fun, theta, step).ravel()
            err = relative_error(hjb_fd, hjb_jac)
            errors.append(err)
            max_fd["HJB"] = max(max_fd["HJB"], err)
        convergence_samples.append(errors)

    # AD should agree near machine precision. The largest listed central-difference
    # error includes the coarsest step; the finest-step threshold is checked below.
    convergence_samples = np.asarray(convergence_samples)
    finest_hjb_errors = convergence_samples[:, -1]
    median_hjb_errors = np.median(convergence_samples, axis=0)
    observed_orders = [
        math.log(median_hjb_errors[j] / median_hjb_errors[j + 1])
        / math.log(steps[j] / steps[j + 1])
        for j in range(len(steps) - 1)
    ]
    ad_pass = max(max_ad.values()) < 5e-12
    fd_pass = max(max_fd.values()) < 1e-6 and float(np.max(finest_hjb_errors)) < 5e-8
    order_pass = min(observed_orders) > 1.7

    print("DS-ANEOC-NNBO derivative verification")
    print(f"seed={seed}, trials={trials}, central-difference steps={steps}")
    print("max relative error, analytic vs forward AD:")
    for key, value in max_ad.items():
        print(f"  {key:>4s}: {value:.3e}")
    print("max relative error, analytic vs central difference (all tested steps):")
    for key, value in max_fd.items():
        print(f"  {key:>4s}: {value:.3e}")
    print(f"max HJB error at finest step: {np.max(finest_hjb_errors):.3e}")
    print("median HJB error by step:",
          ", ".join(f"{value:.3e}" for value in median_hjb_errors))
    print("observed central-difference orders:",
          ", ".join(f"{value:.2f}" for value in observed_orders))
    passed = ad_pass and fd_pass and order_pass
    print("RESULT:", "PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(run())
