# ============================================================
# yeast_shared.py
# Shared parameters and functions for the yeast (S. cerevisiae)
# nanobody production models, Models 1 and 2.
#
# Transcribed verbatim from the team's "Modeling-Yeast" Colab
# notebook. Import in any yeast figure script with:
#     from yeast_shared import *
# ============================================================

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

OUT = os.environ.get("YEAST_OUT", "./outputs")
os.makedirs(OUT, exist_ok=True)

# ---------------- Model 1 / Model 2 shared parameters ----------------
X0   = 0.1       # g/L      initial biomass
mu   = 0.0854    # h^-1     max specific growth rate (FBA-informed)
beta = 0.5       # -        post-induction growth factor
q    = 0.05      # mg/g/h   specific nanobody production rate
kd   = 0.01      # h^-1     product degradation / loss rate
T    = 48        # h        total cultivation time

# ---------------- Model 2 only: Monod constants ----------------
KS_X = 0.5
KO_X = 0.05
KS_P = 0.5
KO_P = 0.05

# ---------------- Model 2 only: yield coefficients ----------------
YX_S = 0.5
YP_S = 0.2
YX_O = 0.5
YP_O = 0.2

# ---------------- Model 2 only: oxygen transfer / reactor ----------------
kLa   = 0.2
alpha = 1.0
V_l   = 1.0
V_g   = 1.0

# ---------------- Model 2 only: initial conditions ----------------
P0  = 0.0
S0  = 10.0
O0  = 1.0
Og0 = 1.0

t_ind_default = 10   # h, example induction time used for single runs

# ---------------- Plot colours used throughout ----------------
blue   = "#1f77b4"
red    = "#d62728"
green  = "#2ca02c"
purple = "#9467bd"
orange = "#ff7f0e"
black  = "#000000"

t = np.linspace(0, T, 600)


# ====================================================================
#  MODEL 1  --  analytical induction model (no resource limitation)
# ====================================================================

def biomass(t, t_i):
    """Analytical biomass X(t) with a growth-rate switch at t_i."""
    X = np.zeros_like(t)
    for j, tj in enumerate(t):
        if tj < t_i:
            X[j] = X0 * np.exp(mu * tj)
        else:
            X[j] = X0 * np.exp(mu * (t_i + beta * (tj - t_i)))
    return X


def nanobody(t, t_i):
    """Analytical nanobody P(t); zero before induction."""
    P = np.zeros_like(t)
    X_i = X0 * np.exp(mu * t_i)
    for j, tj in enumerate(t):
        if tj < t_i:
            P[j] = 0
        else:
            tau = tj - t_i
            P[j] = (q * X_i / (beta * mu + kd)) * (
                np.exp(beta * mu * tau) - np.exp(-kd * tau)
            )
    return P


def final_nanobody(t_i):
    """Harvest titre P(T) as a function of induction time."""
    if t_i >= T:
        return 0
    tau = T - t_i
    X_i = X0 * np.exp(mu * t_i)
    return (q * X_i / (beta * mu + kd)) * (
        np.exp(beta * mu * tau) - np.exp(-kd * tau)
    )


def optimal_induction_time():
    """Closed-form optimum: returns (tau*, t_i*)."""
    tau_star = (1 / (beta * mu + kd)) * np.log((mu + kd) / ((1 - beta) * mu))
    t_i_star = max(0, min(T, T - tau_star))
    return tau_star, t_i_star


# ====================================================================
#  MODEL 2  --  five-state resource-depleting bioreactor model
#  States: [X, P, S, O, Og]
# ====================================================================

def monod(S, O):
    S = max(S, 0)
    O = max(O, 0)
    if S <= 0 or O <= 0:
        return 0.0, 0.0
    Phi_X = (S / (KS_X + S)) * (O / (KO_X + O))
    Phi_P = (S / (KS_P + S)) * (O / (KO_P + O))
    return Phi_X, Phi_P


def full_5ode(t, y, t_ind):
    X, P, S, O, Og = y
    X  = max(X,  0)
    P  = max(P,  0)
    S  = max(S,  0)
    O  = max(O,  0)
    Og = max(Og, 0)

    mu_e = mu if t < t_ind else beta * mu
    q_e  = 0.0 if t < t_ind else q

    Phi_X, Phi_P = monod(S, O)

    dXdt = mu_e * Phi_X * X
    dPdt = q_e * Phi_P * X - kd * P
    dSdt = -((mu_e * Phi_X / YX_S) + (q_e * Phi_P / YP_S)) * X

    oxygen_transfer    = kLa * (alpha * Og - O)
    oxygen_consumption = ((mu_e * Phi_X / YX_O) + (q_e * Phi_P / YP_O)) * X

    dOdt  = oxygen_transfer - oxygen_consumption
    dOgdt = -(kLa * V_l / V_g) * (alpha * Og - O)

    return [dXdt, dPdt, dSdt, dOdt, dOgdt]


def run_5ode(t_ind, rtol=1e-7, atol=1e-9, n_points=600):
    """Integrate Model 2. Single-trajectory tolerances by default."""
    y0 = [X0, P0, S0, O0, Og0]
    t_eval = np.linspace(0, T, n_points)
    sol = solve_ivp(full_5ode, (0, T), y0, t_eval=t_eval,
                    method="RK45", rtol=rtol, atol=atol, args=(t_ind,))
    return sol.t, sol.y
