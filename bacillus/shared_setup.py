# ============================================================
# shared_setup.py
# Import this in every figure file using:
#   from shared_setup import *
# This defines all parameters, ODE functions, and run_simulation.
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os

# Figures are written to an outputs/ folder beside this file.
OUT = "./outputs"
os.makedirs(OUT, exist_ok=True)

# --- Illustrative baseline parameters (not experimentally validated) ---
mu_B    = 0.40
beta_B  = 0.10   # CORRECTED: the zip export had 0.70, which does not
                 # reproduce the published figures. The team Colab and
                 # the Wiki parameter table both use 0.10, which gives
                 # the optimum at t_i = 6.5 h and 198 mg/L as published.
q_B     = 2.00
kd0_B   = 0.010
kprot_B = 0.005

KS_X_B  = 0.10
KN_X_B  = 0.05
KO_X_B  = 0.10

KS_P_B  = 0.10
KN_P_B  = 0.05
KO_P_B  = 0.10

YX_S_B  = 0.50
YX_N_B  = 4.00
YX_O_B  = 0.50

YP_S_B  = 100.0
YP_N_B  = 800.0
YP_O_B  = 100.0

kLa_B   = 0.30
alpha_B = 0.03
V_l_B   = 1.00
V_g_B   = 0.25

X0_B    = 0.05
P0_B    = 0.00
S0_B    = 20.0
N0_B    = 5.00
O0_B    = 7.00
Og0_B   = 270.0

T       = 72.0
t_eval  = np.linspace(0, T, 1000)
EPS     = 1e-9

# --- Induction switch ---
def mu_effective(t, t_ind):
    return mu_B if t < t_ind else beta_B * mu_B

def q_effective(t, t_ind):
    return 0.0 if t < t_ind else q_B

# --- Monod limitation ---
def Phi_X(S, N, O):
    s = max(S, 0.0) / (KS_X_B + max(S, 0.0) + EPS)
    n = max(N, 0.0) / (KN_X_B + max(N, 0.0) + EPS)
    o = max(O, 0.0) / (KO_X_B + max(O, 0.0) + EPS)
    return s * n * o

def Phi_P(S, N, O):
    s = max(S, 0.0) / (KS_P_B + max(S, 0.0) + EPS)
    n = max(N, 0.0) / (KN_P_B + max(N, 0.0) + EPS)
    o = max(O, 0.0) / (KO_P_B + max(O, 0.0) + EPS)
    return s * n * o

# --- ODE system ---
def bacillus_model(t, y, t_ind):
    X, P, S, N, O, Og = y
    X  = max(X,  0.0); P  = max(P,  0.0)
    S  = max(S,  0.0); N  = max(N,  0.0)
    O  = max(O,  0.0); Og = max(Og, 0.0)

    mu_e = mu_effective(t, t_ind)
    q_e  = q_effective(t, t_ind)
    phiX = Phi_X(S, N, O)
    phiP = Phi_P(S, N, O)

    dXdt  = mu_e * phiX * X
    dPdt  = q_e * phiP * X - (kd0_B + kprot_B * X) * P
    dSdt  = -((mu_e * phiX / YX_S_B) + (q_e * phiP / YP_S_B)) * X
    dNdt  = -((mu_e * phiX / YX_N_B) + (q_e * phiP / YP_N_B)) * X
    driving = alpha_B * Og - O
    dOdt  = kLa_B * driving - ((mu_e * phiX / YX_O_B) + (q_e * phiP / YP_O_B)) * X
    dOgdt = -(kLa_B * V_l_B / V_g_B) * driving
    return [dXdt, dPdt, dSdt, dNdt, dOdt, dOgdt]

# --- Solver ---
def run_simulation(t_ind):
    y0  = [X0_B, P0_B, S0_B, N0_B, O0_B, Og0_B]
    sol = solve_ivp(
        fun=lambda t, y: bacillus_model(t, y, t_ind),
        t_span=(0.0, T), y0=y0, t_eval=t_eval,
        method='RK45', rtol=1e-6, atol=1e-9)
    return sol.t, np.maximum(sol.y, 0.0)

def vline(ax, t_ind):
    ax.axvline(t_ind, color='red', linestyle='--', linewidth=1.4,
               label=f'Induction $t_i={t_ind:.0f}$ h')
