# ============================================================
# ode_model_solver.py
# Bacillus subtilis Nanobody Production — ODE Model & Solver
#
# 6 state variables:
#   X   — biomass (g/L)
#   P   — nanobody (mg/L)
#   S   — glucose/carbon source (g/L)
#   N   — nitrogen/amino-acid precursor pool (g/L)
#   O   — dissolved oxygen in liquid (mg/L)
#   Og  — oxygen in gas headspace (mg/L)
#
# Bacillus-specific features:
#   • N(t) state variable (amino-acid precursor pool)
#   • kprot * X * P degradation term (Bacillus extracellular protease)
#   • Multiplicative Monod kinetics for Phi_X and Phi_P
#   • Induction switch: piecewise mu_e and q_e at t = t_ind
# ============================================================

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import os

# ------------------------------------------------------------------
# Output directory
# ------------------------------------------------------------------
OUT = '/mnt/user-data/outputs/figs'
os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------

# Growth / production
mu_B    = 0.40    # max specific growth rate (1/h)
beta_B  = 0.70    # fraction of mu_e active post-induction (–)
q_B     = 2.00    # max specific nanobody production rate (mg/g/h)
kd0_B   = 0.010   # baseline nanobody degradation rate (1/h)
kprot_B = 0.005   # protease degradation coefficient (L/g/h)  [Bacillus-specific †]

# Monod half-saturation constants — growth
KS_X_B  = 0.10   # glucose Ks for growth (g/L)
KN_X_B  = 0.05   # nitrogen Ks for growth (g/L)           [Bacillus-specific †]
KO_X_B  = 0.10   # oxygen Ks for growth (mg/L)

# Monod half-saturation constants — production
KS_P_B  = 0.10   # glucose Ks for production (g/L)
KN_P_B  = 0.05   # nitrogen Ks for production (g/L)        [Bacillus-specific †]
KO_P_B  = 0.10   # oxygen Ks for production (mg/L)

# Yield coefficients — growth (per gram biomass formed)
YX_S_B  = 0.50   # biomass yield on glucose (g X / g S)
YX_N_B  = 4.00   # biomass yield on nitrogen (g X / g N)   [Bacillus-specific †]
YX_O_B  = 0.50   # biomass yield on oxygen (g X / g O)

# Yield coefficients — production (per mg nanobody formed)
YP_S_B  = 100.0  # glucose demand per mg P (mg S / mg P)
YP_N_B  = 800.0  # nitrogen demand per mg P (mg N / mg P)  [Bacillus-specific †]
YP_O_B  = 100.0  # oxygen demand per mg P (mg O / mg P)

# Oxygen transfer
kLa_B   = 0.30   # volumetric O2 transfer coefficient (1/h)
alpha_B = 0.03   # gas-phase O2 depletion rate constant (1/h)
V_l_B   = 1.00   # liquid volume (L)
V_g_B   = 0.25   # gas headspace volume (L)

# Initial conditions
X0_B  = 0.05    # g/L
P0_B  = 0.00    # mg/L
S0_B  = 20.0    # g/L
N0_B  = 5.00    # g/L
O0_B  = 7.00    # mg/L
Og0_B = 270.0   # mg/L

# Simulation time
T      = 72.0                        # total run time (h)
t_eval = np.linspace(0, T, 1000)     # evaluation points

EPS = 1e-9   # numerical floor to avoid division by zero

# ------------------------------------------------------------------
# Resource limitation terms  (multiplicative Monod)
# ------------------------------------------------------------------

def Phi_X(S, N, O):
    """Growth limitation factor (0–1)."""
    return (S / (KS_X_B + S + EPS)) * \
           (N / (KN_X_B + N + EPS)) * \
           (O / (KO_X_B + O + EPS))

def Phi_P(S, N, O):
    """Production limitation factor (0–1)."""
    return (S / (KS_P_B + S + EPS)) * \
           (N / (KN_P_B + N + EPS)) * \
           (O / (KO_P_B + O + EPS))

# ------------------------------------------------------------------
# ODE system
# ------------------------------------------------------------------

def bacillus_ode(t, y, t_ind):
    """
    Right-hand side of the 6-ODE Bacillus nanobody model.

    State vector y = [X, P, S, N, O, Og]

    Parameters
    ----------
    t     : float — current time (h)
    y     : array-like, length 6 — current state
    t_ind : float — induction time (h)
    """
    X, P, S, N, O, Og = y

    # Clamp to non-negative (avoid numerical drift below zero)
    X  = max(X,  0.0)
    P  = max(P,  0.0)
    S  = max(S,  0.0)
    N  = max(N,  0.0)
    O  = max(O,  0.0)
    Og = max(Og, 0.0)

    # Resource limitation
    phiX = Phi_X(S, N, O)
    phiP = Phi_P(S, N, O)

    # Induction switch
    # Before induction: full growth rate, no production
    # After  induction: growth reduced to beta fraction, production active
    if t < t_ind:
        mu_e = mu_B               # uninduced specific growth rate (1/h)
        q_e  = 0.0                # no nanobody synthesis yet
    else:
        mu_e = beta_B * mu_B      # post-induction reduced growth (1/h)
        q_e  = q_B                # nanobody production rate (mg/g/h)

    # ------ ODEs ------

    # dX/dt — biomass
    dX = mu_e * phiX * X

    # dP/dt — nanobody [Bacillus-specific: kprot*X*P degradation term]
    dP = q_e * phiP * X - (kd0_B + kprot_B * X) * P

    # dS/dt — glucose consumption (growth + production sinks)
    dS = -(1.0 / YX_S_B) * mu_e * phiX * X \
         -(1.0 / YP_S_B) * q_e  * phiP * X

    # dN/dt — nitrogen/amino-acid precursor pool [Bacillus-specific †]
    dN = -(1.0 / YX_N_B) * mu_e * phiX * X \
         -(1.0 / YP_N_B) * q_e  * phiP * X

    # dO/dt — dissolved oxygen (transfer in from gas, consumption out)
    dO = kLa_B * (Og - O) \
         -(1.0 / YX_O_B) * mu_e * phiX * X \
         -(1.0 / YP_O_B) * q_e  * phiP * X

    # dOg/dt — gas-phase oxygen (depletes as it transfers to liquid)
    dOg = -alpha_B * Og - kLa_B * (V_l_B / V_g_B) * (Og - O)

    return [dX, dP, dS, dN, dO, dOg]

# ------------------------------------------------------------------
# Simulation runner
# ------------------------------------------------------------------

def run_simulation(t_ind):
    """
    Integrate the ODE system from t=0 to t=T for a given induction time.

    Parameters
    ----------
    t_ind : float — induction time (h); use 0 for immediate induction

    Returns
    -------
    t_s : 1-D array, shape (1000,) — time points (h)
    s_s : 2-D array, shape (6, 1000) — state trajectories
          rows: [X, P, S, N, O, Og]
    """
    y0  = [X0_B, P0_B, S0_B, N0_B, O0_B, Og0_B]
    sol = solve_ivp(
        fun     = bacillus_ode,
        t_span  = (0.0, T),
        y0      = y0,
        args    = (t_ind,),
        method  = 'RK45',
        t_eval  = t_eval,
        rtol    = 1e-6,
        atol    = 1e-9,
        dense_output = False
    )
    return sol.t, sol.y

# ------------------------------------------------------------------
# Quick-check: run with default induction at t=10 h and plot X and P
# ------------------------------------------------------------------

if __name__ == '__main__':
    t_ind_default = 10.0   # h — change as needed

    t_s, s_s = run_simulation(t_ind_default)

    X_t, P_t = s_s[0], s_s[1]
    S_t, N_t = s_s[2], s_s[3]
    O_t, Og_t = s_s[4], s_s[5]

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    labels = ['Biomass X (g/L)', 'Nanobody P (mg/L)',
              'Glucose S (g/L)', 'Nitrogen N (g/L)',
              'Dissolved O₂ (mg/L)', 'Gas O₂ (mg/L)']
    colors = ['steelblue', 'darkorange', 'forestgreen',
              'purple', 'teal', 'crimson']
    data   = [X_t, P_t, S_t, N_t, O_t, Og_t]

    for ax, d, lab, col in zip(axes.flat, data, labels, colors):
        ax.plot(t_s, d, color=col, linewidth=2)
        ax.axvline(t_ind_default, color='red', linestyle=':', linewidth=1.2,
                   label=f'Induction at {t_ind_default} h')
        ax.set_xlabel('Time (h)')
        ax.set_ylabel(lab)
        ax.set_title(lab)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    fig.suptitle(f'Bacillus Nanobody ODE Model (t_ind = {t_ind_default} h)',
                 fontsize=14)
    fig.tight_layout()
    out_path = f'{OUT}/ode_model_overview.png'
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.show()
    print(f"Saved {out_path}")

    # Print final values
    print("\n--- Final state at t = {:.0f} h ---".format(T))
    for name, val, unit in zip(
            ['X', 'P', 'S', 'N', 'O', 'Og'],
            [X_t[-1], P_t[-1], S_t[-1], N_t[-1], O_t[-1], Og_t[-1]],
            ['g/L', 'mg/L', 'g/L', 'g/L', 'mg/L', 'mg/L']):
        print(f"  {name:3s} = {val:.4f} {unit}")
