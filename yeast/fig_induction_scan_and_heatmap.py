# ============================================================
# Produces: Yeast/yeast_induction_scan.png
#           Yeast/yeast_yield_loss_heatmap.png
# Model 2. Sweeps induction time across the whole batch and
# records peak nanobody, yield lost by waiting, and the full
# P(t) landscape.  Takes ~1 min (500 integrations).
# ============================================================
from yeast_shared import *

induction_times = np.linspace(0, T, 500)
t_eval = np.linspace(0, T, 600)

max_P_list, final_P_list, P_matrix = [], [], []
for t_ind in induction_times:
    _, y = run_5ode(t_ind, rtol=1e-6, atol=1e-8)   # relaxed for the sweep
    P_traj = y[1]
    max_P_list.append(np.max(P_traj))
    final_P_list.append(P_traj[-1])
    P_matrix.append(P_traj)

max_P   = np.array(max_P_list)
final_P = np.array(final_P_list)
P_matrix = np.array(P_matrix).T
hours_ticks = [0, 12, 24, 36, 48]

# --- Peak nanobody vs induction time ---------------------------------
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(induction_times, max_P, color=red, linewidth=2.5)
ax.set_xlabel("Induction Time (h)", fontsize=12)
ax.set_ylabel("Maximum Nanobody P", fontsize=12)
ax.set_title("Peak Nanobody Yield vs. Induction Time", fontsize=13)
ax.set_xlim(0, T); ax.set_xticks(hours_ticks); ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(f"{OUT}/yeast_induction_scan.png", dpi=150, bbox_inches="tight")

# --- Percent yield lost by harvesting at the end ----------------------
percent_lost = ((max_P - final_P) / max_P) * 100
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(induction_times, percent_lost, color=red, linewidth=2.5)
ax.set_xlabel("Induction Time (h)", fontsize=12)
ax.set_ylabel("Yield Lost (%)", fontsize=12)
ax.set_title("Percent of Nanobody Lost by Waiting Until End", fontsize=13)
ax.set_xlim(0, T); ax.set_xticks(hours_ticks); ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(f"{OUT}/yeast_percent_yield_lost.png", dpi=150, bbox_inches="tight")

# --- P(t) landscape ---------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5.5))
X_grid, Y_grid = np.meshgrid(induction_times, t_eval)
cp = ax.pcolormesh(X_grid, Y_grid, P_matrix, shading="auto", cmap="viridis")
ax.contour(X_grid, Y_grid, P_matrix, levels=[0.2], linestyles="--",
           colors="black")
colorbar = fig.colorbar(cp)
colorbar.set_label("Nanobody Concentration P(t)", fontsize=11)
ax.set_xlabel("Induction Time (h)", fontsize=12)
ax.set_ylabel("Time (h)", fontsize=12)
ax.set_title("Nanobody P(t) over Varying Induction Times", fontsize=13)
ax.set_xticks(hours_ticks); ax.set_yticks(hours_ticks)
ax.set_aspect("equal", "box")
fig.tight_layout()
fig.savefig(f"{OUT}/yeast_yield_loss_heatmap.png", dpi=150, bbox_inches="tight")

print(f"Peak nanobody max = {max_P.max():.4f} mg/L "
      f"at t_i = {induction_times[np.argmax(max_P)]:.1f} h")
print("Saved yeast_induction_scan.png, yeast_percent_yield_lost.png, "
      "yeast_yield_loss_heatmap.png")
