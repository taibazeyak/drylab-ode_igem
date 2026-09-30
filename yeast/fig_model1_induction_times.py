# ============================================================
# Produces: Yeast/yeast_model1_induction_times.png
# Model 1. Biomass (left) and nanobody (right) for five
# induction times, including the analytical optimum.
# ============================================================
from yeast_shared import *

tau_star, t_i_star = optimal_induction_time()
induction_times = [0, 12, 24, t_i_star, 40]

cmap = plt.get_cmap("plasma")
norm = plt.Normalize(0, T)
hours_ticks = [0, 12, 24, 36, 48]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5))

for t_i in reversed(induction_times):
    color      = cmap(norm(t_i))
    is_optimal = (t_i == t_i_star)
    linestyle  = "--" if is_optimal else "-"
    label = (f"$t_i^*$ = {t_i:.1f} h" if is_optimal
             else f"$t_i$ = {t_i:.1f} h")
    ax1.plot(t, biomass(t, t_i),  color=color, linestyle=linestyle,
             linewidth=3, label=label)
    ax2.plot(t, nanobody(t, t_i), color=color, linestyle=linestyle,
             linewidth=3)

ax1.set_xlabel("Time (h)", fontsize=12)
ax1.set_ylabel("Biomass, X(t) (g/L)", fontsize=12)
ax1.set_xlim(0, T); ax1.set_xticks(hours_ticks); ax1.grid(True, alpha=0.3)

ax2.set_xlabel("Time (h)", fontsize=12)
ax2.set_ylabel("Nanobody, P(t) (mg/L)", fontsize=12)
ax2.set_xlim(0, T); ax2.set_xticks(hours_ticks); ax2.grid(True, alpha=0.3)

handles, labels = ax1.get_legend_handles_labels()
handles, labels = handles[::-1], labels[::-1]
plt.tight_layout()
fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1, 0.5),
           fontsize=11, frameon=True)
fig.savefig(f"{OUT}/yeast_model1_induction_times.png", dpi=150,
            bbox_inches="tight")
print("Saved yeast_model1_induction_times.png")
