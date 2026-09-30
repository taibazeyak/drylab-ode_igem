# ============================================================
# fig9_multiple_induction_comparison.py
# Figure 9 — Biomass & nanobody for 5 induction strategies
# NOTE: best_t_final is hard-coded as 6.5 h from the scan.
#       Re-run fig6 first if you changed parameters.
# ============================================================

from shared_setup import *

best_t_final   = 6.5   # from fig6 scan — update if parameters change
compare_t_inds = [0.0, 12.0, 24.0, best_t_final, 40.0]
colors_cmp     = ['royalblue', 'forestgreen', 'darkorange', 'red', 'purple']
lw_cmp         = [1.5, 1.5, 1.5, 2.5, 1.5]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

for ti, col, lw in zip(compare_t_inds, colors_cmp, lw_cmp):
    t_s, s_s = run_simulation(ti)
    label = (f'$t_i={ti:.0f}$ h'
             + (' ← optimal' if np.isclose(ti, best_t_final) else ''))
    ax1.plot(t_s, s_s[0], color=col, linewidth=lw, label=label)
    ax2.plot(t_s, s_s[1], color=col, linewidth=lw, label=label)

for ax, ylabel, title in zip(
        [ax1, ax2],
        ['Biomass X (g/L)', 'Nanobody P (mg/L)'],
        ['Biomass Trajectories', 'Nanobody Trajectories']):
    ax.set_xlabel('Time (h)', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

fig.suptitle('Multiple Induction Time Comparison', fontsize=14, y=1.01)
fig.tight_layout()
fig.savefig(f'{OUT}/fig9_multiple_induction_comparison.png',
            dpi=150, bbox_inches='tight')
plt.show()
print("Saved fig9_multiple_induction_comparison.png")
