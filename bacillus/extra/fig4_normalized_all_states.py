# ============================================================
# fig4_normalized_all_states.py
# Figure 4 — All 6 state variables normalised to their max
# ============================================================

from shared_setup import *

t_ind_example = 10.0
t_out, sol = run_simulation(t_ind_example)

labels  = ['X (biomass)', 'P (nanobody)', 'S (sugar)',
           'N (nitrogen)', 'O (dissolved O₂)', 'Og (headspace O₂)']
colors  = ['steelblue', 'darkorange', 'forestgreen',
           'mediumpurple', 'deepskyblue', 'navy']
linests = ['-', '-', '-', '-', '--', '--']

fig, ax = plt.subplots(figsize=(10, 5))

for i in range(6):
    mx = np.max(sol[i])
    if mx > 0:
        ax.plot(t_out, sol[i] / mx, color=colors[i],
                linestyle=linests[i], linewidth=1.8, label=labels[i])

vline(ax, t_ind_example)
ax.set_xlabel('Time (h)', fontsize=12)
ax.set_ylabel('Normalised value (–)', fontsize=12)
ax.set_title('All State Variables — Normalised', fontsize=13)
ax.legend(fontsize=9, ncol=2)
ax.grid(alpha=0.3)

fig.tight_layout()
fig.savefig(f'{OUT}/fig4_normalized_all_states.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved fig4_normalized_all_states.png")
