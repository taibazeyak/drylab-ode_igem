# ============================================================
# fig3_oxygen.py
# Figure 3 — Dissolved O2 O(t) and Headspace O2 Og(t)
# ============================================================

from shared_setup import *

t_ind_example = 10.0
t_out, sol = run_simulation(t_ind_example)

fig, ax = plt.subplots(figsize=(9, 5))

ax.plot(t_out, sol[4], color='deepskyblue', linewidth=2,
        label='Dissolved O₂ (mg/L)')
ax.plot(t_out, sol[5], color='navy',        linewidth=2, linestyle='--',
        label='Headspace O₂ (mg/L)')
vline(ax, t_ind_example)
ax.set_xlabel('Time (h)', fontsize=12)
ax.set_ylabel('Oxygen Concentration (mg/L)', fontsize=12)
ax.set_title('Oxygen Dynamics', fontsize=13)
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

fig.tight_layout()
fig.savefig(f'{OUT}/fig3_oxygen.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved fig3_oxygen.png")
