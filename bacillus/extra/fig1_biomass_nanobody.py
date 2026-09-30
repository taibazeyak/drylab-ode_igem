# ============================================================
# fig1_biomass_nanobody.py
# Figure 1 — Biomass X(t) and Nanobody P(t) for t_ind = 10 h
# ============================================================

from shared_setup import *

t_ind_example = 10.0
t_out, sol = run_simulation(t_ind_example)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(t_out, sol[0], color='steelblue', linewidth=2)
vline(ax1, t_ind_example)
ax1.set_xlabel('Time (h)', fontsize=12)
ax1.set_ylabel('Biomass X (g/L)', fontsize=12)
ax1.set_title('Bacillus Biomass', fontsize=13)
ax1.legend(fontsize=10)
ax1.grid(alpha=0.3)

ax2.plot(t_out, sol[1], color='darkorange', linewidth=2)
vline(ax2, t_ind_example)
ax2.set_xlabel('Time (h)', fontsize=12)
ax2.set_ylabel('Nanobody P (mg/L)', fontsize=12)
ax2.set_title('Nanobody Concentration', fontsize=13)
ax2.legend(fontsize=10)
ax2.grid(alpha=0.3)

fig.suptitle('Bacillus-Adapted Model — Biomass & Nanobody', fontsize=14, y=1.01)
fig.tight_layout()
fig.savefig(f'{OUT}/fig1_biomass_nanobody.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved fig1_biomass_nanobody.png")
