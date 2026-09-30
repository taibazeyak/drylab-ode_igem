# ============================================================
# fig2_sugar_nitrogen.py
# Figure 2 — Sugar S(t) and Nitrogen N(t) for t_ind = 10 h
# ============================================================

from shared_setup import *

t_ind_example = 10.0
t_out, sol = run_simulation(t_ind_example)

fig, ax = plt.subplots(figsize=(9, 5))

ax.plot(t_out, sol[2], color='#1f77b4', linewidth=2, label='Sugar S (g/L)')
ax.plot(t_out, sol[3], color='#ff7f0e', linewidth=2,
        label='Nitrogen / AA pool N (g/L)')
vline(ax, t_ind_example)
ax.set_xlabel('Time (h)', fontsize=12)
ax.set_ylabel('Concentration (g/L)', fontsize=12)
ax.set_title('Carbon and Nitrogen Resources', fontsize=13)
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

fig.tight_layout()
fig.savefig(f'{OUT}/fig2_sugar_nitrogen.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved fig2_sugar_nitrogen.png")
