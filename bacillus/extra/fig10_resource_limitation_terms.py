# ============================================================
# fig10_resource_limitation_terms.py
# Figure 10 — Phi_X(t) and Phi_P(t) resource limitation terms
# Uses best_t_final = 6.5 h (update if parameters change)
# ============================================================

from shared_setup import *

best_t_final = 6.5   # from fig6 scan — update if parameters change
t_s, s_s     = run_simulation(best_t_final)

phiX_t = np.array([Phi_X(s_s[2, i], s_s[3, i], s_s[4, i]) for i in range(len(t_s))])
phiP_t = np.array([Phi_P(s_s[2, i], s_s[3, i], s_s[4, i]) for i in range(len(t_s))])

fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(t_s, phiX_t, color='steelblue',  linewidth=2,
        label='$\\Phi_X$ — growth limitation')
ax.plot(t_s, phiP_t, color='darkorange', linewidth=2, linestyle='--',
        label='$\\Phi_P$ — production limitation')
ax.axvline(best_t_final, color='red', linestyle=':', linewidth=1.4,
           label=f'Induction at $t_i={best_t_final:.1f}$ h')
ax.set_xlabel('Time (h)', fontsize=12)
ax.set_ylabel('Resource limitation factor (–)', fontsize=12)
ax.set_ylim(-0.05, 1.10)
ax.set_title('Resource Limitation Terms $\\Phi_X$ and $\\Phi_P$', fontsize=13)
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

fig.tight_layout()
fig.savefig(f'{OUT}/fig10_resource_limitation_terms.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved fig10_resource_limitation_terms.png")
