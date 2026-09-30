# ============================================================
# fig6_final_nanobody_vs_induction.py
# Figure 6 — Final nanobody at harvest vs induction time
# ============================================================

from shared_setup import *

scan_t_ind   = np.linspace(0, 48, 97)
final_P_arr  = np.array([run_simulation(ti)[1][1, -1] for ti in scan_t_ind])

best_t_final = scan_t_ind[np.argmax(final_P_arr)]
print(f"Best induction time (final nanobody at harvest): {best_t_final:.1f} h")

fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(scan_t_ind, final_P_arr, color='steelblue', linewidth=2)
ax.axvline(best_t_final, color='red', linestyle='--', linewidth=1.4,
           label=f'Optimal $t_i={best_t_final:.1f}$ h')
ax.set_xlabel('Induction time (h)', fontsize=12)
ax.set_ylabel('Final nanobody at harvest (mg/L)', fontsize=12)
ax.set_title('Final Nanobody at Harvest vs Induction Time', fontsize=13)
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

fig.tight_layout()
fig.savefig(f'{OUT}/fig6_final_nanobody_vs_induction.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved fig6_final_nanobody_vs_induction.png")
