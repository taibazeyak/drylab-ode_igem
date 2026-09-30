# ============================================================
# fig5_peak_nanobody_vs_induction.py
# Figure 5 — Peak nanobody concentration vs induction time
# Runs the full induction scan (0–48 h) — takes ~30s in Colab
# ============================================================

from shared_setup import *

scan_t_ind = np.linspace(0, 48, 97)
peak_P_arr = np.array([np.max(run_simulation(ti)[1][1]) for ti in scan_t_ind])

best_t_peak = scan_t_ind[np.argmax(peak_P_arr)]
print(f"Best induction time (peak nanobody): {best_t_peak:.1f} h")

fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(scan_t_ind, peak_P_arr, color='darkorange', linewidth=2)
ax.axvline(best_t_peak, color='red', linestyle='--', linewidth=1.4,
           label=f'Optimal $t_i={best_t_peak:.1f}$ h')
ax.set_xlabel('Induction time (h)', fontsize=12)
ax.set_ylabel('Peak nanobody P (mg/L)', fontsize=12)
ax.set_title('Peak Nanobody vs Induction Time', fontsize=13)
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

fig.tight_layout()
fig.savefig(f'{OUT}/fig5_peak_nanobody_vs_induction.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved fig5_peak_nanobody_vs_induction.png")
