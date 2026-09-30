# ============================================================
# fig7_percent_yield_lost.py
# Figure 7 — Percent nanobody yield lost from peak to harvest
# ============================================================

from shared_setup import *

scan_t_ind = np.linspace(0, 48, 97)
peak_P_arr  = []
final_P_arr = []

for ti in scan_t_ind:
    _, s = run_simulation(ti)
    peak_P_arr.append(np.max(s[1]))
    final_P_arr.append(s[1, -1])

peak_P_arr  = np.array(peak_P_arr)
final_P_arr = np.array(final_P_arr)

with np.errstate(divide='ignore', invalid='ignore'):
    pct_lost = np.where(peak_P_arr > 0,
                        (peak_P_arr - final_P_arr) / peak_P_arr * 100,
                        0.0)

fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(scan_t_ind, pct_lost, color='crimson', linewidth=2)
ax.set_xlabel('Induction time (h)', fontsize=12)
ax.set_ylabel('Yield lost from peak to harvest (%)', fontsize=12)
ax.set_title('Percent Nanobody Lost (Peak → Harvest)', fontsize=13)
ax.grid(alpha=0.3)

fig.tight_layout()
fig.savefig(f'{OUT}/fig7_percent_yield_lost.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved fig7_percent_yield_lost.png")
