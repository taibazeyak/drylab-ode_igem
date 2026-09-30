# ============================================================
# fig8_heatmap_nanobody.py
# Figure 8 — Heatmap of P(t) across induction times
# ============================================================

from shared_setup import *

scan_t_ind = np.linspace(0, 48, 97)
P_heatmap  = np.zeros((len(scan_t_ind), len(t_eval)))

for idx, ti in enumerate(scan_t_ind):
    _, s = run_simulation(ti)
    P_heatmap[idx, :] = s[1]

global_max = np.max(P_heatmap)

fig, ax = plt.subplots(figsize=(11, 6))
im = ax.pcolormesh(t_eval, scan_t_ind, P_heatmap,
                   shading='auto', cmap='viridis')
cbar = fig.colorbar(im, ax=ax)
cbar.set_label('Nanobody P(t) (mg/L)', fontsize=11)

# 80% contour
try:
    cs = ax.contour(t_eval, scan_t_ind, P_heatmap,
                    levels=[0.80 * global_max],
                    colors='white', linewidths=1.5)
except Exception:
    pass

ax.set_xlabel('Cultivation time (h)', fontsize=12)
ax.set_ylabel('Induction time (h)', fontsize=12)

fig.tight_layout()
fig.savefig(f'{OUT}/fig8_heatmap_nanobody.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved fig8_heatmap_nanobody.png")
