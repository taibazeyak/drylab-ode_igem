# ============================================================
# Produces: Yeast/yeast_constitutive_vs_inducible.png
# Model 1. Constitutive expression (t_i = 0) versus inducible
# expression at the analytically optimal induction time.
# ============================================================
from yeast_shared import *

tau_star, t_i_star = optimal_induction_time()
print("Optimal production duration tau*:", round(tau_star, 2), "h")
print("Optimal induction time t_i*:", round(t_i_star, 2), "h")

P_const = nanobody(t, 0)
P_ind   = nanobody(t, t_i_star)

plt.figure(figsize=(7, 5))
plt.plot(t, P_const, color=purple, linewidth=2.5, label="Constitutive: $t_i=0$")
plt.plot(t, P_ind,   color=green,  linewidth=2.5,
         label=f"Inducible: $t_i^*$ = {t_i_star:.1f} h")
plt.axvline(t_i_star, color=orange, linestyle="--", linewidth=2,
            label="Optimal induction")
plt.xlabel("Time (h)", fontsize=13)
plt.ylabel("Nanobody, P(t) (mg/L)", fontsize=13)
plt.xticks(fontsize=11); plt.yticks(fontsize=11)
plt.legend(fontsize=10); plt.grid(True, alpha=0.3); plt.tight_layout()
plt.savefig(f"{OUT}/yeast_constitutive_vs_inducible.png", dpi=150,
            bbox_inches="tight")
print("Saved yeast_constitutive_vs_inducible.png")
