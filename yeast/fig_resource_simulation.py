# ============================================================
# Produces: Yeast/yeast_resource_simulation.png
# Model 2, five-state resource-depleting run at t_ind = 10 h.
# All states normalised to their own maximum.
# ============================================================
from yeast_shared import *

t_out, sol = run_5ode(t_ind_default)
X, P, S, O, Og = sol

def norm_(y):
    ymax = np.max(y)
    return y / ymax if ymax > 0 else y

plt.figure(figsize=(8, 5))
plt.plot(t_out, norm_(X),  color=blue,   linewidth=2.5, label="Biomass")
plt.plot(t_out, norm_(P),  color=red,    linewidth=2.5, label="Nanobody")
plt.plot(t_out, norm_(S),  color=green,  linewidth=2.5, label="Sugar")
plt.plot(t_out, norm_(O),  color=orange, linewidth=2.5, label="Dissolved oxygen")
plt.plot(t_out, norm_(Og), color=purple, linewidth=2.5, label="Headspace oxygen")
plt.axvline(t_ind_default, color=black, linestyle="--", linewidth=1.8,
            label="Induction")
plt.xlabel("Time (h)", fontsize=13)
plt.ylabel("Normalized value", fontsize=13)
plt.xticks(fontsize=11); plt.yticks(fontsize=11)
plt.legend(fontsize=10); plt.grid(True, alpha=0.3); plt.tight_layout()
plt.savefig(f"{OUT}/yeast_resource_simulation.png", dpi=150,
            bbox_inches="tight")

print("Final biomass:",          X[-1])
print("Final nanobody:",         P[-1])
print("Final sugar:",            S[-1])
print("Final dissolved oxygen:", O[-1])
print("Final headspace oxygen:", Og[-1])
print("Saved yeast_resource_simulation.png")
