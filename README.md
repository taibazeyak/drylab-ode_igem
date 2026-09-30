# iGEM KAUST 2026 — Dry Lab: ODE Models of Nanobody Production

Dynamic (time-resolved) models of recombinant nanobody production in
*Saccharomyces cerevisiae* and *Bacillus subtilis*, built to answer a question
that flux balance analysis cannot: **when should recombinant expression be
switched on?**

FBA is steady-state by construction. It reports rates, not trajectories, and
has no notion of a batch that starts, depletes its medium and is harvested.
These models supply that time dimension.

---

## The three models

| | Model 1 | Model 2 | Model 3 |
|---|---|---|---|
| host | *S. cerevisiae* | *S. cerevisiae* | *B. subtilis* |
| type | analytical, closed form | 5-state ODE | 6-state ODE |
| states | X, P | X, P, S, O, O_g | X, P, S, **N**, O, O_g |
| resources | none (unlimited) | sugar, oxygen | sugar, **nitrogen**, oxygen |
| product loss | first-order `kd·P` | first-order `kd·P` | `(kd0 + kprot·X)·P` |
| cultivation time | 48 h | 48 h | 72 h |
| burden factor β | 0.50 | 0.50 | 0.10 |

**Model 1** assumes resources are never limiting, which makes it solvable in
closed form. The optimal production-phase duration is

```
τ* = 1/(β·μ + kd) · ln( (μ + kd) / ((1 − β)·μ) )      t_i* = T − τ*
```

and the optimal induction time follows. Three limiting cases (β→1, β→0, kd→0)
each recover a biologically sensible answer, which is the validation available
to a purely analytical model.

**Model 2** adds a depletable sugar pool, dissolved oxygen and a headspace
oxygen reservoir, each entering through multiplicative Monod terms. Once the
batch can run out of things, the recommendation inverts: the ideal model
favours late induction, the resource-limited one favours early induction.

**Model 3** adds two features specific to *B. subtilis*: an explicit
nitrogen/amino-acid precursor pool (a nanobody is a protein, so its synthesis
needs more than carbon and energy), and a biomass-proportional protease loss
term, since *B. subtilis* secretes extracellular proteases that degrade
heterologous product. That second term means peak titre and harvest titre are
not the same quantity.

---

## Repository layout

```
yeast/
  yeast_shared.py                      the model: Models 1 and 2
  fig_constitutive_vs_inducible.py     constitutive vs inducible expression
  fig_model1_induction_times.py        biomass + nanobody, 5 induction times
  fig_resource_simulation.py           all 5 states, t_ind = 10 h
  fig_induction_scan_and_heatmap.py    induction sweep + P(t) landscape

bacillus/
  shared_setup.py                      the model: Model 3
  fig2_sugar_nitrogen.py               carbon and nitrogen dynamics
  fig6_final_nanobody_vs_induction.py  final titre vs induction time
  fig8_heatmap_nanobody.py             P(t) landscape
  fig9_multiple_induction_comparison.py  biomass + nanobody, 5 strategies
  extra/                               further figures not used on the Wiki
```

The two `*_shared*.py` / `shared_setup.py` files contain the parameters, the
ODE systems and the solvers. Everything else imports them and plots.

---

## Running

```bash
pip install -r requirements.txt

cd yeast    && python fig_model1_induction_times.py
cd bacillus && python fig9_multiple_induction_comparison.py
```

Each script writes its figure into an `outputs/` folder beside it.
`yeast/fig_induction_scan_and_heatmap.py` runs 500 integrations and takes
about a minute; the rest are instant.

Solver: `scipy.integrate.solve_ivp`, RK45.
Yeast uses `rtol=1e-7, atol=1e-9` for single trajectories, relaxed to
`rtol=1e-6, atol=1e-8` for the 500-run sweep. Bacillus uses
`rtol=1e-6, atol=1e-9`. All state variables are floored at zero at each
evaluation so the Monod terms stay well defined.

---

## Expected output

If the code is working, these are the numbers it should print:

| | value |
|---|---|
| Model 1, optimal production duration τ* | 15.25 h |
| Model 1, optimal induction time t_i* | 32.75 h |
| Model 2 at t_ind = 10 h, final biomass | 0.4404 g/L |
| Model 2 at t_ind = 10 h, final nanobody | 0.2020 mg/L |
| Model 2, peak nanobody across the sweep | 0.2383 mg/L at t_i = 0 h |
| Model 3, optimal induction time | 6.5 h |
| Model 3, final titre at that optimum | 198.0 mg/L |

---

## A note on the parameters

Only the yeast maximum specific growth rate (μ = 0.0854 h⁻¹) is
FBA/GEM-informed. Every other numerical value is an illustrative simulation
setting: no parameter here was fitted to experimental data.

The *structure* is literature-supported (Monod kinetics for nutrient-limited
growth; a k_L·a term for oxygen transfer; protease-mediated product loss as a
known *B. subtilis* bottleneck), but the numbers are not strain-specific.

These models should therefore be read as **qualitative demonstrations of the
induction-timing trade-off**, not as calibrated predictions. The optimal
induction times above are model-dependent and should be refined against
experimental data or FBA-derived estimates before being used as process
set-points.

---

## License

<!-- iGEM requires an open license. MIT or CC-BY-4.0 are both common.
     Add a LICENSE file before submission. -->

---

*iGEM KAUST 2026 — Dry Lab. See the team Wiki for the full write-up, including
the genome-scale metabolic modelling that these dynamic models sit alongside.*
