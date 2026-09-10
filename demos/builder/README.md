# Builder demos

Quick-start tutorials that use the fluent
[`FermenterBuilder`](../../models/vlmodels/fermenter/config/builder.py)
API. Each demo configures a fermenter in a handful of chained calls
and runs a `Simulation` shipped in
[`simulation-class`](../../docs/dev/implementation/shipped/SIMULATION_CLASS.md).

Read these if you want to **run** a fermenter — you don't care how
the `Phase` / `Link` / `ControlVolume` plumbing fits together, just
that it does. For that, read the
[`model_api/D2Cworkshop/`](../model_api/D2Cworkshop/) demos.

## Demos

| File | Topology / features |
|---|---|
| [batch_fermenter.py](batch_fermenter.py) | 0-D batch, sparged. PI pH control. The simplest topology — read first. |
| [cstr_fermenter.py](cstr_fermenter.py) | 0-D chemostat (feed + drain, constant volume). PI pH + DO-driven agitation control. |
| [fed_batch_fermenter.py](fed_batch_fermenter.py) | 0-D fed-batch (feed only; working volume grows). PI pH + DO agitation. Includes a substrate-balance sanity check. |
| [microplate_fermenter.py](microplate_fermenter.py) | 96-well plate well, membrane gas exchange (no sparging). Minimal chemistry — useful for seeing what falls out when you strip most of the moving parts. |

## Running

From the repo root:

```bash
python demos/builder/batch_fermenter.py
python demos/builder/cstr_fermenter.py
python demos/builder/fed_batch_fermenter.py
python demos/builder/microplate_fermenter.py
```

Each prints a short report (initial → final concentrations, pH where
applicable, controller diagnostics).

## Reading order

1. **`batch_fermenter.py`** — the canonical build → wire → run →
   report pattern under the simplest topology.
2. **`cstr_fermenter.py`** — adds inlet/outlet boundaries
   (`LiquidFeed` + `LiquidDrain`) and a DO controller alongside the
   pH one. Same chemistry as batch.
3. **`fed_batch_fermenter.py`** — adds the volume-grows-over-time
   wrinkle and the substrate-balance sanity check.
4. **`microplate_fermenter.py`** — strips chemistry to the bare
   minimum and swaps sparging for a `MembraneGasBoundary`.

## Canonical shape

```python
from vlmodels.fermenter.config import FermenterBuilder
from PyOMES.core import Simulation
from PyOMES.control.cv_loops import PHController

cv = (
    FermenterBuilder()
    .vessel(V_total_L=2000, T_K=305.15)
    .gas_feed(vvm_min=1.0, composition={"O2": 0.21, "N2": 0.79})
    .transfer_kinetic(kLa_O2=150.0)
    .chemistry(speciation_level=1)
    .organism("Yeast")
    .substrate("AceticAcid", mu_max=0.5, Ks=5e-3, yield_gX_gS=0.36)
    .build()
)
cv.boundaries.append(PressureReliefVent(P_set_atm=1.10, mode="instant"))

sim = Simulation(cvs={"main": cv}, controllers=[PHController(setpoint=5.0)])
result = sim.run(tau_h=5.0, n_steps=1000)
```

`FermenterBuilder.build_simulation_and_run(tau_h, n_steps)` is the
one-shot equivalent of `build()` → boundaries → `Simulation(...).run(...)`
when you don't need to inspect the intermediate CV or attach custom
boundaries. `batch_fermenter.py` uses that form.
