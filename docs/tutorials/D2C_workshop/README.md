# D2C Workshop

Manual assembly of a `Simulation` from `Phase` + `Link` + `ControlVolume`
primitives — no builder. An engineering-first approach: derive the design
basis analytically first, then build the simulation and check it against
that prediction. The three systems form the DARE2CYCLE scale-up
progression: sealed microplate well → sparged batch fermenter → continuous
CSTR.

## Notebooks

| Notebook | System | Central question | Learning objectives |
|---|---|---|---|
| [Example1_mtp_well.ipynb](Example1_mtp_well.ipynb) | Sealed MTP well | Does the headspace hold enough O₂ to sustain growth before the substrate runs out? | Define species and reactions (`Species`, `ReactionBuilder`, `ReactionSystem`); select a gas-liquid transfer model (`EquilibriumTransferModel`); build a closed `ControlVolume` with no external boundaries |
| [Example2_batch_fermenter.ipynb](Example2_batch_fermenter.ipynb) | Sparged batch fermenter | Does aeration meet oxygen demand, and is pH control necessary? | Switch to kinetic gas transfer (`KineticTransferModel`); attach gas inlet and pressure vent boundaries (`GasFeed`, `PressureReliefVent`); add a pH controller (`PHController`) |
| [Example3_CSTR.ipynb](Example3_CSTR.ipynb) | CSTR | How should the dilution rate be chosen to balance productivity against washout risk? | Add liquid feed and drain boundaries (`LiquidFeed`, `LiquidDrain`) for continuous operation; configure steady-state operation via dilution rate; validate convergence to the analytical steady state |

Launch from the repo root with `jupyter lab docs/tutorials/D2C_workshop/`.

## `raw_construction.py`

A separate, standalone comparison script — mirrors
[`../templates/batch_fermenter.py`](../templates/batch_fermenter.py)'s
topology (0-D sparged batch, aerobic growth on acetic acid, PI pH control)
built explicitly from `PyOMES.core` primitives instead of
`StirredTankBuilder`, to show what the builder hides. Not used by the three
notebooks above. Its chemistry is declared inline rather than imported from
[`../reactions/reaction_system.ipynb`](../reactions/reaction_system.ipynb), so
this folder stays self-sufficient.

```bash
python docs/tutorials/D2C_workshop/raw_construction.py
```
