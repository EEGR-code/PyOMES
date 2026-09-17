# Protocols

Single-class, single-feature deep dives: you already know which class
you're using and want its instantiation/call conventions and gotchas,
rather than a scenario-first worked example. Two independent protocol
families, each its own subfolder:

| Subfolder | Covers |
|---|---|
| [`ChemicalEquilibriumProtocol/`](ChemicalEquilibriumProtocol/) | `ChemicalEquilibriumEngineProtocol` — the three peer engines (`BisectionChemicalEquilibriumEngine`, `NRChemicalEquilibriumEngine`, `PHREEQCChemicalEquilibriumEngine`) that solve aqueous equilibrium chemistry: construction, `solve()` conventions, and the black/gray/white-box capability tiers. Start at [`0_README.ipynb`](ChemicalEquilibriumProtocol/0_README.ipynb) for the architecture overview and engine comparison table. |
| [`SolverProtocols/`](SolverProtocols/) | `StepSolver`/`SystemSolver` — the two orthogonal solver axes (per-CV physics vs. whole-system orchestration) and how to write your own. Start at [`0_README.ipynb`](SolverProtocols/0_README.ipynb). |

For the numerical *accuracy* side of the equilibrium engines (does the NR
solution actually match a closed-form formula, or PHREEQC?) see
[`tests/validation/speciation/`](../../../tests/validation/speciation/)
instead — that's benchmark/validation content, not mechanics.

## Running

```bash
jupyter lab docs/tutorials/protocols/
```
