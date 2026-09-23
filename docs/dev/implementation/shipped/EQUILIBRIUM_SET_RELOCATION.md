# EquilibriumSet Relocation — Design Note

> **Status: Shipped 2026-09-23** — merged into `main` via `git merge --no-ff`
> as commit `852bff2`, tagged `equilibrium-set-relocation-shipped`. Two
> checkpoints: the move of `EquilibriumSet`/`EquilibriumDef`/`WaterDef` to
> `chemical_equilibrium/engines/bisection/equilibria.py`, then the deletion of
> `EquilibriumSet.bsm2_default()`. All three open questions were settled as the
> note leaned (no re-export, the `from_reactions()` import stays lazy, one
> checklist with two checkpoints). The note below is the audit as originally
> written; where it disagrees with the checklist's "Re-verification" section
> (it listed three consumer files, and missed comment and docstring references
> in `acid_base.py`, `bsm2.py`, `test_speciation.py` and `OPEN_WORK.md`), the
> checklist is correct. Full suite green post-merge: 2086 passed, 0 failed. See
> [`EQUILIBRIUM_SET_RELOCATION_CHECKLIST.md`](EQUILIBRIUM_SET_RELOCATION_CHECKLIST.md).

> Written before implementation, from a conversational audit. Surfaced 2026-09-23 during a conversational investigation of
> `chemistry/species_check.py`'s package placement (a separate, already
> resolved question — see that conversation's conclusion: stays in
> `chemistry/`, merges into `species.py`). While auditing `chemistry/`'s
> layering more broadly, `chemistry/equilibria.py` turned out to already
> be a known but unactioned finding: `OPEN_WORK.md`'s "Three small loose
> ends" section (lines 681-684, logged 2026-09-21/22 during
> `chemistry-reactions-kinetics-cleanup`) names this exact relocation.
> This note develops that one-line finding into a concrete plan, plus a
> second finding surfaced in the same conversation while re-reading the
> file's contents: `EquilibriumSet.bsm2_default()` is vestigial and
> should be deleted, not just relocated — see "Also found" below.

## Motivation

`chemistry/equilibria.py` defines `EquilibriumSet`/`EquilibriumDef`/
`WaterDef` — a declarative charge-balance format: named systems with a
`category` (`"acid"`, `"cation_acid"`, `"inorganic_acid"`,
`"strong_ion"`), an ordered pKa ladder, an `n_active` cutoff, van't Hoff
temperature correction, and (for the VFA rows in `bsm2_default()`) a
deprecated `{name}_HA`/`{name}_A-` synthesized-key fallback for entries
with no `species_refs`. Despite living in `chemistry/` — the subpackage
whose other contents (`species.py`, `common_species.py`, `partition.py`)
are consumed across the whole engine — this format is read by exactly
one consumer: the Bisection chemical-equilibrium engine.

## Verified: consumer audit

Repo-wide search for `EquilibriumSet`/`EquilibriumDef` (all file types)
returns exactly three files:

- `PyOMES/chemistry/equilibria.py` — the definition itself.
- `PyOMES/chemical_equilibrium/engines/bisection/engine.py` — imports it
  lazily inside `from_reactions()` (`engine.py:202`:
  `from ....chemistry.equilibria import EquilibriumSet`) to build the
  set from a `ReactionSystem`'s declared `EquilibriumConstraint`s, and
  cross-references it in three docstrings (lines 18, 92, 156).
- `PyOMES/chemical_equilibrium/engines/bisection/acid_base.py` — its
  `solve_from_equilibrium_set()` takes `equilibrium_set` as a duck-typed
  parameter (no import needed; the docstring names `EquilibriumSet`/
  `EquilibriumDef` in prose only, line 847-863).

Nothing in `engines/nr/` or `engines/phreeqc/` touches it.
`chemical_equilibrium/protocols.py` — the engine-agnostic shared
interface both engines implement — has zero references to it either,
confirming it never entered the cross-engine contract. `nr/engine.py`'s
own docstring frames itself as the generalized successor ("handles
**arbitrary reaction networks**... once the chemistry is declared") and
sources its thermodynamics from `thermo/equilibrium_constants.py`'s
`vant_hoff_log_K` — a structurally different, independent path.

`chemistry/__init__.py` does not re-export `EquilibriumSet`/
`EquilibriumDef` at all (checked directly — no match). So every
consumer already reaches it via the deep path
`PyOMES.chemistry.equilibria.EquilibriumSet`, not a package-level name;
relocating it doesn't touch any re-export surface, only that deep path.

The one code call site outside `chemical_equilibrium/` that touches it
directly is `tests/standalone/test_cv_compute_interface.py:278,294`
(`from PyOMES.chemistry.equilibria import EquilibriumSet`, in
`test_bsm2_equilibrium_set_includes_carbonate_and_nh` and
`test_bsm2_vfa_entries_not_included` — both construct an
`EquilibriumSet.bsm2_default()` directly and hand it to a
`BisectionChemicalEquilibriumEngine` to check `algebraic_species()`).
Every other repo hit for `PyOMES.chemistry import *Equilibrium*` found
during the audit was `HenryEquilibrium`/`RaoultEquilibrium`/
`KspEquilibrium` from the unrelated `chemistry/partition.py` — a
different file, not part of this note.

## Also found: `EquilibriumSet.bsm2_default()` is vestigial

Independent of where the file lives, `bsm2_default()` itself (the
Rosen & Jeppsson 2006 BSM2/ADM1 acid-base preset —
`equilibria.py:432-463`) doesn't belong baked into what is otherwise a
generic charge-balance *format* module — and unlike the relocation
question, this one isn't just a layering nit, it's dead weight with a
live duplication risk.

**It's not the actual source of BSM2's numbers.**
`models/vlmodels/adm1/bsm2.py:467-739` independently re-declares the
identical pKa/van't-Hoff constants (`_BSM2_PKW`, `_BSM2_DH_W`,
`_BSM2_PKA_CO2_1`, `_BSM2_DH_CO2_1`, ...) as its own module-level
literals feeding `EquilibriumReaction` declarations — with a comment
at `bsm2.py:472,734` stating *"Matches EquilibriumSet.bsm2_default()
bit-for-bit."* That's two independent copies of the same domain
constants: one live (feeds the actual BSM2 model and the
`test_bsm2_reference.py` golden regression), one dormant
(`bsm2_default()`), with nothing but that comment enforcing they stay
in sync.

**History confirms this is leftover, not intentional.** The shipped
`CHEMISTRY_UNIFICATION_PLAN.md`/`CHEMISTRY_UNIFICATION.md` already
call it *"the legacy `bsm2_default()` factory path"* —
`chemistry-unification-3`/`3b` migrated the real BSM2 model off it
onto declared reactions. Its three sibling presets
(`bsm2_diprotic_co2`, `bsm2_with_sulfide`, `adm1_full`) were already
deleted outright in `chemistry-reactions-kinetics-cleanup` checkpoint
8 for having "zero external callers." `bsm2_default()` only survived
that sweep because it still has callers — but:

**Repo-wide search for actual invocations (`bsm2_default()`, not just
mentions) finds exactly two:** its own doctest (`equilibria.py:23`)
and `tests/standalone/test_cv_compute_interface.py:279,295`
(`test_bsm2_equilibrium_set_includes_carbonate_and_nh`,
`test_bsm2_vfa_entries_not_included`) — both checking that
`algebraic_species()` includes entries with `species_refs` and
excludes VFA rows without them. Neither test asserts anything
BSM2-specific (pH, concentrations); they exercise generic
`species_refs` mechanics and just use BSM2's ladder as convenient
multi-entry sample data. The same test file already has a smaller,
non-branded synthetic fixture for the identical purpose —
`_make_speciation_reactions()` (water + CO2 + acetate,
`test_cv_compute_interface.py:55-79`) — used by adjacent tests in the
same class. A BSM2-branded preset isn't functionally necessary for
what these two tests use it for.

**Proposed fix:** delete `bsm2_default()`, matching how its three
siblings were already deleted, and replace the two
`test_cv_compute_interface.py` call sites with a small synthetic
`EquilibriumSet` built inline — enough entries with and without
`species_refs` to exercise both assertions, without reaching for real
BSM2 numbers. This also removes the silent-drift risk between
`bsm2.py`'s live constants and `equilibria.py`'s dormant copy.

## Proposed change

Move `chemistry/equilibria.py` to
`chemical_equilibrium/engines/bisection/equilibria.py`, alongside its
one real consumer.

**Inside the moved file — two import rewrites:**
- `from ..units import R_J_PER_MOL_K as _R_J` →
  `from ....units import R_J_PER_MOL_K as _R_J` (matches the depth
  `acid_base.py`/`engine.py` already use for `from ....thermo import ...`).
- `from .common_species import CO2, HCO3_minus, CO3_2minus, NH4_plus, NH3`
  → `from ....chemistry.common_species import CO2, HCO3_minus,
  CO3_2minus, NH4_plus, NH3`.

**Consumer updates:**
1. `engine.py:202` — `from ....chemistry.equilibria import EquilibriumSet`
   → `from .equilibria import EquilibriumSet`. Now a same-package
   sibling import; see Open question 2 for whether to also hoist it out
   of the lazy `from_reactions()` body.
2. `engine.py` docstrings (lines 18, 92, 156) — repoint from
   `PyOMES.chemistry.equilibria.EquilibriumSet` to
   `PyOMES.chemical_equilibrium.engines.bisection.equilibria.EquilibriumSet`.
3. `acid_base.py` — no import change (duck-typed parameter); prose
   mentions are unqualified already.
4. `tests/standalone/test_cv_compute_interface.py:278,294` — update the
   import path the same way.
5. `chemistry/__init__.py` — no export change (nothing to remove — see
   audit above), but its module docstring ("species declarations,
   phase-partition models, and **acid-base equilibrium sets**") should
   drop that last clause, since it would no longer be true.
6. `docs/architecture.md` — move `equilibria.py` out of the `chemistry/`
   file listing (line 357) into the `bisection/` listing (line 386).
7. `OPEN_WORK.md:681-684` — close out this entry once shipped (move to
   the shipping phase's checklist, per this folder's usual convention).
8. Delete `EquilibriumSet.bsm2_default()` (`equilibria.py:432-463`) and
   its doctest reference (`equilibria.py:23`) — see "Also found" above.
9. `test_cv_compute_interface.py:279,295` — replace the two
   `EquilibriumSet.bsm2_default()` calls with an inline synthetic
   `EquilibriumSet` (entries with and without `species_refs`, matching
   what each test actually asserts).

**Left alone on purpose:** the shipped-phase checklists
(`CHEMISTRY_REACTIONS_KINETICS_CLEANUP_CHECKLIST.md`,
`CHEMICAL_EQUILIBRIUM_ENGINES_SUBFOLDER_CHECKLIST.md`) that mention
`chemistry.equilibria` in passing — those are historical audit logs of
a past state, not live references, and this repo's convention is not to
retro-edit shipped docs.

## Import graph analysis

**No interaction with the logged `chemistry`↔`reactions` package cycle**
(`OPEN_WORK.md`'s "Package-level layering" entry, caused by
`chemistry/partition.py` importing `reactions.equilibrium`/
`reactions.stoichiometry`). `equilibria.py` imports nothing from
`reactions/`, and nothing in `reactions/` imports `equilibria.py` —
this is a fully independent dependency thread from that cycle. Worth
being explicit about, since both are "chemistry/ layering" findings
from the same audit and are easy to conflate.

**New edge:** `chemical_equilibrium → chemistry` (via `common_species`,
module-level in the moved file, where today it's an intra-package
`chemistry → chemistry` edge). Checked both directions for a cycle
risk:
- `chemistry/` → zero references to `chemical_equilibrium` anywhere
  (checked directly).
- `reactions/`'s only references to `chemical_equilibrium` are
  docstring cross-refs plus two lazy in-method imports inside
  `reaction_system.py` (`from ..chemical_equilibrium.engines.nr.engine
  import ...` / `...bisection.engine import ...`) — not module-level.

So the new edge stays one-directional; no cycle introduced. It's also
not a new *kind* of dependency: `engine.py` already documents itself
against `chemistry.partition.HenryEquilibrium`/`KspEquilibrium`/
`RaoultEquilibrium` and lazily imports
`reactions.equilibrium.EquilibriumConstraint`. `chemical_equilibrium/`
is the consumer-of-everything engines layer by design, so a forward
dependency on `chemistry/` is the shape it's supposed to have — this
move makes an already-accepted relationship explicit rather than
introducing new coupling.

## Open questions

1. **Re-export from `chemical_equilibrium/__init__.py`?**
   `chemical_equilibrium/__init__.py` currently re-exports
   `BisectionChemicalEquilibriumEngine`, `ChemicalEquilibriumEngine`,
   `EquilibriumResult`, `ionic_strength_from_speciation`,
   `solve_acid_base`. `EquilibriumSet`/`EquilibriumDef` could join that
   list, or stay a deep-only import the way `chemistry.equilibria` is
   today (nothing currently re-exports it, so leaving it deep-only is
   the zero-risk default — this is a polish question, not a blocker).
2. **Keep the `EquilibriumSet` import lazy in `engine.py`, or hoist to
   module level now that it's a same-package sibling?** It's lazy today
   alongside the `reactions.equilibrium` import in the same
   `from_reactions()` body — that adjacent import still needs to stay
   lazy (the `reactions → chemical_equilibrium` lazy-only pattern is
   untouched by this move). Hoisting just the `EquilibriumSet` half
   asymmetrically is possible now that it's risk-free, but splitting
   one lazy-import block into "one eager, one still lazy" may read as
   more confusing than leaving both lazy for consistency. Lean toward
   leaving both as-is unless there's a concrete reason to change import
   timing.
3. **Ship the relocation and the `bsm2_default()` deletion together, or
   split them?** They touch the same file and were found in the same
   pass, so bundling is natural — but they're logically independent
   (you could delete `bsm2_default()` without moving the file, or move
   the file without touching it). Leaning toward bundling: the deletion
   shrinks what has to move, and splitting into two checkpoints inside
   one checklist (rather than two separate phases) captures both the
   convenience and the independence.

## Trigger conditions

No hard trigger — this is a layering/cosmetic cleanup already logged
in `OPEN_WORK.md`, not a bug or blocker. Small and independent enough
to pick up standalone, or to bundle with any future pass that also
resolves the `chemistry`↔`reactions` package cycle (`partition.py`'s
`*Equilibrium` constraint classes) — the two are unrelated in
mechanism (see "Import graph analysis" above) but both are "finish the
`chemistry-reactions-kinetics-cleanup` audit's layering findings"
work, so a single future "layering cleanup" phase could reasonably
carry both.

## How to start one

Per this folder's usual convention
([README.md](../upcoming/README.md#how-to-start-one)): write a checklist file
(`EQUILIBRIUM_SET_RELOCATION_CHECKLIST.md`), cut a branch off `main`
(suggested name: `equilibrium-set-relocation`). Single-phase — no
sequencing, no open API decision blocking a start (Open questions 1-3
above are polish/scoping, not prerequisites). Estimated size: 1 file
moved (~470 lines after the `bsm2_default()` deletion), 2 import
rewrites inside it, 1 import site + 3 docstring refs in `engine.py`, 2
call sites replaced in 1 test file, 1 line in `architecture.md`, 1
docstring clause in `chemistry/__init__.py` — comparable to a couple of
checkpoints elsewhere in this repo's phase docs, not a multi-checkpoint
phase on its own.
