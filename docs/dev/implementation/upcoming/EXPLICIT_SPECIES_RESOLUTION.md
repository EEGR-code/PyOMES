# Explicit Species Resolution — Design Note

> Status: design note, not yet started. No branch, no checklist, no code
> yet. Surfaced 2026-09-22 during a conversational investigation of
> whether `chemistry/common_species.py` should move to `PyOMES/databases/`
> (see `CHEMISTRY_REACTIONS_KINETICS_CLEANUP.md`'s package-layering
> thread and `OPEN_WORK.md`'s "Package-level layering" entry for that
> original question). The relocation question turned out to be secondary
> to a more important finding made along the way: species-id resolution
> in three internal call sites is not model-scoped at all — it falls
> back to scanning `common_species.py`'s module namespace, silently, for
> any id it doesn't otherwise recognize. This note is about removing
> that fallback, not about where the file lives; see "Relationship to
> the relocation question" below for how the two connect.

## Motivation

`chemistry/common_species.py` presents itself as a convenience: shared
`Species` objects a model file can import instead of redeclaring
(so cross-file `ReactionSystem` compositions share objects by identity,
per its own docstring). That framing implies it's opt-in — you import
what you need, like any other module.

In practice, three places in the engine internals do not import
specific names from it. They introspect the *entire module namespace*
with `vars(common_species).values())`, building an implicit
id → `Species` catalog that is consulted automatically whenever an id
string shows up that the caller's own explicit declarations don't
cover. A model author who never imports `common_species` at all can
still have its contents silently participate in their model — for
string-based stoichiometry parsing, for gas/liquid partition-model
species resolution, and for a `ControlVolume`'s charge-conservation
accounting. That is the "hidden backend" problem: nothing about a
model's own files tells you which of its species came from where the
user wrote them vs. from this one hardcoded module, and whether an id
resolves at all depends on whether it happens to collide with
something declared there.

## How species resolution works today: three call sites

### 1. `reactions/stoichiometry.py` — string stoichiometry parsing

`_get_common_species()` (`stoichiometry.py:174-180`) builds
`{v.id: v for v in vars(_cs_mod).values() if isinstance(v, Species)}`
from `chemistry.common_species` (imported as `_cs_mod`,
`stoichiometry.py:25`). `_parse_stoichiometry()` (`stoichiometry.py:183`)
seeds its lookup dict with this catalog *first*, then overlays the
caller's own `species=` dict on top (`stoichiometry.py:259-262`:
`lookup = dict(common); lookup.update(species)`). This is reached from
the public API whenever `EquilibriumReaction(stoichiometry="...", ...)`
(`reactions/equilibrium.py:221,231`) or
`KineticReaction(stoichiometry="...", ...)` (`reactions/kinetic.py:90,96`)
is given a string instead of a list of `StoichiometryEntry` objects.

This one already has the right shape for a fix: `species=` is a real,
documented caller override. The only problem is the ambient pre-seed
that runs before it.

### 2. `chemistry/partition.py` — partition-model species fields

`_resolve_species()` (`partition.py:54-71`) takes
`Union[str, Species, None]`. If given a `Species` object it returns it
unchanged; if given a string, it builds the *same* kind of
`vars(common_species).values()` catalog inline (`partition.py:63-65`,
duplicating `_get_common_species()`'s logic — the docstring says as
much) and looks the id up there, with no caller-supplied override path
at all. `HenryEquilibrium.gas_species`/`.liquid_species`
(`partition.py:145-146`) and `RaoultEquilibrium.gas_species`/
`.liquid_species` (`partition.py:310-311`) are typed
`Union[str, Species, None]` and resolved through this function
(`partition.py:245-246, 380-381`). `RaoultEquilibrium.liquid_species`
and `.gas_species` both *default* to the literal string `"H2O"`
(`partition.py:310-311`) — meaning a `RaoultEquilibrium()` constructed
with no arguments at all resolves its species from the ambient module
by design, not as a fallback for an edge case.

This is the hardest of the three: there is no `species=` dict anywhere
in reach to override, so fixing it needs a real API decision, not just
deleting a pre-seed step.

### 3. `core/control_volume.py` — charge-conservation registry

`ControlVolume._collect_species_registry()` (`control_volume.py:840-878`)
builds its species registry primarily from what's actually explicit —
every `Species` object referenced by a `StoichiometryEntry` in the
attached `reaction_system` (`control_volume.py:861-867`). Then, for any
id present in a phase's `n_mol` that *isn't* covered by that (the
docstring's own example: `Cl-`, `Na+`, `K+`, "which never appear in
reaction stoichiometry but must be included for charge conservation
accounting"), it calls `_common_species_catalog()`
(`control_volume.py:880-887`) — another `vars(common_species)` scan —
and includes the species only if the id happens to match an entry
there (`control_volume.py:871-877`). The method's own docstring admits
the alternative: "Species not appearing in any reaction stoichiometry
... are absent from the registry and silently skipped during
conservation accounting" (`control_volume.py:846-849`) — i.e. today's
actual behavior is a three-way split (explicit via reactions / silently
recovered via ambient-catalog name collision / silently dropped) that
depends on a hardcoded module the model author may never have looked
at.

Important correction found while scoping this: `ControlVolume.__init__`
already accepts `chemistry_db: Optional[Any] = None` and stores it as
`self.chemistry_db` (`control_volume.py:191, 219`) — but nothing else in
the file reads that attribute. `ChemistryDatabase.species` (the
explicit, per-model species dict a database author writes —
`databases/database.py:52`) is sitting right there, unused, at exactly
the point where `_common_species_catalog()` reaches for the ambient
module instead.

## Why this matters — concrete failure modes

- **Silent inclusion from the wrong source.** A model that never
  imports `common_species` can still have a bare id like `"CO3--"`
  resolved from it (via any of the three sites above) instead of
  raising "unknown species" — the id happening to match is enough.
  Because `Species` equality is structural, not by identity
  (`chemistry/species.py`'s own docstring: two independently-built
  `Species("H2O", ...)` compare equal), a model author who declares
  their *own* `H2O`/`H_plus`/etc. with a different `MW` override or
  different `atoms` convention gets no signal that a *different*
  object of the same id was available and might be silently preferred
  by a code path they didn't call directly (e.g. `RaoultEquilibrium()`'s
  bare `"H2O"` default).
- **Inconsistent, name-collision-dependent charge accounting.**
  Checked directly: `databases/aqueous.py`'s `_SPECIES` dict
  (`aqueous.py:38-47`) declares exactly the 8 species that appear in its
  own reactions — no `Cl-`/`Na+`/`K+`. So today, whether a spectator ion
  a user puts in `phase.n_mol` is counted for charge conservation is
  governed entirely by whether its id happens to be one of the ~24
  names in `common_species.py` — not by anything the model's own
  `ChemistryDatabase` declared.

## Verified: not currently load-bearing for shipped content

The blast radius of removing the ambient fallback is smaller than it
might look, because the framework's own shipped domain databases don't
rely on it. `databases/aqueous.py`, `bioprocess_basic.py` and
`anaerobic_digestion.py` all build `StoichiometryEntry(species=H2O, ...)`
from `Species` objects imported by name (`aqueous.py:53-79` and
equivalents) — never through the string-stoichiometry parser and never
through a bare id string. The ambient fallback is exercised by: (a)
tests and tutorials using the string-stoichiometry convenience API
without passing `species=`, (b) `RaoultEquilibrium`'s bare-string
default, and (c) any model whose `phase.n_mol` contains a spectator ion
that happens to match a `common_species` name. None of that is the
shipped `databases/*.py` content itself.

## Phase 0: `reactions/stoichiometry.py` (cheap, decided shape)

Remove the `common = _get_common_species()` seed
(`stoichiometry.py:259`) so `lookup` is built from the caller's
`species=` dict alone; delete `_get_common_species()` and the `_cs_mod`
import (`stoichiometry.py:25`) once nothing else uses them. Any
string-stoichiometry construction that currently omits `species=` and
relies on implicit resolution will start raising — needs a repo-wide
audit of `EquilibriumReaction(stoichiometry="...")` /
`KineticReaction(stoichiometry="...")` call sites without a `species=`
argument (expected: tests, tutorials, notebook generators — see
"Verified" above for why the shipped databases aren't at risk here) and
updating each to pass its species explicitly.

## Phase 1: `core/control_volume.py` (moderate — the plumbing already exists)

Change `_common_species_catalog()` to consult `self.chemistry_db.species`
instead of scanning `common_species`. Since `chemistry_db` is already an
optional constructor argument that several call sites already pass
(`templates/stirred_tank/factory.py:120`,
`tests/standalone/test_chemistry_database.py:339,391`), this is mostly
wiring, not new API surface — but it raises one real design question:
what happens when `chemistry_db` is `None` (the common case in
hand-built unit tests that construct a bare `ControlVolume` directly)?
Options: (a) fall through to today's documented "silently skipped"
behavior with no catalog at all — consistent with the docstring's
existing caveat, just without the name-collision recovery path; or (b)
warn when an `n_mol` id can't be resolved and no `chemistry_db` was
given, so the gap is visible instead of silent. Needs its own audit:
which existing tests currently rely on the ambient catalog recovering a
spectator ion and would need `chemistry_db=` added to keep passing.

## Phase 2: `chemistry/partition.py` (hardest — needs an API decision)

No caller-facing override path exists today, so this isn't a
delete-the-seed fix. Two directions, not yet decided:

- **Drop string-id convenience entirely.** Type
  `HenryEquilibrium`/`RaoultEquilibrium`'s `gas_species`/`liquid_species`
  fields as `Species` only (no `Union[str, ...]`, no `"H2O"` default),
  matching the precedent `StoichiometryEntry.species` already sets
  (`stoichiometry.py` — always a real `Species` object, never a bare
  id). Pushes resolution entirely to the call site, consistent with the
  principle behind this whole note. Cost: every current caller passing
  a bare string (need a repo-wide audit; `RaoultEquilibrium()`'s
  zero-argument default is the most visible one) has to start passing
  an explicit `Species` object — likely `H2O` from `common_species`,
  imported explicitly, which is exactly the "opt in by importing" shape
  the module's docstring already claims for itself.
- **Keep string-id convenience, but resolve against a caller-supplied
  dict instead of the ambient module.** Requires threading a species
  dict into `HenryEquilibrium`/`RaoultEquilibrium` construction or into
  whatever attaches them to a `ChemistryDatabase`/`ReactionSystem`
  (unclear where that dict would come from at construction time, since
  these are plain frozen dataclasses built independently of any
  `ChemistryDatabase` today — needs design work, not just plumbing).

First option is simpler and consistent with Phase 0/1's direction; it's
the one to default to unless the audit turns up a strong reason
call-by-string needs to survive.

## Open questions

1. **Phase 0 audit scope.** Exact list of `stoichiometry="..."` call
   sites without `species=` — needed before Phase 0 can ship without
   breaking tests/tutorials silently.
2. **Phase 1's `chemistry_db=None` behavior.** Silent-drop (status quo
   minus the recovery path) vs. a warning. Affects how many existing
   tests need `chemistry_db=` added to keep today's charge-conservation
   coverage.
3. **Phase 2's direction.** Drop string convenience vs. thread an
   explicit dict through — leaning toward "drop," pending the call-site
   audit for `RaoultEquilibrium`'s bare default and any other bare-string
   `HenryEquilibrium`/`RaoultEquilibrium` construction in the repo.
4. **Does any minimal ambient fallback survive on purpose?** Raised and
   rejected in the conversation that produced this note: even water
   autoionization (`H_plus`/`OH_minus`/`H2O`) was considered as a
   "safe" minimal ambient core and rejected, on the grounds that a user
   could declare their own water species under a different convention
   and be unaware the engine was silently reaching past it. This note
   assumes **no** implicit fallback survives in any of the three sites
   — flagged here as a decision, not just an inherited assumption, in
   case a narrower "just water" compromise is reconsidered later.

## Relationship to the `common_species.py` relocation question

The original question (move `common_species.py` to `PyOMES/databases/`?)
is downstream of this note, not parallel to it. Today, `common_species.py`
is partly a convenience library and partly load-bearing ambient engine
state (the three sites above) — that dual role is itself an argument
against moving it anywhere until the ambient-fallback role is gone, since
relocating an implicit dependency doesn't fix the implicitness, and (per
the same investigation) moving it into `databases/` specifically would
have introduced a new `databases` \<-\> `reactions` package cycle, given
`databases/database.py` already imports `reactions.reaction_system` for
real. Once Phases 0-2 land, `common_species.py` reverts to being exactly
what its docstring already claims — an opt-in convenience import, nothing
scans it automatically — and the relocation question can be revisited
on that cleaner premise (at that point it is a `chemistry/`-internal
shared-declarations module, structurally closest to `species.py`, and
the case for leaving it in `chemistry/` gets stronger, not weaker).

## Trigger conditions

No hard trigger — this is a correctness/predictability cleanup, not a
blocking bug. Reasonable pickup points: (a) before or alongside any
future work that adds more species to `common_species.py` (e.g.
precipitation/complexation additions per
`MULTICOMPONENT_COMPLEXATION_AND_PRECIPITATION_PLAN.md`), since every
addition silently grows the ambient catalog's reach; or (b) if a real
model hits the silent-drop or silent-collision behavior in practice.

## How to start one

Per this folder's usual convention
([README.md](README.md#how-to-start-one)): resolve the open questions
above (especially Phase 2's direction), run the Phase 0/1 call-site
audits, write a checklist file
(`EXPLICIT_SPECIES_RESOLUTION_CHECKLIST.md`), cut a branch off `main`
(suggested name: `explicit-species-resolution`). Phases are independently
shippable and increasing in cost/risk (0 < 1 < 2); Phase 0 alone is
low-risk enough to pick up without resolving Phase 1/2's open questions
first.
