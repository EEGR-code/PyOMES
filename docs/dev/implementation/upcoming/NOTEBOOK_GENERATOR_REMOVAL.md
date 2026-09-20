# Notebook generator removal — design discussion

> Status: pre-phase design discussion, 2026-09-15. No branch, no checklist,
> no code yet. Written up from a planning conversation (chat, not a design
> session against code) that started from auditing `demos/aerobic_fermentation/`
> and widened to all 31 notebooks in the repo and their output/generation
> conventions. When this is picked up: follow `README.md`'s "How to start
> one" — write a checklist file from `PHASE_KICKOFF_TEMPLATE.md`, branch,
> implement, ship.
>
> **Update 2026-09-17:** `demos/usecases/` (and its `_generate_notebooks.py`,
> one of the 5 scripts below) was deleted outright, outside this phase.
> That generator's 2 notebooks (`0_README`, `03`) are gone rather than
> converted to standalone files, so this phase now only has 4 generator
> scripts / 19 notebooks left in scope, not 5/21. The duplication this doc
> accepted as a known cost (the `water`/`p1`/`p2`/`p3`/`nh4` network
> retyped in `demos/usecases/_generate_notebooks.py` among other places) is
> correspondingly reduced by one copy. Separately, `demos/` itself (having
> no content left besides its own README) was deleted the same day;
> `demos/aerobic_fermentation_stoichiometry.ipynb` moved to
> `docs/tutorials/reactions/aerobic_fermentation_stoichiometry.ipynb`
> (path below updated to match).
>
> **Update 2026-09-20 (re-audit against `main`):** scope had drifted further.
> Current state of the 32 tracked notebooks (the gitignored `scratch/`
> folder is not counted):
>
> - **2 generator scripts remain**, not 4: `docs/tutorials/ArXiv_preprint/
>   _generate_notebooks.py` (3 notebooks) and `tests/validation/speciation/
>   _generate_notebooks.py` (10 notebooks) — 13 notebooks in scope for
>   conversion. The `ChemicalEquilibriumProtocol` and `SolverProtocols`
>   generators were deleted 2026-09-16 (`tutorials-reorg`); their 7
>   notebooks now live as standalone files in `docs/tutorials/protocols/`,
>   all with outputs embedded.
> - **6 notebooks have no executed cells at all**: `ArXiv_preprint` 01 and
>   03, `D2C_workshop` Example2 and Example3,
>   `reactions/aerobic_fermentation_stoichiometry`, and
>   `templates/batch_fermenter`. Two more (`ArXiv_preprint` 02, `D2C_workshop`
>   Example1) contain a few cells that were never executed.
> - The D2C workshop has 4 notebooks, not 7.
> - The CI rule proposed below ("a code cell has source but empty outputs")
>   would false-positive on cells that ran but print nothing — see Open
>   question 3.
>
> The counts quoted in the "Why", "Alternatives considered" and "Known cost"
> sections below are the September audit's and are left as originally
> written, as a record of that audit.

## Commit discipline for this phase

**The repo owner runs every `git add`, `git commit`, and `git push` for
this phase personally — an assistant executing this plan should not run
those commands.** At each checkpoint: make the file edits, report exactly
what changed (and run the checkpoint's sanity check), then stop and wait
for the owner to review, stage, and commit before moving to the next
checkpoint.

## Goal

Make every notebook in the repo its own single source of truth, committed
with its outputs embedded — matching how packages like biosteam display
plots inline on GitHub without needing local execution. Concretely:

1. Every committed `.ipynb` carries its outputs (execute-and-save before
   commit is normal practice, not an afterthought).
2. The remaining 2 `_generate_notebooks.py` scripts (`ArXiv_preprint` and
   `tests/validation/speciation`) are retired — notebooks they currently
   build become standalone, hand-edited files like any other notebook in
   the repo.
3. A CI check is added that fails a push/PR if any code cell has
   non-empty source but empty `outputs` — catches "edited but forgot to
   execute-and-save" without needing CI to write back to the branch.

## Why

- **Current state is inconsistent.** Of 31 notebooks, 19 have zero
  embedded outputs (nothing renders on GitHub without running them
  locally) and 12 already have outputs committed. `0_README.ipynb` files
  are correctly always empty (pure markdown, no code).
- **The generator scripts make notebooks disposable build artifacts, not
  editable files.** Each `_generate_notebooks.py` builds its notebooks'
  JSON directly and hardcodes `"outputs": []` on every cell, every run.
  Hand-editing a generated notebook in Jupyter is a trap — the next
  `python _generate_notebooks.py` run silently discards it.
- **Drift is already happening, not hypothetical.** `usecases/03_grow_ecoli_on_acetic_acid.ipynb`
  and `docs/tutorials/ArXiv_preprint/01_predict_ph_simple_liquid.ipynb` (moved there, and
  `usecases/04_compare_runtime_by_usecase.ipynb` deleted entirely, by the `tutorials-reorg`
  phase — see `TUTORIALS_REORG_CHECKLIST.md`) contain plotting code but currently have zero
  outputs, while sibling notebooks from the same generator scripts do have outputs. Nothing
  catches this today — no CI touches notebooks at all.
- **No docs-build pipeline exists** (no Sphinx/mkdocs/jupyter-book
  config anywhere in the repo), so there's no separate mechanism
  rendering these notebooks elsewhere — GitHub's own notebook viewer,
  reading the committed `.ipynb` directly, is the only rendering path
  that exists today.

## Alternatives considered and rejected

- **Keep the generators; add execute-and-save discipline plus a CI
  check, nothing else.** Rejected — this permanently carries the
  two-source-of-truth risk (Python script vs. `.ipynb`); a CI check only
  catches the mistake after the fact, it doesn't remove the reason it
  keeps happening.
- **Keep generation, but extract shared boilerplate into an importable
  module instead of inlined strings ("Option 3" in the originating
  conversation).** Investigated concretely, not dismissed on vibes: the
  five-reaction `water`/`p1`/`p2`/`p3`/`nh4` acid-base network is
  retyped verbatim in at least 3 places across
  `docs/tutorials/_generate_notebooks.py` and
  `demos/usecases/_generate_notebooks.py`, and the `_find_repo()`/`_e()`
  setup header is independently duplicated 9 times across the 5
  generator files. That's real, not imagined, duplication. Rejected
  anyway in favor of the simpler option below — recorded here so it
  isn't re-litigated from scratch, only reopened deliberately (see
  Open question 1).
- **Sphinx/mkdocs/jupyter-book.** Out of scope by design — that solves
  "PyOMES should have a documentation website," a separate, bigger
  decision. This phase only fixes how notebooks render when browsed
  directly on GitHub.

## Known cost being accepted

Removing the generators gives up the real code-reuse identified above:

- The `water`/`p1`/`p2`/`p3`/`nh4` network duplication will remain (or
  worsen) once notebooks are hand-authored — correcting a shared
  parameter (e.g. a p$K_a$) will mean manually finding and fixing every
  copy, not editing one source.
- Cross-folder notebook reshuffling (e.g. usecases 01/02/02b/05 moving
  to `docs/tutorials/`, documented in `usecases/_generate_notebooks.py`'s
  own docstring) becomes manual notebook surgery instead of moving a
  Python object between files.

This is an accepted trade-off, not an oversight — recorded so a future
reader doesn't rediscover the duplication and assume it was missed.

## Scope

2 generator scripts left to retire (5 at the September 2026 audit). The
other three are already gone: the `ChemicalEquilibriumProtocol` and
`SolverProtocols` generators were deleted 2026-09-16 by `tutorials-reorg`
(see `TUTORIALS_REORG_CHECKLIST.md`) and `demos/usecases/_generate_notebooks.py`
on 2026-09-17. Counts as of the 2026-09-20 re-audit:

| Generator script | Notebooks produced |
|---|---|
| `tests/validation/speciation/_generate_notebooks.py` | `01`–`08`, `0_README`, `10` (10) |
| `docs/tutorials/ArXiv_preprint/_generate_notebooks.py` | `01`, `02`, `03` (3) |

13 notebooks are currently generator-produced. The other 19 tracked
notebooks are already standalone/hand-authored and unaffected by the
generator removal, but are in scope for the "every notebook carries its
outputs" policy: the 7 in `docs/tutorials/protocols/` (the former
generator-produced ones), the 4 in `docs/tutorials/D2C_workshop/`, 6 in
`docs/tutorials/reactions/` (including
`aerobic_fermentation_stoichiometry.ipynb`), `docs/tutorials/results/
01_exporting_results.ipynb`, and `docs/tutorials/templates/
batch_fermenter.ipynb`.

## Policy adopted

1. Every committed `.ipynb` carries its outputs. No exceptions besides
   genuinely output-free notebooks (pure-markdown `0_README.ipynb`
   files).
2. The remaining 2 `_generate_notebooks.py` scripts are deleted once their
   notebooks are converted to standalone files.
3. CI gains a check that fails the push/PR if any code cell has
   non-empty source but empty `outputs`.
4. No docs-build pipeline is introduced as part of this phase.

## Open questions for whoever picks this up

1. **Duplication handling post-removal.** Leave the repeated
   `water`/`p1`/`p2`/`p3`/`nh4`-style blocks duplicated across notebooks
   (each notebook fully self-contained — easiest to read standalone), or
   introduce a lightweight shared import module later if the
   duplication turns out to bite in practice? Option 2 (this doc) was
   chosen over Option 3 for simplicity; nothing prevents revisiting
   Option-3-style extraction later if drift becomes a real maintenance
   problem.
2. **Conversion order.** Convert folder-by-folder (smallest first —
   `ArXiv_preprint`, 3 notebooks, then `tests/validation/speciation`, 10)
   or both generator folders in one pass? Not decided.
3. **CI mechanism.** A small custom script comparing `outputs` against
   `source` per cell, or an existing tool (`nbval` / `pytest --nbval`,
   which re-executes and diffs — stronger but slower)? Not decided.
   **Caveat found 2026-09-20:** the rule as stated in Goal 3 ("non-empty
   source but empty `outputs`") false-positives on cells that ran but print
   nothing — e.g. `07_iron_oxidation` cell 12 is a bare assignment with
   `execution_count` 7 and no output. Checking for `execution_count: null`
   instead flags only never-run cells, which is the drift this check is
   for (e.g. the `import` cells in `ArXiv_preprint/02`).
4. **`aerobic_fermentation_stoichiometry.ipynb` needs re-running** so
   its committed outputs actually include the Section 6 dynamic-simulation
   plots. ~~Once outputs are embedded, the standalone
   `plt.savefig("aerobic_fermentation/simulation_timeseries.png")` call
   becomes redundant — drop it, or keep it if a linked PNG is still
   wanted for reference elsewhere (not decided).~~ Resolved ahead of this
   phase (2026-09-15): the `savefig` call was dropped from the Section 6
   plotting cell, so this is no longer open. The notebook still needs
   re-running to embed its outputs — as of 2026-09-20 it has no executed
   cells at all, not just a missing Section 6 — that part of this item
   stands.

## Related, separate item from the same conversation

`demos/aerobic_fermentation/` contained 3 stale, unreferenced PNGs
(`mass_balance.png`, `overview.png`, `yield_sensitivity.png`) left over
from an earlier version of the notebook — orphaned, not linked from
anywhere in the repo. Deleted ahead of this phase (2026-09-15), along
with the `plt.savefig` call in `aerobic_fermentation_stoichiometry.ipynb`
that used to target the same folder (see Open question 4).

## Checkpoints

Not yet defined — this is a pre-phase design discussion. Checkpoints
belong in a checklist file (modelled on `PHASE_KICKOFF_TEMPLATE.md`)
once a branch is started, per this folder's convention.
