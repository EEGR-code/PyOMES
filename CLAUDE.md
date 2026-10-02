# PyOMES — working notes for Claude

Standing rules for working in this repo. Phase checklists add phase-specific
rules on top of these; they do not need to repeat them.

## Layout

- `PyOMES/` is the package. `models/` (`vlmodels`) is a separate package in the
  same repo with its own `pyproject.toml`; the tests make it importable.
- `pytest` and CI run `tests/standalone/` and `tests/validation/` (the
  `testpaths` in `pyproject.toml`); `tests/performance/` is run by hand.
  `tests/standalone/test_package_layering.py` enforces which subpackages may
  import which.
- `tools/` holds developer scripts: committed, not packaged, not collected by
  pytest.
- `/scratch/` and `notes/` are gitignored. One-off scripts go in `/scratch/` or the
  session scratchpad, never in the repo.
- `docs/dev/implementation/` holds planning docs: `upcoming/` for open work,
  `shipped/` and `obsolete/` for history, `OPEN_WORK.md` for small logged items.

## Tests

- Full suite: `python -m pytest -p no:cacheprovider -q` (about 2 minutes in the
  owner's terminal, about 15 under Claude's CPU cap).
- It runs before and after any change to `PyOMES/` or `models/`, and the result
  is reported as it is. The owner runs it: Claude hands over
  `python -m pytest -p no:cacheprovider -q | Select-Object -Last 1` and asks for
  the line it prints. Claude may run single test files or `-k` selections that
  finish in a few minutes.
- Claude's tool processes run in a Windows job object capped at about 15% of one
  core; nothing in the repo or in process priority lifts it. Runs take roughly
  6-7 times longer than in the owner's terminal. Do not run notebooks while the
  suite is running.
- If a run is expected to take more than about 5 minutes under the cap (long
  notebooks, sweeps, notebook surveys), do not start it in the background and
  wait. Give the owner the exact command, one line, ready to paste into
  PowerShell 5.1 from the repo root (scripts in the session scratchpad by full
  path), say what output to paste back, and stop anything of Claude's that
  would compete for CPU. A run that could not finish is reported as not run.

## Git

- The repo owner runs every git command that changes state (branch, add, commit,
  push, merge, tag). Read-only local git (`status`, `diff`, `log`, `show`,
  `ls-files`) is fine; nothing that touches the network.
- Hand over git commands one per line. The owner uses PowerShell 5.1, so never
  chain with `&&`.
- Commit messages: exactly two `-m` flags, a short strapline (about 50 characters)
  and a one- or two-sentence body on one line. No attribution or
  `Co-Authored-By` lines.
- The repo lives in OneDrive. A commit may ask "Rename ... index.lock ... failed.
  Should I try again? (y/n)": answer `y`. Run no git command while a commit may be
  waiting.
- Larger work follows the phase convention in
  `docs/dev/implementation/upcoming/README.md` (feature branch, `--no-ff` merge,
  `<name>-shipped` tag). Small changes can use a plain branch without a design
  note or checklist. Older `*-shipped` tags cited in the docs predate GitHub and
  do not exist in git; do not flag them.

## Line endings

- `core.autocrlf=true`; most working-tree files are CRLF. Edits must keep each
  file's existing line endings and trailing bytes. A hook runs
  `tools/check_line_endings.py` after every Edit/Write and reports mixed files.
- Move files with a plain filesystem move, not `git mv`. Before staging a moved
  file, check the old path with `git ls-files --eol <old path>`: if it is
  `i/crlf`, stage the new path with `git -c core.autocrlf=false add <path>`;
  otherwise a plain `git add`. Check `git diff --cached -M --stat` shows a
  rename, not a rewrite.

## Code and docs

- The package has no outside users yet. When something moves or is renamed, old
  import paths stop working: no shims, aliases or re-exports. Update every
  importer, including notebooks and `models/`.
- Inside `PyOMES/`, import from a module's own folder with `from .x import y`
  and from anywhere above it with `from PyOMES.x import y`; never two or more
  dots. `tests/standalone/test_package_layering.py` enforces this.
- Before calling anything unused or safe to delete, search the whole repo across
  all file types (`.py`, `.ipynb`, `.md`, config, runner scripts), check for
  dynamic lookups, and say what was searched.
- Docstrings, comments and error messages describe current behaviour in the
  present tense: no phase names, checkpoint labels, "NEW:" markers or pointers
  to `upcoming/` docs. Check each claim against the code. Tests match some
  message fragments, so search before rewording runtime strings.
- In overview docs, describe classes by their role; name one only when the
  reader needs it to find or use the API.
- Leave `docs/dev/implementation/shipped/` and `docs/dev/ideas/` untouched unless
  asked.
- Bugs or dead code found outside the task's scope go in
  `docs/dev/implementation/OPEN_WORK.md`, not fixed on the spot.
