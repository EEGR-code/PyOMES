# Phase Kickoff Checklist — &lt;PHASE_NAME&gt;

> Copy this file to `docs/upcoming/<PHASE_NAME>_CHECKLIST.md` at the
> start of any phase-shaped work — before writing implementation code. See
> `docs/upcoming/README.md`'s "Branching and tagging convention" for
> the full rationale; this is the fill-in-the-boxes version. Modelled on the
> shape of existing checklists, e.g.
> [`../shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md`](../shipped/CHEMISTRY_UNIFICATION_4_CHECKLIST.md).

## Pre-flight

- [ ] `git status -sb` clean (no stray uncommitted work left over from a
      previous task)
- [ ] `git log origin/main..main --oneline` empty (nothing unpushed sitting
      around from earlier work)
- [ ] Branch created off current `main`: `git checkout -b <phase-name>`
- [ ] This checklist file committed on that branch as the first commit
      (even as a stub) — so partial work is never silent

## During

- [ ] Plan/design doc exists in `docs/upcoming/` (or `docs/design/`)
      describing scope
- [ ] Checkpoints tracked below as they land, one commit per checkpoint
- [ ] **If work stalls or is paused before shipping:** add a status banner
      to the top of the plan doc *immediately* — what's built, what's
      tested, why it stopped, which branch/commit it's on. Don't leave this
      for a future audit to discover.

### Checkpoints

<!-- - [ ] N. <checkpoint description> — file(s) touched, sanity check -->

## Shipping

- [ ] Full test suite green on the branch
- [ ] `git checkout main`
- [ ] `git merge --no-ff <phase-branch> -m "Merge <phase-branch>: <summary>"`
- [ ] `git tag <phase-name>-shipped <commit-hash>`
- [ ] `git push && git push --tags` (both — tags are not pushed by default)
- [ ] `git branch -d <phase-branch>` and
      `git push origin --delete <phase-branch>`
- [ ] Move the plan doc + this checklist to `docs/shipped/`, add a
      "Shipped" banner to both
- [ ] Update `docs/upcoming/README.md`'s "Recently shipped" list
