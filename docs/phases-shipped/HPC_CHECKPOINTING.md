# HPC Checkpointing — Planning Note

> **Status:** Shipped 2026-06-08. Tag `hpc-checkpointing-shipped`.
> C1 (snapshot audit), C2 (start_t_h + BatchResult.concat),
> C3 (save/load mode='data'), C4 (mode='data+env') all landed.
> Original trigger-gated note below for reference.

> **Original status:** Trigger-gated. Surfaced 2026-05-28 during the
> post-simulation-class exploration session, motivated by
> resuming long simulations across HPC time-limit boundaries.
> The work splits into a **verification half** (audit what
> `Simulation.snapshot()` / `ControlVolume.snapshot()`
> actually preserve) and an **API half** (design the
> `save_checkpoint` / `load_checkpoint` surface). Either can
> be picked up independently; the verification half is the
> cheaper starter.

## Context

[`Simulation.snapshot()`](../../src/core/simulation.py#L277-L295)
and [`ControlVolume.snapshot()`](../../src/core/control_volume.py#L787-L796)
exist and deep-copy CV state (phases / `n_mol`). They share
links, controllers, profiles, the solver, and the recorder —
those are treated as behaviour/config, not state.

For HPC time-limit recovery the question is different: not
"give me a copy of the CV at this instant" but "if my job
hits the wallclock limit, can I resume from disk on the next
job allocation and continue where I left off — bit-for-bit
when the env is preserved, loudly-erroring when it drifts."

The current `snapshot()` machinery is a building block but
not the answer:

- Snapshot is in-memory only. HPC resume needs persistence to
  disk.
- Snapshot intentionally shares the recorder. HPC resume
  needs the mid-array recorder state too (or a streaming
  recorder that's already on disk — see
  [RUN_HISTORY.md](RUN_HISTORY.md)).
- Snapshot has no environment fingerprint. A resume against a
  changed Python / numpy / scipy is silent today (the pickle
  either loads-but-wrong, or fails with a cryptic
  `ModuleNotFoundError`).

## The verification half — what does `snapshot()` capture?

Audit task; no design decisions. Read each `snapshot()`
implementation and produce a coverage table for every piece
of mid-run state. The known suspects are:

| Mid-run state | Lives on | Captured by `snapshot()` today? |
|---|---|---|
| `_t_h` accumulator | `Simulation` | unknown — verify |
| `RunContext.is_running` flag | `Simulation` | unknown — verify (probably needs explicit handling) |
| Controller `_controller_sample_state` | `Simulation` | unknown — verify |
| Controller `last_action` cache | `Simulation` | unknown — verify |
| `BatchRecorder` mid-array state | `Simulation` | shared (not copied); needs explicit handling for resume |
| `AccuracyMonitor` per-CV state | each `ControlVolume` | unknown — verify |
| `SpeciationEngine._last_logH` warm-start | each `ReactionSystem`'s engine | unknown — verify |
| `Phase.n_mol` arrays | each `Phase` | yes (deep-copied) |
| Step-counter / iteration index | none today | needs new field if HPC resume requires it |

The verification half produces a small test suite asserting
exactly which fields snapshot copies vs shares, plus a
short doc summarising the coverage. Anything in the
"unknown" rows that turns out to be shared needs to be
explicitly deep-copied (or explicitly chosen as
config-not-state with a justification).

Estimated: ~50 LOC of new tests + a doc table + likely 1-3
small fixes if gaps are found.

## The API half — what does `save_checkpoint` look like?

Drawing on the discussion summary from the 2026-05-28
session, the levels of fidelity available to Python code are:

| `mode=` | What's bundled | Detection of drift on resume | Bit-for-bit when env matches? |
|---|---|---|---|
| `"data"` | pickle of `Simulation` (object graph) | none | **silent corruption possible** if classes renamed or numpy dtypes shift |
| `"data+env"` | above + `pip freeze` + Python version + platform string | refuse-to-load on version mismatch | yes on same env; **loud refusal** on drift |
| `"source"` | above + tarball of every VLsim module the run imported + the user script | content-hash mismatch on any code edit | yes; code-side drift caught at load |
| `"full"` | above + SHA-256 of every imported `.so`/`.pyd` (numpy, scipy, BLAS) | detects compiled-extension drift | yes; C-extension drift caught |

The CRIU-style memory-image option is **explicitly
out-of-scope** — it requires OS / kernel cooperation, is
Linux-only, is not portable across hardware, and is not
something a Python library can deliver. Users wanting that
level of fidelity use CRIU directly or a VM/container
snapshot, both of which are external to VLsim.

### Numerical reproducibility floor

No checkpoint level changes the floor; the floor is set by
the runtime. Two persistent gotchas to surface in the API:

- **Threaded BLAS summation order is non-deterministic across
  thread counts.** Same code, same Python, same numpy ABI,
  different `OMP_NUM_THREADS` ⇒ different ULPs. The
  checkpoint metadata should stamp the active thread count
  and the API should refuse (or loudly warn) on resume under
  a different value.
- **Checkpoint must be taken at step boundaries.** Mid-step
  state (inside a scipy LP solve, inside a speciation Newton
  iteration) is not resumable; the snapshot's atomicity is
  one `Simulation` step. The API should expose only
  `sim.save_checkpoint()` between `sim.run()` calls, and
  document the constraint.

### Proposed API

```python
class Simulation:
    def save_checkpoint(
        self,
        path: str | Path,
        *,
        mode: Literal["data", "data+env", "source", "full"] = "data+env",
        include_recorder: bool = True,
        compression: Literal["none", "gzip", "zstd"] = "gzip",
    ) -> None: ...

    @classmethod
    def load_checkpoint(
        cls,
        path: str | Path,
        *,
        verify: Literal["strict", "warn", "skip"] = "strict",
    ) -> "Simulation": ...

    def run(
        self,
        tau_h: float,
        n_steps: int,
        *,
        start_t_h: float = 0.0,   # ← added for checkpoint continuation
    ) -> Any: ...
```

`start_t_h` shifts the time grid: the run covers
`[start_t_h, start_t_h + tau_h]` instead of `[0, tau_h]`.
Default `0.0` preserves all existing behaviour.
`load_checkpoint` exposes the saved wall-clock time as
`sim.checkpoint_t_h` so the caller can pass it directly:

```python
sim = Simulation.load_checkpoint("my_checkpoint/")
result2 = sim.run(tau_h=5.0, n_steps=500, start_t_h=sim.checkpoint_t_h)
full = BatchResult.concat(result1, result2)
```

> **Decision pinned 2026-06-04:** `run()` receives a
> `start_t_h` parameter rather than accumulating time
> internally (which would reverse decision 10). This is
> backward-compatible (`start_t_h=0.0` by default), avoids
> reworking `BatchRecorder` (single-use by design), and keeps
> controller-reset semantics unchanged. One residual: after
> resume with `start_t_h=T`, a controller with
> `sample_period_h=P` fires at `T + ε, T + P + ε, …`
> (offset by one dt at the boundary). Acceptable for HPC use.

Layout on disk (a directory rather than a single file, so
the source bundle and env metadata are inspectable):

```
my_checkpoint/
  simulation.pkl                  ← always present
  manifest.json                   ← always present (includes t_h at save)
  environment/                    ← mode >= "data+env"
    pip_freeze.txt
    python_version.txt
    platform.txt
    omp_num_threads.txt
  source/                         ← mode >= "source"
    vlsim/                        ← tarball of every imported VLsim module
    user_script.py                ← copy of __main__ if available
  binaries.sha256                 ← mode == "full"
```

`compression` controls only `simulation.pkl` (the dominant
size); the source bundle is uncompressed for inspectability.

### Open design questions

> **Decisions pinned 2026-06-04:**
>
> - **Object graph scope:** `save_checkpoint` pickles the full
>   `Simulation` object graph — controllers, links, profiles,
>   and solver are included. `load_checkpoint` returns a
>   self-contained simulation ready to resume with `run()`.
>   This differs deliberately from `Simulation.snapshot()`,
>   which shares config; checkpointing is a different use case
>   (disk-based HPC resume, not in-memory state isolation).
>   Document that objects with unpicklable external state
>   (file descriptors, database handles) need
>   `__getstate__`/`__setstate__`; typical pure-Python
>   VLsim controllers never hit this.
>
> - **Mid-run checkpoint trigger:** `save_checkpoint` is called
>   explicitly by the user between `sim.run()` calls (user
>   drives the loop). No `checkpoint_every_n_steps` parameter
>   is added to `run()` in this iteration — that is a natural
>   follow-on once the basic API exists. The verification half
>   must confirm `t_h` accumulates correctly across sequential
>   `run()` calls.

1. **`include_recorder=`.** **Decision:** always pickle the
   recorder (`include_recorder=True` default). The
   `StreamingFileRecorder` optimisation (skip its mid-array
   buffer) is deferred until that recorder ships; a follow-up
   note is in [RUN_HISTORY.md](RUN_HISTORY.md).
2. **Resume-side validation policy.** **Decision (pinned
   2026-06-04):** `verify="strict"` (default) raises on
   Python version and VLsim version mismatch; warns on
   third-party package versions, platform string, and
   OMP_NUM_THREADS mismatch. `verify="warn"` downgrades
   everything to warnings. `verify="skip"` bypasses all
   checks. Rationale: Python and VLsim version mismatches
   risk silent pickle corruption or missing class attributes;
   third-party and platform drift is usually safe for
   pure-Python code.
3. **Deterministic seed handling.** **Resolved (2026-06-04):**
   No `numpy.random` usage anywhere in `src/` — the codebase
   is fully deterministic. No action required.
4. **Cross-machine reproducibility.** **Decision (pinned
   2026-06-04):** Record `numpy.__config__` BLAS information
   in the manifest. In `"data+env"` mode, warn (do not raise)
   on BLAS mismatch between save and load. Document clearly
   that bit-for-bit identity requires identical compiled
   extensions; `"full"` mode (second iteration) is the path
   for binary-hash-level fidelity.

## Trigger conditions

Any one of:

1. **A real HPC time-limit failure.** A user runs a long
   simulation, hits the wallclock, and discovers there's no
   resume path. Ship the verification half + `mode="data"` /
   `mode="data+env"` as the first iteration.
2. **A multi-day batch study.** Anyone running >24h
   simulations as a parameter sweep or Monte-Carlo benefits
   from periodic checkpoints regardless of HPC. Same first
   iteration covers this.
3. **An audit / reproducibility requirement.** A
   publication, regulatory submission, or
   industrial-validation context needs `mode="source"` or
   `mode="full"`. That's the second iteration's territory.
4. **Bundled with [RUN_HISTORY.md](RUN_HISTORY.md).** The
   `StreamingFileRecorder` is the time-series side of HPC
   resume; this note is the orchestrator side. Shipping
   them together gives a complete checkpoint story.

## Out of scope

- **Memory-image / CRIU-style checkpointing.** OS-level,
  not Python-library territory.
- **Cross-Python-version resume.** Pickle protocol is
  versioned and the API will refuse cross-version loads;
  upgrading Python requires re-running from `t=0` or
  manually replaying via a higher-level mechanism (which
  doesn't exist today).
- **Encrypted / signed checkpoints.** Not a current use
  case.
- **Resume mid-step.** Atomicity is per-step.

## Relationship to other phases

- **Bundles naturally with**
  [RUN_HISTORY.md](RUN_HISTORY.md) — the
  `StreamingFileRecorder` is the time-series persistence
  half; this note is the orchestrator-state persistence
  half. Both together produce a complete crash-resilient
  recorder.
- **Verification half is independent** — audit the current
  `snapshot()` coverage today; doesn't block on anything.
- **Independent of** every other shipped or trigger-gated
  phase. The API half does not need any source-level
  refactor to land — `Simulation` is already the orchestrator
  surface the API hangs off.
- **Follows** `simulation-class` (shipped 2026-05-27), which
  installed `Simulation.snapshot()` in its current shape.
