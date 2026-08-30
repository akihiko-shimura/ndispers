# AGENTS.md — working on ndispers itself

For agents *using* the library, read
[docs/llms.txt](docs/llms.txt) instead (full version:
[llms-full.txt](https://ndispers.readthedocs.io/en/latest/llms-full.txt)).
This file is for changing the package — most commonly, adding a crystal.

## Ground rules

- **Never write Sellmeier or d coefficients from memory.** Transcribe them
  from a paper you can actually read (PDF at hand, or Crossref-verifiable),
  and cite it in the docstring. The typical bug in this repo's history is a
  transcription error, not physics. After writing, re-extract the numbers
  from the source independently and diff (see docs/dev/materials_plan.md).
- Computed results changing = minor version bump, even for a bug fix.
- Keep diffs minimal; match the style of the file you are in.

## Adding a crystal (the 1-file recipe)

1. Pick the point-group base class from `ndispers/groups` (Uniax_3m, Uniax_32,
   Uniax_42m, Uniax_4mm, Biax_mm2, Biax_2). No class for your group yet →
   that is a bigger change, read `groups/_base.py` first.
2. Copy the closest existing medium as a template —
   `ndispers/media/crystals/_betaBBO_Eimerl1987.py` (uniaxial) or
   `_KTP.py` (biaxial, one class per principal plane) — into
   `_<Name>_<Source><Year>.py`. Fill in:
   - docstring: composition, point group, transparency range, Sellmeier
     equation as written in the source, validity range, full references;
   - `_wl_range = (lo, hi)` in µm, from the validity statement;
   - the coefficients; `dndT = 0` with a docstring note if the source has no
     temperature term;
   - `_d_ref` (pm/V at the measurement wavelengths) and `_d_note` (how well
     Miller scaling is validated), if the crystal is non-centrosymmetric.
3. Export it in `ndispers/media/crystals/__init__.py` (class name states the
   source: `BetaBBO_Eimerl1987`, not `BBO`). Isotropic media go to
   `media/glasses` instead — the split is by optical isotropy.
4. `uv run python tools/compile_media.py` — regenerates the pre-compiled
   numpy functions (`test_compiled_modules_are_current` fails until you do).
5. Add the new class to the media list in `docs/llms.txt`, then
   `uv run python tools/gen_llms_full.py` (tests enforce both).
6. Tests: add a regression pin comparing `n` at a few wavelengths against
   values printed in the source paper (see `ndispers/tests/test_LBO.py` for
   the pattern). Run `uv run pytest ndispers/tests -q`.

## Releasing

See docs/RELEASING.md. Versions are stamped with
`uv run python tools/set_version.py X.Y.Z` — never edited by hand.

## Repo map

- `ndispers/_baseclass.py` — Medium: all dispersion/phase-matching methods.
- `ndispers/groups/` — point-group classes carrying d_eff closed forms.
- `ndispers/media/{crystals,glasses}/` — one file per medium+source.
- `ndispers/_compiled/` — generated; never edit, regenerate (step 4).
- `docs/dev/` — plans and theory notes; `docs/dev/todo.md` is the queue.
- `apps/` — marimo web apps (GitHub Pages).
