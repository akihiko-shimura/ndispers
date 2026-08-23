# Releasing

Releases are published to PyPI by the `publish` workflow when a version tag is
pushed. Authentication uses PyPI's trusted publishing (OpenID Connect), so no
API token is stored in the repository, in GitHub secrets or on any developer's
machine.

## One-time PyPI setup

On <https://pypi.org/manage/project/ndispers/settings/publishing/>, add a new
publisher with exactly these values:

| Field | Value |
|---|---|
| Owner | `akihiko-shimura` |
| Repository name | `ndispers` |
| Workflow name | `publish.yml` |
| Environment name | `pypi` |

Then, in this repository's **Settings → Environments**, create an environment
named `pypi`. Adding yourself as a required reviewer there is optional but
recommended: it makes every upload pause for a manual approval, which is the
last chance to stop a release that should not go out.

Once this is done, delete `~/.pypirc` — the token in it is no longer needed and
a long-lived credential on disk is the thing trusted publishing exists to
remove.

## Cutting a release

1. If any medium, `_baseclass.py` or `groups/_base.py` changed, regenerate the
   pre-generated numpy functions and commit them:

   ```
   uv run python tools/compile_media.py
   ```

   (`test_compiled_modules_are_current` fails until this is done.)
2. Update `__version__` in `ndispers/__init__.py`. Follow the project's
   practice: bump the minor version when computed results change, even if the
   change is a bug fix, so the change cannot reach a pinned install silently.
3. Commit and push to `main`.
4. Tag and push:

   ```
   git tag -a v1.2.3 -m "v1.2.3 — short summary"
   git push origin v1.2.3
   ```

5. The workflow then verifies that the tag matches `ndispers.__version__`, runs
   the test suite and the docs build, builds the sdist and wheel, checks them
   with `twine check`, and uploads. If the `pypi` environment requires a
   reviewer, approve the run in the Actions tab.
6. Write the release notes at
   <https://github.com/akihiko-shimura/ndispers/releases>. When results change,
   state the before and after values for the affected media — someone who has
   already computed with an earlier version needs to know whether to recompute.

## Publishing by hand

Should the workflow be unavailable, `twine` reads `~/.pypirc` with
`username = __token__` and an API token as the password:

```
uv build && uvx twine upload dist/*
```

Note that `uv publish` does not recognise the `__token__` username convention
in `.pypirc`; it expects `--token` or `UV_PUBLISH_TOKEN`.
