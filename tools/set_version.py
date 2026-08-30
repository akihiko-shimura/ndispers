"""
Stamp a release version everywhere it lives, so no file is updated by hand:

    uv run python tools/set_version.py 1.2.3

Rewrites __version__ in ndispers/__init__.py, and version / date-released
(today) in CITATION.cff. The publish workflow refuses a tag whose version
disagrees with either file, so a forgotten stamp cannot reach PyPI.
"""
import re
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def stamp(path, pattern, repl):
    p = ROOT / path
    s = p.read_text()
    new, n = re.subn(pattern, repl, s)
    assert n == 1, f'{path}: expected exactly one match for {pattern!r}, got {n}'
    p.write_text(new)
    print(f'{path}: {repl}')


if __name__ == '__main__':
    if len(sys.argv) != 2 or not re.fullmatch(r'\d+\.\d+\.\d+', sys.argv[1]):
        sys.exit(f'usage: {sys.argv[0]} X.Y.Z')
    v = sys.argv[1]
    stamp('ndispers/__init__.py', r'__version__ = "[^"]+"', f'__version__ = "{v}"')
    stamp('CITATION.cff', r'(?m)^version: .+$', f'version: {v}')
    stamp('CITATION.cff', r'(?m)^date-released: .+$', f'date-released: "{date.today()}"')
