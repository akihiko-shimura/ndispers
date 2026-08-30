"""
Regenerate docs/llms-full.txt, the expanded agent reference: llms.txt, the
full conventions page, and a machine-generated catalog of every medium.

Run after editing docs/llms.txt, docs/conventions.md or any medium:

    uv run python tools/gen_llms_full.py

test_llms_full_txt_is_current fails until you do.
"""
from pathlib import Path

import ndispers as nd

ROOT = Path(__file__).resolve().parent.parent


def catalog_table():
    lines = ['| name | kind | point group | plane | validity (um) | d components |',
             '|---|---|---|---|---|---|']
    for e in nd.catalog():
        rng = f"{e['wl_range'][0]:g}-{e['wl_range'][1]:g}" if e['wl_range'] else 'see docstring'
        lines.append(f"| {e['name']} | {e['kind']} | {e['point_group'] or '-'} | "
                     f"{e['plane'] or '-'} | {rng} | {', '.join(e['d_components']) or '-'} |")
    return '\n'.join(lines)


def build():
    llms = (ROOT / 'docs' / 'llms.txt').read_text()
    conv = (ROOT / 'docs' / 'conventions.md').read_text()
    return (llms
            + '\n\n---\n\n'
            + conv
            + '\n\n---\n\n# Machine-generated media catalog\n\n'
            + 'The output of `ndispers.catalog()`, one row per medium class.\n\n'
            + catalog_table() + '\n')


if __name__ == '__main__':
    out = ROOT / 'docs' / 'llms-full.txt'
    out.write_text(build())
    print(f'wrote {out.relative_to(ROOT)} ({len(out.read_text().splitlines())} lines)')
