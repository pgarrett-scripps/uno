#!/usr/bin/env python3
"""Fail if any text in the compiled PDF sits outside the printable area.

This exists because of a real defect. Tables S2 and S3 are 24 rows each and sat
inside non-breakable `#figure` blocks, so instead of paging they overflowed the
text block: the last data row printed on top of the caption and the page number,
and was unreadable. Nothing in `typst compile` warns about it, the source looks
correct, and the only way to see it is to look at where the glyphs landed.

Usage:
    check_overflow.py paper.pdf [--margin PT]

Exits non-zero and lists the offending words if any are found. Requires
`pdftotext` (poppler-utils); skips silently with exit 0 if it is unavailable, so
the check never blocks a machine that lacks it.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# How close to the physical page edge counts as "outside". The text block ends
# well before this; anything past it is a layout escape, not tight typesetting.
DEFAULT_MARGIN_PT = 25.0

PAGE_RE = re.compile(r'<page width="([\d.]+)" height="([\d.]+)">(.*?)</page>', re.S)
WORD_RE = re.compile(
    r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">'
    r"([^<]*)</word>"
)


def overflowing(pdf: Path, margin: float) -> list[tuple[int, float, float, str]]:
    """Return (page, ymax, page_height, text) for every word past the margin."""
    if shutil.which("pdftotext") is None:
        return []
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "bbox.html"
        subprocess.run(
            ["pdftotext", "-bbox", str(pdf), str(out)],
            check=True, capture_output=True,
        )
        html = out.read_text(encoding="utf-8", errors="replace")

    bad: list[tuple[int, float, float, str]] = []
    for page_no, m in enumerate(PAGE_RE.finditer(html), start=1):
        height = float(m.group(2))
        for w in WORD_RE.finditer(m.group(3)):
            ymax = float(w.group(4))
            if ymax > height - margin:
                bad.append((page_no, ymax, height, w.group(5)))
    return bad


def main() -> int:
    args = [a for a in sys.argv[1:]]
    margin = DEFAULT_MARGIN_PT
    if "--margin" in args:
        i = args.index("--margin")
        margin = float(args[i + 1])
        del args[i:i + 2]
    if not args:
        print(__doc__.strip().splitlines()[0], file=sys.stderr)
        print("usage: check_overflow.py paper.pdf [--margin PT]", file=sys.stderr)
        return 2

    pdf = Path(args[0])
    if not pdf.exists():
        print(f"check_overflow: {pdf} not found", file=sys.stderr)
        return 1
    if shutil.which("pdftotext") is None:
        print("check_overflow: pdftotext not installed, skipping")
        return 0

    bad = overflowing(pdf, margin)
    if not bad:
        print(f"check_overflow: no text outside the printable area ({pdf.name})")
        return 0

    print(f"OVERFLOW: {len(bad)} words sit outside the printable area of {pdf.name}")
    seen_pages = []
    for page_no, ymax, height, text in bad[:20]:
        print(f"  p{page_no}  y={ymax:.1f} (page height {height:.1f})  {text!r}")
        if page_no not in seen_pages:
            seen_pages.append(page_no)
    if len(bad) > 20:
        print(f"  ... and {len(bad) - 20} more")
    print(f"  affected pages: {', '.join(str(p) for p in seen_pages)}")
    print("  usual cause: a long table in a non-breakable #figure block.")
    print("  fix: #show figure: set block(breakable: true)")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
