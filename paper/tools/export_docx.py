#!/usr/bin/env python3
"""paper.resolved.typ -> paper.docx, through pandoc's real Typst reader.

WHY THIS PATH EXISTS. The HTML route (tools/typst2docx.py, `just docx-html`)
rasterizes every equation into a PNG, because Typst's HTML export drops math
outright and images were the only way back. Pandoc reads Typst natively --
a real evaluator -- and writes NATIVE, editable Word equations. The resolver
has already replaced every project helper with plain Typst, so what this
script feeds pandoc is exactly what a person would read in the source.

Two adaptations pandoc needs, both made here rather than in the resolver,
because they are pandoc's quirks and not properties of the manuscript:

  - `#bibliography(...)` is real Typst (the resolved file compiles standalone,
    references and all), but pandoc's reader parses the call without wiring it
    into citeproc. The call is swapped for a `= <title>` heading and the .bib
    paths and CSL style are handed to pandoc as flags; citeproc then sets the
    reference list under that heading.
  - The CSL style: Typst bundles styles by name ("american-chemical-society");
    pandoc wants a .csl FILE. If <style>.csl or csl/<style>.csl exists in the
    manuscript root it is used; otherwise pandoc's default (Chicago
    author-date) applies, with a printed note -- visibly different citations,
    not silently different numbers.

Every citation key is checked against the .bib file(s) BEFORE conversion.
Citeproc renders a missing key as bold text plus a warning and still exits 0,
which is exactly the kind of shipped-anyway failure this pipeline exists to
refuse.

Usage: uv run python tools/export_docx.py     (via `just docx`)
"""
from __future__ import annotations

import re
import sys
import zipfile
from xml.dom import minidom
from pathlib import Path

# The manuscript root, one level up: this file lives in tools/.
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))

from resolve_typst import FLOAT_PREFIX, _call_span  # noqa: E402
from typst_prose import CITE  # noqa: E402

SRC = ROOT / "paper.resolved.typ"
OUT = ROOT / "paper.docx"


def format_tables(path: Path) -> None:
    """Restore readable scientific tables after pandoc drops Typst sizing.

    Use the 6.5-inch Word text area, content-specific column widths, and
    intact rows. Preserve the XML namespaces and every text/math node.
    """
    with zipfile.ZipFile(path) as archive:
        entries = [(info, archive.read(info.filename)) for info in archive.infolist()]
    doc = minidom.parseString(dict((i.filename, b) for i, b in entries)[
        "word/document.xml"])
    ns = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"

    def children(parent, name):
        return [n for n in parent.childNodes if n.nodeType == n.ELEMENT_NODE
                and n.tagName == "w:" + name]

    def set_child(parent, name, **attrs):
        found = children(parent, name)
        node = found[0] if found else doc.createElementNS(ns, "w:" + name)
        if not found:
            if name.endswith("Pr"):
                parent.insertBefore(node, parent.firstChild)
            else:
                parent.appendChild(node)
        for key, value in attrs.items():
            node.setAttributeNS(ns, "w:" + key, str(value))
        return node

    for table in doc.getElementsByTagName("w:tbl"):
        rows = children(table, "tr")
        header = " ".join(n.firstChild.data for n in
                          rows[0].getElementsByTagName("w:t") if n.firstChild)
        if header.startswith("Component"):
            weights = [20, 35, 45]
        elif header.startswith("Name in text"):
            weights = [13, 17, 35, 35]
        elif "Draws" in header:
            weights = [12, 20, 8, 30, 30]
        elif "Features" in header:
            weights = [10, 12, 17, 10, 10, 14, 13, 14]
        elif "Entrapment FDP" in header:
            weights = [12, 23, 12, 10, 25, 18]
        elif "Definition" in header:
            weights = [6, 25, 55, 14]
        else:
            raise ValueError(f"Set Word column widths for new table: {header}")
        widths = [round(9360 * w / sum(weights)) for w in weights]
        widths[-1] += 9360 - sum(widths)
        props = set_child(table, "tblPr")
        set_child(props, "tblW", w=9360, type="dxa")
        set_child(props, "tblLayout", type="fixed")
        margins = set_child(props, "tblCellMar")
        for side in ("top", "bottom", "left", "right"):
            set_child(margins, side, w=60, type="dxa")
        grid = set_child(table, "tblGrid")
        for col, width in zip(children(grid, "gridCol"), widths, strict=True):
            col.setAttributeNS(ns, "w:w", str(width))
        for index, row in enumerate(rows):
            row_props = set_child(row, "trPr")
            set_child(row_props, "cantSplit")
            if index == 0:
                set_child(row_props, "tblHeader")
            for cell, width in zip(children(row, "tc"), widths, strict=True):
                cell_props = set_child(cell, "tcPr")
                set_child(cell_props, "tcW", w=width, type="dxa")
                set_child(cell_props, "vAlign", val="center")
                for paragraph in cell.getElementsByTagName("w:p"):
                    pp = set_child(paragraph, "pPr")
                    set_child(pp, "spacing", before=0, after=0, line=240,
                              lineRule="auto")
                for run in cell.getElementsByTagName("w:r"):
                    rp = set_child(run, "rPr")
                    set_child(rp, "sz", val=20)
                    set_child(rp, "szCs", val=20)
    # Typst block spacing becomes empty Word paragraphs. Keep those with the
    # following content so a heading cannot be stranded before a figure.
    body = doc.getElementsByTagName("w:body")[0]
    for paragraph in children(body, "p"):
        pp = set_child(paragraph, "pPr")
        styles = children(pp, "pStyle")
        style = styles[0].getAttribute("w:val") if styles else ""
        drawings = paragraph.getElementsByTagName("w:drawing")
        empty = not paragraph.getElementsByTagName("w:t") and not drawings
        if style.startswith("Heading") or drawings or empty:
            set_child(pp, "keepNext", val=1)
        if style == "ImageCaption":
            set_child(pp, "keepLines", val=1)
    xml = doc.toxml(encoding="utf-8")
    with zipfile.ZipFile(path, "w") as archive:
        for info, data in entries:
            archive.writestr(info, xml if info.filename == "word/document.xml" else data)


def split_bibliography(src: str) -> tuple[str, list[str], str | None]:
    """Swap the #bibliography call for its heading; return (src, paths, style).

    The heading is followed by a `#block[]<refs>` anchor: citeproc sets the
    reference list inside a Div with id "refs", and WITHOUT one it appends
    the list at the very end of the document -- which, now that the resolver
    keeps the call in the PDF's position, would strand the references after
    the entire SI. Pandoc's Typst reader turns the anchor into a Span, and
    tools/refs_div.lua promotes it to the Div citeproc looks for.

    Pure and separate from the pandoc run so the tests can hold it still.
    A resolved file without a bibliography passes through unchanged.
    """
    call = _call_span(src, "#bibliography(")
    if call is None:
        return src, [], None
    paths = re.findall(r'"([^"]+\.(?:bib|yml|yaml|json))"', call)
    style = re.search(r'style:\s*"([^"]+)"', call)
    title = re.search(r"title:\s*\[([^\]]*)\]", call)
    heading = (f"= {title.group(1) if title else 'Bibliography'}"
               "\n\n#block[]<refs>")
    return (src.replace(call, heading), paths,
            style.group(1) if style else None)


def check_citations(src: str, bib_paths: list[Path]) -> list[str]:
    """Every @key cited but not defined by the bibliography files.

    Float and section prefixes are cross-references, already resolved to
    literal text by the resolver; anything else surviving as `@key` is a
    citation and must have an entry, or citeproc ships it as bold prose.
    typst_prose.CITE decides what a key looks like -- one authority, so a
    trailing period is prose here exactly as it is everywhere else. One
    guard on top of it: an `@` that is escaped (`\\@scripps` in the author
    email) or mid-word (`"mailto:pgarrett@scripps.edu"`) is not citation
    syntax to Typst and must not be one here -- both live in the back
    matter, and each briefly failed this check as "@scripps not in the
    bibliography".
    """
    cited = {m.group(0)[1:] for m in re.finditer(r"(?<![\w\\:./-])" + CITE, src)
             if m.group(0)[1:].split(":", 1)[0] not in FLOAT_PREFIX}
    known: set[str] = set()
    for p in bib_paths:
        known |= {m.group(1) for m in
                  re.finditer(r"(?m)^@\w+\{([^,\s]+)\s*,",
                              p.read_text(errors="replace"))}
    return sorted(cited - known)


def main() -> int:
    if not SRC.is_file():
        print("error: paper.resolved.typ is missing -- run: just resolve",
              file=sys.stderr)
        return 1
    src, bib, style = split_bibliography(SRC.read_text())

    args = ["--resource-path", str(ROOT)]
    if bib:
        missing = check_citations(src, [ROOT / b for b in bib])
        if missing:
            print("error: cited but not in the bibliography: "
                  + ", ".join(f"@{k}" for k in missing)
                  + " -- citeproc would ship each as bold prose and exit 0.",
                  file=sys.stderr)
            return 1
        # The Lua filter must precede --citeproc: pandoc applies filters in
        # command-line order, and the refs anchor has to be a Div before
        # citeproc goes looking for one.
        args += ["--lua-filter", str(ROOT / "tools" / "refs_div.lua")]
        args += ["--citeproc"]
        args += [f"--bibliography={ROOT / b}" for b in bib]
        csl = next((p for p in (ROOT / f"{style}.csl",
                                ROOT / "csl" / f"{style}.csl")
                    if style and p.is_file()), None)
        if csl:
            args += ["--csl", str(csl)]
        elif style:
            print(f'note: no {style}.csl in the manuscript root; citations '
                  f"use pandoc's default style (Chicago author-date). Drop "
                  f"the CSL file there to match the PDF.")

    import pypandoc
    pypandoc.convert_text(src, "docx", format="typst",
                          outputfile=str(OUT), extra_args=args)
    format_tables(OUT)

    # Trust, then verify: count what actually landed in the file.
    with zipfile.ZipFile(OUT) as z:
        doc = z.read("word/document.xml").decode("utf-8", errors="replace")
    n_math = doc.count("<m:oMath>")
    n_tbl = doc.count("<w:tbl>")
    n_img = doc.count("<pic:pic ") + doc.count("<pic:pic>")
    mb = OUT.stat().st_size / 1e6
    print(f"wrote {OUT.name} ({mb:.1f} MB) -- {n_math} native equations, "
          f"{n_tbl} tables, {n_img} images"
          + (f", references set from {', '.join(bib)}" if bib else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
