#!/usr/bin/env python3
"""Check the manuscript prose against the mechanical rules in STYLE.md.

Only the rules a machine can judge. "Simpler is better" is not checkable and is
not checked; a British spelling is, and so is an em dash. Everything here runs on
the same cleaned prose the word count and readability report use, so citations,
figures, tables, math, and code are already out of the way and cannot trigger a
false hit.

Two severities. ERRORS are rules with no legitimate exception in this manuscript,
and they exit non-zero so `just prose-check` can gate a commit. WARNINGS are
judgement calls, reported with counts and locations so they can be skimmed and
ignored. A style checker that fails the build over the word "very" gets disabled
within a week, so it does not.

Usage:
    python3 prose_check.py               # main text + SI
    python3 prose_check.py --strict      # warnings become errors too
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import readability

HERE = Path(__file__).resolve().parent

# --- ERRORS: no legitimate exception --------------------------------------

# British -> American.
#
# An explicit list, NOT a general "-ise -> -ize" rule. The general rule looks
# tempting and is wrong: it fires on every domain word that happens to end in
# those letters, which in the manuscript this came from meant flagging the name of
# the software itself (dnoise, denoise) on every mention. A checker that cries
# wolf on the project's own vocabulary gets switched off. Add words as you hit
# them rather than trying to be clever.
_ISE = [
    "normalise", "analyse", "optimise", "characterise", "summarise", "recognise",
    "organise", "standardise", "minimise", "maximise", "utilise", "emphasise",
    "visualise", "realise", "categorise", "prioritise", "generalise",
    "specialise", "stabilise", "polymerise", "ionise", "oxidise", "neutralise",
    "localise", "randomise", "digitise", "hydrolyse", "catalyse", "paralyse",
    "criticise", "memorise", "familiarise", "harmonise", "synchronise",
    "hypothesise", "apologise", "colourise",
]


def _americanize(w: str) -> str:
    return w.replace("yse", "yze").replace("ise", "ize").replace("isation", "ization")


# Forms that collide with a correct American word and must never be flagged.
# "analyses" is the ordinary plural of "analysis"; only the verb sense is British,
# and nothing here can tell them apart.
_COLLIDES = {"analyses", "practises", "practising", "practised", "programmes"}


def _british_words() -> dict[str, str]:
    out = {}
    for w in _ISE:
        stem = w[:-1]                       # normalise -> normalis
        for form in (w, stem + "ed", stem + "ing", w + "s", stem + "ation"):
            if form in _COLLIDES:
                continue
            out[form] = _americanize(form)
    out.update({
        "colour": "color", "colours": "colors", "coloured": "colored",
        "behaviour": "behavior", "behaviours": "behaviors",
        "favour": "favor", "favours": "favors", "favoured": "favored",
        "labour": "labor", "centre": "center", "centres": "centers",
        "metre": "meter", "metres": "meters", "fibre": "fiber",
        "catalogue": "catalog", "grey": "gray", "artefact": "artifact",
        "artefacts": "artifacts", "modelling": "modeling", "modelled": "modeled",
        "labelling": "labeling", "labelled": "labeled",
        "signalling": "signaling", "signalled": "signaled",
        "towards": "toward", "whilst": "while", "amongst": "among",
        "learnt": "learned", "practise": "practice", "defence": "defense",
        "programme": "program", "sulphur": "sulfur", "haemoglobin": "hemoglobin",
        "oedema": "edema", "foetal": "fetal",
    })
    return out


BRITISH = _british_words()

# --- WARNINGS: judgement calls --------------------------------------------

VERBOSE = {
    "utilize": "use", "utilizes": "uses", "utilized": "used",
    "in order to": "to", "due to the fact that": "because",
    "it should be noted that": "(cut it)", "it is worth noting that": "(cut it)",
    "the fact that": "that", "a number of": "several",
    "in the event that": "if", "at this point in time": "now",
    "each and every": "every", "first and foremost": "first",
    "very": "(cut it)", "quite": "(cut it)", "clearly": "(cut it)",
    "obviously": "(cut it)", "importantly": "(cut it)",
    "basically": "(cut it)", "essentially": "(cut it)",
}
DOUBLE_HEDGE = re.compile(
    r"\b(may|might|could|can)\s+(possibly|potentially|perhaps|conceivably)\b", re.I
)
MAX_SENTENCE_WORDS = 40      # a hard "this is a run-on" line, not the 25-word aim
OPENER_RUN = 3               # N consecutive sentences opening with the same word
# Words common enough that repeating them is invisible; only flag beyond these.
COMMON = set("""
the a an and or but of in on at to for with from by as is are was were be been
this that these those it its we our they their he she them his her him us you your
not no if then than so such which who whom whose what when where while during
can could may might will would shall should must do does did have has had
one two three all both each any some more most other another same different
into over under between within across after before above below through
""".split())


def sentences(text: str) -> list[str]:
    """Split into sentences, protecting the abbreviations readability knows about."""
    t = text
    for ab in readability._ABBR:
        t = t.replace(ab, ab.replace(".", "\x00"))
    t = re.sub(r"(?<=\d)\.(?=\d)", "\x00", t)
    parts = re.split(r"(?<=[.!?])\s+", t)
    return [p.replace("\x00", ".").strip() for p in parts if len(p.split()) >= 2]


def _ctx(text: str, i: int, w: int = 45) -> str:
    return "..." + re.sub(r"\s+", " ", text[max(0, i - w):i + w]) + "..."


def check(label: str, text: str, spellable: str | None = None
          ) -> tuple[list[str], list[str]]:
    """`text` is the cleaned prose. `spellable` is the same prose with inline code
    removed rather than unwrapped, and is what the spelling check runs on.

    They have to differ. readability.clean() turns `Ms1.Normalised` into a bare
    word, because journals count an inline-code term as a word. Spell-checking
    that word then flags a DIA-NN column name as a British spelling, which is not
    something the author can act on.
    """
    spellable = text if spellable is None else spellable
    errors, warnings = [], []
    sents = sentences(text)

    # --- errors ---
    for m in re.finditer(r"—", text):
        errors.append(f"{label}: em dash  {_ctx(text, m.start())}")

    for m in re.finditer(r"\b[A-Za-z]+\b", spellable):
        fix = BRITISH.get(m.group(0).lower())
        if fix:
            errors.append(f"{label}: British spelling {m.group(0)!r} -> {fix}  "
                          f"{_ctx(spellable, m.start())}")

    # NB: the doubled-word check is NOT here. It runs on the raw source in
    # check_house_style(), because clean() replaces each dropped construct with a
    # space, and `and $N_"human"$ and` then collapses to a literal `and and`.
    # A real doubled word is adjacent before the stripping, not after it.

    # --- warnings ---
    for s in sents:
        n = len(s.split())
        if n > MAX_SENTENCE_WORDS:
            warnings.append(f"{label}: {n}-word sentence  \"{s[:90]}...\"")

    for phrase, fix in VERBOSE.items():
        for m in re.finditer(rf"\b{re.escape(phrase)}\b", text, re.I):
            warnings.append(f"{label}: {phrase!r} -> {fix}  {_ctx(text, m.start())}")

    for m in DOUBLE_HEDGE.finditer(text):
        warnings.append(f"{label}: double hedge {m.group(0)!r}  {_ctx(text, m.start())}")

    # repeated sentence openers
    run, first = 1, 0
    for i in range(1, len(sents) + 1):
        same = (i < len(sents)
                and sents[i].split()[:1] == sents[i - 1].split()[:1]
                and sents[i].split()[:1])
        if same:
            run += 1
        else:
            if run >= OPENER_RUN:
                word = sents[first].split()[0]
                warnings.append(f"{label}: {run} sentences in a row open with {word!r}")
            run, first = 1, i

    # a distinctive word repeated inside one sentence
    for s in sents:
        seen: dict[str, int] = {}
        for w in re.findall(r"[A-Za-z][A-Za-z'-]{3,}", s.lower()):
            if w in COMMON:
                continue
            seen[w] = seen.get(w, 0) + 1
        for w, c in seen.items():
            if c >= 3:
                warnings.append(f"{label}: {w!r} appears {c}x in one sentence  "
                                f"\"{s[:80]}...\"")

    semis = text.count(";")
    if semis:
        warnings.append(f"{label}: {semis} semicolon(s); STYLE.md asks for sparing use")

    return errors, warnings


def no_code(src: str) -> str:
    """Drop inline-code spans outright, so identifiers are not spell-checked."""
    return re.sub(r"`[^`]*`", " ", src)


# Acronyms that are never expanded because expanding them would be absurd.
ACRONYM_OK = {
    # file formats and computing
    "PDF", "CSV", "TSV", "JSON", "HTML", "XML", "URL", "API", "CPU", "GPU",
    "RAM", "SSD", "OS", "ID", "IDS", "AI", "MIT", "BSD", "GNU",
    # units and quantities
    "GB", "MB", "KB", "TB", "MS", "NS", "PPM", "RPM", "SD", "SE", "CI", "CV",
    # places, orgs, identifiers
    "USA", "UK", "EU", "CA", "NY", "ORCID", "DOI", "ISO", "UTC", "PHD", "NIH",
    # near-universal in science
    "DNA", "RNA", "PCR", "FDR", "PCA",
    # --- this manuscript's field: mass spectrometry / proteomics ---------
    # Standard in every JPR paper, and a reader of this journal who needed
    # "MS1" expanded would not be reading it. Expanded once at first use where
    # the manuscript does so; this list is what must NOT be nagged about.
    "MS", "MS1", "MS2", "LC", "RT", "IM", "CCS", "FAIMS", "CV", "PSM", "PFM",
    "FDP", "TDC", "AUC", "ROC", "MBR", "LFQ", "TMT", "DIA", "DDA", "UPS", "UPS1",
    "BH", "KS", "PPM", "FASTA", "TSV", "MZML", "PRIDE", "PXD",
    # methods and models named by acronym throughout
    "LDA", "MLP", "CNN", "SVR", "GBDT", "RANSAC", "RNG", "ML", "CLI", "RSS",
    # repo files and section pointers that survive the inline-code stripper
    "RESULTS", "TODO", "SI",
}


KIND = {"fig": "Figure", "tbl": "Table", "eq": "Equation"}


def check_reference_order(sources: dict[str, str]) -> list[str]:
    """Figures and tables should be first cited in numerical order.

    Typst numbers them by order of appearance in the source, so the number a
    reader sees is fixed by where the #figure sits. Journals then require the text
    to reach them in that order, and citing Figure 3 before Figure 2 is a
    copy-editing return at many publishers. Reordering prose during revision is
    exactly how it happens.

    A WARNING rather than an error, because there is one defensible exception: a
    conventions or overview paragraph that legitimately forward-references a later
    figure. Take it seriously anyway, since a journal will.

    Two details make this correct rather than approximately right. A reference
    inside a figure's own caption does not count as the text reaching it, so
    figure blocks are stripped before looking for citations. And each document is
    scored on its own sequence, because the SI restarts at S1.

    Reported one line per offending EARLY citation, not one per figure it jumps
    ahead of. A single early mention of Figure 7 puts six later figures out of
    order, and six near-identical messages describing one edit is noise.
    """
    out = []
    for name, src in sources.items():
        # Numbering order: where each #figure/#table actually sits.
        defined = [m.group(1) for m in
                   re.finditer(r"<((?:fig|tbl):[A-Za-z0-9_-]+)>", src)]
        number, label_of = {}, {}
        for kind in ("fig", "tbl"):
            seq = [d for d in defined if d.startswith(kind + ":")]
            for i, label in enumerate(seq):
                number[label] = i + 1
                label_of[(kind, i + 1)] = label

        # Citation order, ignoring cross-references made from inside a caption.
        prose = readability._strip_balanced(src, "#figure(")
        cited: list[str] = []
        for m in re.finditer(
                r"(?:@|#refn\(\s*<)((?:fig|tbl):[A-Za-z0-9_-]+)", prose):
            if m.group(1) not in cited:
                cited.append(m.group(1))

        for kind in ("fig", "tbl"):
            seq = [c for c in cited if c.startswith(kind + ":") and c in number]
            jumped: dict[int, list[int]] = {}
            highest = 0
            for label in seq:
                n = number[label]
                if n < highest:
                    jumped.setdefault(highest, []).append(n)
                else:
                    highest = n
            for early, skipped in jumped.items():
                lo, hi = min(skipped), max(skipped)
                span = f"{KIND[kind]} {lo}" if lo == hi else \
                    f"{KIND[kind]}s {lo}–{hi}"
                out.append(
                    f"{name}: {KIND[kind]} {early} (<{label_of[(kind, early)]}>) "
                    f"is cited before {span}; either move its first mention later "
                    f"or move the {KIND[kind].lower()} earlier")
    return out


def check_structure(sources: dict[str, str]) -> tuple[list[str], list[str]]:
    """Checks that need the Typst source rather than the extracted prose.

    A figure or table nobody points to is the one defect here with no honest
    defence: most journals require every one to be cited in the text, in order,
    and a reader who is never sent to a figure will not look at it. So that is an
    error. Undefined acronyms are a warning, because deciding what counts as
    common knowledge in a given field is not something this script can do.
    """
    errors, warnings = [], []
    joined = "\n".join(sources.values())

    labels, refs = {}, set()
    for name, src in sources.items():
        for m in re.finditer(r"<((?:fig|tbl|eq):[A-Za-z0-9_-]+)>", src):
            labels.setdefault(m.group(1), name)
    for m in re.finditer(r"@((?:fig|tbl|eq):[A-Za-z0-9_-]+)", joined):
        refs.add(m.group(1))
    for m in re.finditer(r"#refn\(\s*<((?:fig|tbl|eq):[A-Za-z0-9_-]+)>", joined):
        refs.add(m.group(1))

    for label, where in sorted(labels.items()):
        if label not in refs:
            kind = {"fig": "figure", "tbl": "table", "eq": "equation"}[
                label.split(":")[0]]
            errors.append(f"{where}: {kind} <{label}> is never referenced in the "
                          f"text (most journals require every figure and table "
                          f"to be cited)")

    warnings += check_reference_order(sources)

    # An acronym used more than once but never followed by, or preceded by, a
    # parenthetical expansion anywhere in the manuscript.
    prose = re.sub(r"`[^`]*`", " ", joined)
    counts: dict[str, int] = {}
    # Hyphenated acronyms (DIA-NN, LC-MS) count as one token, not two fragments.
    for m in re.finditer(r"\b[A-Z]{2,}[0-9]*(?:-[A-Z0-9]{2,})*\b", prose):
        counts[m.group(0)] = counts.get(m.group(0), 0) + 1
    for acr, n in sorted(counts.items()):
        if n < 2 or acr.upper() in ACRONYM_OK:
            continue
        defined = (re.search(rf"\(\s*{re.escape(acr)}s?\s*\)", prose)
                   or re.search(rf"\b{re.escape(acr)}s?\s*\([A-Za-z]", prose))
        if not defined:
            warnings.append(f"acronym {acr!r} used {n}x but never expanded")

    return errors, warnings


# --- uno house style: rules that need the raw source ----------------------
#
# Ported from the retired scripts/paper_stats.py, which parsed the whole file and
# therefore needed elaborate exemptions for the ACS back matter. Here they are
# cheaper: prose_check only ever sees the body between the BODY markers plus the
# SI, so the run-in bold labels of Associated Content / Author Information /
# Acknowledgments are out of scope by construction, and table cells go with the
# #figure blocks that contain them.


def _prose_source(src: str) -> str:
    """Raw Typst with the things that are not running prose removed: comments,
    whole figures (caption, table cells and all), and block code / config dumps."""
    src = re.sub(r"(?m)^\s*//.*$", " ", src)
    src = re.sub(r"```.*?```", " ", src, flags=re.S)
    src = readability._strip_balanced(src, "#raw(")
    return readability._strip_balanced(src, "#figure(")


def check_house_style(label: str, src: str) -> tuple[list[str], list[str]]:
    """ACS house style, judged on the raw source because `clean()` has already
    unwrapped the markup these rules are about."""
    errors, warnings = [], []
    body = _prose_source(src)

    # Doubled words, judged before any construct is stripped (see check()).
    for m in re.finditer(r"\b(\w+)\s+\1\b", body, re.I):
        if m.group(1).lower() in {"had", "that"}:   # legitimately doubles
            continue
        errors.append(f"{label}: doubled word {m.group(0)!r}  {_ctx(body, m.start())}")

    # Rule 2: lists. ACS wants running prose in the body; the SI may enumerate.
    for m in re.finditer(r"(?m)^\s*([-+])\s+\S", body):
        kind = "bulleted" if m.group(1) == "-" else "numbered"
        errors.append(f"{label}: {kind} list item in running prose  "
                      f"{_ctx(body, m.start())}")

    # Rule 3: bold. A whole-paragraph bold span is a run-in label and is fine;
    # bold used for emphasis inside a sentence is not.
    for m in re.finditer(readability._markup("*"), body):
        line = body[body.rfind("\n", 0, m.start()) + 1:
                    (body.find("\n", m.end()) + 1 or len(body)) - 1]
        if line.strip() == m.group(0).strip():
            continue                      # run-in label on its own line
        errors.append(f"{label}: bold span {m.group(1)[:40]!r} in running prose  "
                      f"{_ctx(body, m.start())}")

    return errors, warnings


def check_bibliography(sources: dict[str, str]) -> list[str]:
    """Bib entries nothing cites. Typst drops them silently, so they accumulate:
    a reference list that no longer matches what the manuscript argues."""
    bib = HERE / "references.bib"
    if not bib.exists():
        return []
    keys = set(re.findall(r"(?m)^\s*@\w+\s*\{\s*([^,\s]+)\s*,", bib.read_text()))
    # The key stops at a trailing `:`, as Typst's own parser does. Allowing one
    # swallows the colon in `... Klammer et al. @klammer2007: amino-acid counts`
    # and reports a cited entry as an orphan.
    cited = set(re.findall(r"@([A-Za-z][A-Za-z0-9_-]*(?::[A-Za-z0-9_-]+)*)",
                           "\n".join(sources.values())))
    orphans = sorted(k for k in keys if k not in cited)
    return [f"bib entry {k!r} is in references.bib but never cited" for k in orphans]


def inventory(sources: dict[str, str]) -> str:
    """One line of counts, so a change in scale is visible without a diff."""
    joined = "\n".join(sources.values())
    bib = HERE / "references.bib"
    refs = len(re.findall(r"(?m)^\s*@\w+\s*\{", bib.read_text())) if bib.exists() else 0
    xrefs = len(re.findall(r"@(?:fig|tbl|eq|sec):[A-Za-z0-9_-]+", joined))
    # (cross-references are counted in the body + SI only, so a pointer repeated
    # in the back matter does not inflate the count)
    colons = readability.clean(joined).count(":")
    return (f"  inventory  {refs} references, {xrefs} cross-references, "
            f"{colons} colons in prose (informational; keep list and label "
            f"colons, cut rhetorical ones)")


def main() -> int:
    strict = "--strict" in sys.argv
    paper_src = (HERE / "paper.typ").read_text()
    body = readability.slice_body(paper_src)
    si = (HERE / "si-body.typ").read_text()
    targets = {"main": body, "SI": si}

    errors, warnings = [], []
    for label, src in targets.items():
        e, w = check(label, readability.clean(src), readability.clean(no_code(src)))
        errors += e
        warnings += w
        e, w = check_house_style(label, src)
        errors += e
        warnings += w

    e, w = check_structure(targets)
    errors += e
    warnings += w

    # Citations are scanned over the whole file, not just the body: the back
    # matter cites too.
    warnings += check_bibliography({"paper": paper_src, "SI": si})

    for e in errors:
        print(f"  ERROR   {e}")
    for w in warnings:
        print(f"  warn    {w}")

    if not errors and not warnings:
        print("  prose check clean")
    else:
        print(f"\n  {len(errors)} error(s), {len(warnings)} warning(s)")
        print("  rules: STYLE.md   (warnings are judgement calls, not gates)")
    print(inventory({"main": body, "SI": si}))

    return 1 if errors or (strict and warnings) else 0


if __name__ == "__main__":
    raise SystemExit(main())
