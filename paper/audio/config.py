#!/usr/bin/env python3
"""Narration configuration. Everything project-specific about the audiobook.

The identity of the manuscript (title, authors, institution, cover wordmark) is
NOT duplicated here. It is read out of ../config.typ, because hardcoding a title
in the narrator is exactly how a finished audiobook ends up announcing a title
the paper no longer has. Only genuinely audio-specific settings live below.

Edit the FILL THIS IN block. The reader at the bottom needs no changes.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

# =============================================================================
# FILL THIS IN
# =============================================================================

# Piper voice. `just audio-setup` resolves this name against piper's own voice
# index and downloads it; browse the options at
# https://huggingface.co/rhasspy/piper-voices
#
# The name alone is enough. There used to be a second setting here holding the
# path under the voices repo ("en/en_US/lessac/medium"), because setup built a
# download URL by hand -- so the two had to be kept in step, and a mismatch
# arrived as a 404 that read like a network problem.
VOICE_NAME = "en_US-lessac-medium"

YEAR = "2026"
GENRE = "Audiobook"

# Blurbs written into the .m4b metadata. Keep them to a couple of sentences.
MAIN_DESC = (
    "Audiobook narration of the uno manuscript: a Rust engine for MS1-only "
    "protein identification and quantitation, benchmarked against ms1searchpy. "
    "Prose only; figures, tables, and block code are omitted. Narrated with "
    "Piper text-to-speech."
)
SI_DESC = (
    "Audiobook narration of the Supporting Information: benchmark protocol, "
    "per-dataset results, entrapment calibration, and the retired measurements "
    "this revision supersedes. Prose sections only (tables and figures "
    "omitted). Narrated with Piper text-to-speech."
)

# Words the voice mangles. Regex -> what to say instead. Applied to the title and
# to the body. A coined name, an acronym meant to be spelled out, and any term
# whose spelling defeats the phonemizer belong here.
#
# The two that matter most in this manuscript are the engine names: "uno" is a
# real English-adjacent word the voice handles, but "ms1searchpy" and "koth" are
# coined and get run together into noise.
PRONUNCIATION: dict[str, str] = {
    r"\bms1searchpy\b": "M S one search pie",
    r"\bbiosaur2\b": "bio saur two",
    r"\bkoth_ff\b": "koth F F",
    r"\bkoth\b": "koth",
    r"\bDeepLC\b": "Deep L C",
    r"\bLightGBM\b": "Light G B M",
    r"\bsequoia-boost\b": "sequoia boost",
    r"\bDirectMS1Quant\b": "Direct M S one Quant",
    r"\bDirectMS1\b": "Direct M S one",
    r"\bMS1\b": "M S one",
    r"\bMS2\b": "M S two",
    r"\bLC-MS\b": "L C mass spec",
    r"\bm/z\b": "m over z",
    r"\bFAIMS\b": "faims",
    r"\bUPS1\b": "U P S one",
    r"\bPXD(\d+)\b": r"P X D \1",
    r"\bq-value\b": "q value",
    r"\bRANSAC\b": "ran sack",
    r"\bmzML\b": "m z M L",
    r"\bFASTA\b": "fasta",
    r"\bTSV\b": "T S V",
    r"\bMBR\b": "M B R",
    r"\btop-\$?K\$?\b": "top K",
}

# Inline math, verbatim as it appears in the source, mapped to spoken English.
# Longest keys are applied first, so a specific expression can override a general
# one. Anything left over has its dollar signs stripped and is read as-is, which
# is usually wrong for real equations, so add every equation that appears in
# running prose. Display equations (the FDP definition) are dropped, not read.
MATH: dict[str, str] = {
    r'$-log_10 P(X >= k | n, p)$': (
        "minus log base ten of the probability that X is at least k, "
        "given n and p"
    ),
    r'$R^2_"rev" = 0.96$': "R squared on reversed sequences of zero point nine six",
    r'$R^2_"fwd" = 0.98$': "R squared on forward sequences of zero point nine eight",
    r'$R^2_"rev"$': "R squared on reversed sequences",
    r'$R^2_"fwd"$': "R squared on forward sequences",
    r"$R^2 = 0.9824$": "R squared of zero point nine eight two four",
    r"$R^2 = 0.980$": "R squared of zero point nine eight zero",
    r"$R^2 = 0.972$": "R squared of zero point nine seven two",
    r"$n_\"yeast\" (k)$": "n yeast of k",
    r"$n_\"human\" (k)$": "n human of k",
    r"$N_\"yeast\"$": "capital N yeast",
    r"$N_\"human\"$": "capital N human",
    r"$(d + 1) \/ t$": "d plus one over t",
    r"$r approx 0.71$": "r of about zero point seven one",
    r"$rho = 0.5$": "rho of zero point five",
    r"$r = 0.68$": "r of zero point six eight",
    r"$q <= 0.01$": "q at most zero point zero one",
    r"$plus.minus 1$": "plus or minus one",
    r"$m slash z$": "m over z",
    r"$< 0.01$": "below zero point zero one",
    r"$+0.71$": "plus zero point seven one",
    r"$+0.25$": "plus zero point two five",
    r"$+0.02$": "plus zero point zero two",
    r"$R^2$": "R squared",
    r"$Delta$": "delta",
    r"$sigma$": "sigma",
    r"$rho$": "rho",
    r"$K$": "K",
    r"$k$": "k",
    r"$n$": "n",
    r"$p$": "p",
    r"$q$": "q",
    r"$t$": "t",
    # Not in the manuscript. These two are exercised by tests/fixture.typ, which
    # is how `just test` proves inline math is still verbalized rather than read
    # as raw Typst. Delete them and the golden narration changes.
    r'$t_"obs" <= t_"max"$': "t observed is at most t max",
    r"$alpha$": "alpha",
}

# Typst symbol tokens. These are generic; extend rather than replace.
SYM: dict[str, str] = {
    r"#sym.minus": "minus ",
    r"#sym.plus.minus": "plus or minus ",
    r"#sym.gt.eq": "greater than or equal to ",
    r"#sym.lt.eq": "less than or equal to ",
    r"#sym.dash.en": "to ",
    r"#sym.tilde": "approximately ",
    r"#sym.space": " ",
    r"#sym.times": "times ",
    r"#sym.arrow.r": "to ",
    r"#sym.arrow.l": "from ",
    r"#sym.arrow.lr": "versus ",
}

# Cover art palette. The default is a dark slate with a cyan accent.
COVER_BG = "#0d1b2a"
COVER_ACCENT = "#38bdf8"
COVER_ACCENT2 = "#5eead4"
COVER_FG = "#f8fafc"
COVER_MUTED = "#94a3b8"

# =============================================================================
# Derived from ../config.typ. Nothing to edit below here.
# =============================================================================

HERE = Path(__file__).resolve().parent
PAPER_DIR = HERE.parent
CONFIG_TYP = PAPER_DIR / "config.typ"
PAPER_TYP = PAPER_DIR / "paper.typ"
SI_TYP = PAPER_DIR / "si-body.typ"

_cfg = CONFIG_TYP.read_text() if CONFIG_TYP.exists() else ""


def typst_str(name: str) -> str:
    """Read a `#let <name> = "..."` binding out of config.typ, or fail loudly."""
    m = re.search(rf'#let\s+{re.escape(name)}\s*=\s*"([^"]*)"', _cfg)
    if not m:
        sys.exit(f'error: could not find `#let {name} = "..."` in {CONFIG_TYP}')
    return m.group(1)


def _author_names() -> list[str]:
    """Every `name: "..."` inside the paper-authors block, in order."""
    m = re.search(r"#let\s+paper-authors\s*=\s*\((.*?)\n\)", _cfg, re.S)
    if not m:
        sys.exit(f"error: could not find `#let paper-authors = (...)` in {CONFIG_TYP}")
    return re.findall(r'name:\s*"([^"]+)"', m.group(1))


def speakable(text: str) -> str:
    """Apply the pronunciation fixes to a piece of text."""
    for pattern, spoken in PRONUNCIATION.items():
        text = re.sub(pattern, spoken, text)
    return text


TITLE = typst_str("paper-title")
WORDMARK = typst_str("paper-wordmark")
COVER_SUBTITLE = typst_str("paper-cover-subtitle").replace("\\n", "\n")
INSTITUTION = typst_str("paper-institution")

AUTHOR_NAMES = _author_names()

# The last token is not the surname when the name carries a generational suffix:
# "John R. Yates III" -> "Yates", not "III". ../config.typ keeps the same list.
NAME_SUFFIXES = {"jr", "sr", "ii", "iii", "iv"}


def surname_of(name: str) -> str:
    parts = name.split()
    while len(parts) > 1 and parts[-1].rstrip(".").lower() in NAME_SUFFIXES:
        parts.pop()
    return parts[-1]


SURNAMES = [surname_of(n) for n in AUTHOR_NAMES]
AUTHOR = ", ".join(SURNAMES)                 # .m4b artist tag
COVER_AUTHORS = "  ·  ".join(SURNAMES)  # "Lovelace · Hopper"

SI_TITLE = f"{TITLE}: Supporting Information"


def spoken_title() -> str:
    """The title read aloud: pronunciation fixes applied, and punctuation that no
    voice renders as a pause turned into a sentence break."""
    title = speakable(TITLE)
    title = title.replace(":", ".").replace("–", " ").replace("-", " ")
    return re.sub(r"\s+", " ", title).strip().rstrip(".") + "."
