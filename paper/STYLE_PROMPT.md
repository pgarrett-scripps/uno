# Reusable prompt: copy-edit a manuscript into our house style

Paste everything below the line into Claude Code, adjusting the two bracketed
paths at the top. It is written to be self-contained and safe to re-run on an
already-edited paper (it is idempotent: a clean paper produces no changes).

---

Copy-edit `[PATH/TO/paper.typ]` and `[PATH/TO/si-body.typ]` into our house
style. The reference for voice and formatting is
`/home/patrick-garrett/Repos/d_noise/paper/paper.typ` and its `si-body.typ`.
Read a Results section and a Methods section of that paper before you start, and
match how it reads. We are submitting to the Journal of Proteome Research (ACS).

The goal is prose that reads like a careful human scientist wrote it. The current
draft reads like an LLM produced it. Fix that without changing the science.

## Formatting rules, in priority order

1. **Zero em-dashes (—).** Not one may remain in prose. Do not mechanically swap
   in commas. Per instance, either split into two sentences (usually best),
   recast as a comma-bounded clause, or delete the aside if it is filler. Keep
   en-dashes in numeric ranges (`8–11×`, `A–E`), hyphens in compounds
   (`target-decoy`), and minus signs. If a table uses `—` as a "not applicable"
   marker, convert those cells to `n/a` and check the table still fits.

2. **Zero bulleted lists in the body.** ACS papers do not use them. Convert every
   list to running prose with connective tissue ("Two gates bound extension.
   First, ... Second, ..."). Do not just delete the `- ` and leave stubs. Watch
   for text that referenced the list ("the second bullet below") and reword it.
   You may keep bold run-in labels that open a paragraph like a mini-heading
   (`*Noise-floor early-stop.* Extension halts when...`) — that is correct ACS
   style and the reference paper uses it. What must go is the list formatting.

3. **Almost no bold in running prose.** The reference paper bolds *only* figure
   panel labels (`*(a)*`, `*Top left:*`), back-matter headings (`*Data.*`,
   `*Notes*`, author names), table headers/cells, and paragraph run-in labels.
   It never bolds a number, a tool name, or an emphasis word. Strip bold from
   statistics (`*79.9 %*`), emphasis words (`*not*`, `*only*`, `*every*`,
   `*same*`, `*raw*`), tool/dataset names (`*Dinosaur*`, `*PXD070049*`), and
   rhetorical phrases. If a sentence leaned on the bold for contrast, reword so
   word order carries it.

4. **Avoid semicolons.** Replace with a period or a rewrite. They may survive
   only where they separate list items that already contain commas (grant lists,
   version/flag enumerations), inside table cells, and in the ACS address block.

5. **Limit colons.** Keep those introducing a genuine list or a table/figure
   label. Cut rhetorical ones ("The result is stark:", "Two things follow:").

6. **Break up run-on sentences.** Three or more clauses, or past ~40 words,
   means split it. Short declarative sentences are the house style.

7. **Thin dense numbers in the main text.** Do not stack six tools × four
   statistics into one sentence. Keep the headline number in prose and let the
   table carry the rest, adding a pointer to it. **This rule does NOT apply to
   the SI**, whose tables are the canonical home for statistics — keep every
   number there. Before dropping any number, confirm it appears in a table or
   the SI; if you cannot confirm it, keep it.

8. **Reduce parenthetical content.** Long parentheses carrying real content
   become their own sentences. Parentheses holding cross-references, units, or
   short glosses are fine.

9. **Cut redundancy, including across files.** The same fact stated twice in
   different words means keep the stronger one. Check specifically for: a
   mechanism explained in full in both the main text and the SI (the SI should
   defer with "as described in Methods" or a section reference); a method or
   caveat restated in Methods, Results, and Discussion; captions restating the
   body paragraph beside them; captions re-teaching standard statistics (what an
   IQR is, what whiskers mean). Search across the whole document, not
   section-by-section — the worst duplicates span distant sections.

10. **Cut hedging and LLM tics.** "It is worth noting", "importantly",
    "notably", "crucially", "strikingly", "not just X but Y", "we now turn to",
    heavy signposting, and tricolons built for rhythm are usually deletable.

11. **Italicize species binomials** consistently: `_E. coli_`,
    `_S. cerevisiae_`, `_Homo sapiens_`. Not bold, not plain. Apply this inside
    bold run-in labels and table headers too (nesting is fine). Do not touch
    names inside code spans or data labels like `ECOLI`.

## Hard constraints

- **Never change, add, or round a numeric value, unit, percentage, CI, p-value,
  count, version string, accession, grant number, commit hash, seed, file path,
  or URL.** You may drop a number under rule 7; you may never edit one.
- Preserve every `<label>`, `@citation`, `#ref(...)`, `#refn(...)`, `#figure`,
  `table(...)`, `$math$`, backtick code span, and `\u{...}` escape exactly.
- Do not add, delete, reorder, or renumber any figure, table, heading, or
  citation. Do not delete a paragraph carrying unique content. Do not add new
  claims, caveats, or interpretation.
- Do not alter author names, affiliations, emails, ORCIDs, keywords, the ethics
  or AI-use disclosure, or any availability statement beyond punctuation.
- Keep the existing ~76-column hard wrap so the diff stays reviewable.
- This is a copy-edit, not a rewrite. Most sentences should survive with light
  surgery. Leave clean passages alone.

## How to run it

Split the work across parallel subagents to keep it fast, but **give each agent
its own copy of a disjoint chunk in a scratch directory rather than letting
several agents edit one file** — concurrent read-modify-write on the same file
silently clobbers edits. Split on section boundaries, verify the split
reassembles byte-identically first, then reassemble when they finish.

Handle cross-chunk redundancy (rule 9) yourself after reassembly, since no
per-chunk agent can see it.

## Verification before you report back

Write a checker and actually run it. Confirm, for both files:

- The set of `<labels>`, `#ref` targets, and `@citations` is unchanged.
- The count of `#figure` blocks and headings is unchanged.
- **Every numeric token in the edited file existed in the original** (numbers may
  be dropped, never invented or altered).
- Em-dashes in prose: 0. Bullets: 0. Report the surviving semicolons and bold
  spans and justify each category.
- The document compiles (`typst compile paper.typ`), and report the page count
  before and after.

Beware false positives in your own checker: a greedy `[\d,]*` will swallow a
sentence comma into a numeric token, and a citation moved to the end of a
sentence picks up a trailing period. Fix the checker, do not "fix" the paper.

Report any pre-existing factual or numeric inconsistencies you find (a value
that disagrees between the text and its own table, for instance) but **do not
change them** — flag them for me to arbitrate, since that is a data question,
not a copy-edit.
