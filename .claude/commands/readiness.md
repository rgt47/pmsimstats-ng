---
description: Grade every manuscript under analysis/report/ on readiness for journal submission (A-E), checking each paper's latest referee findings against later commits; writes a dated docs/NN readiness assessment in the docs/66 format
argument-hint: '[paper numbers, e.g. 01 02 15] [--referee]'
---

# Submission readiness assessment

Produce a paper-by-paper statement of how far each manuscript under
`analysis/report/NN-*/` is from journal submission, what blocks it,
and what it would cost to clear. This repeats the method of
`docs/66-manuscript-readiness-assessment.md` (2026-10-29). Read that
document first: its structure, grade scale and tone are the template.

Arguments: `$ARGUMENTS`

- Paper numbers (for example `01 02 15`) restrict the assessment to
  those papers. With none, assess every numbered folder.
- `--referee` runs a full referee pass (the `review-manuscript`
  procedure) on any paper whose evidence base is stale (step 2).
  Without it, stale papers are graded from the evidence available and
  flagged as needing a re-referee. A referee pass is expensive; never
  run one without this flag.

## What this assessment does not do

- It does not edit any manuscript, driver or data file.
- It does not rerun simulations and does not render manuscripts.
- It does not commit. Report the output file and ask.

## Step 1. Inventory

1. Get the date: `date "+%Y-%m-%d %H:%M %Z"`.
2. List the paper folders: `ls -d analysis/report/[0-9][0-9]-*/`.
3. Find the most recent previous assessment:
   `ls docs/ | grep -i readiness`. Its grades are the baseline for
   the change column.
4. For each paper, record:
   - the manuscript source (`report.Rmd`, or the main `.Rmd`/`.tex`);
   - the newest `referee-report-*.md` in the folder and its date;
   - any remediation log, README, `title-page*`, `cover-letter*`;
   - commits touching the folder since that referee report:
     `git log --oneline --since=<report date> -- <folder>`.

## Step 2. Classify the evidence base for each paper

- **Current**: a referee report exists and the manuscript has not
  been rewritten since (commits after it are targeted fixes).
- **Stale**: no referee report, or the manuscript was substantially
  rewritten after the latest one (a restructure, a new framing, a
  rerun that replaced the results). Say which.

## Step 3. Assess each paper

Delegate to subagents (`Agent`, `subagent_type: general-purpose`),
about three papers per agent so that a full run uses five agents.
Launch them in one message. Give each agent this brief, filled in:

> Assess the submission readiness of papers <list> in
> /Users/zenn/prj/res/36-pmsimstats-ng/pmsimstats-ng, using the
> method of docs/66-manuscript-readiness-assessment.md (read it).
> For each paper:
>
> 1. Read the manuscript source in full and its newest referee report
>    (<path>, <date>).
> 2. List every major finding of that report. For each, inspect the
>    commits since the report (`git log -p --since=<date> -- <folder>`)
>    and the current source. Mark it **addressed** only if the fix is
>    visible in the source or a commit diff, **partly open**, **open**,
>    or **disclosed only** (acknowledged in the text, not fixed). Cite
>    commit hashes.
> 3. Check for blockers: `RESULTS-PENDING` or `TODO` markers,
>    unresolved references, hard-coded numbers that no longer match
>    stored results (spot-check the headline numbers against the files
>    the chunks read), data files the Rmd reads that no longer exist,
>    stale DGP settings relative to the refined process in CLAUDE.md
>    (c_1 = 0.2, anchored carryover, stop on non-PD), and citations of
>    unpublished companion papers.
> 4. Check front matter and side files: target journal (and whether
>    the cover letter, title page and README agree on it and on the
>    title), author roster, declarations, data and code availability,
>    a code DOI, word count against the target journal's limit.
> 5. Propose a grade on the scale below and an effort estimate, and
>    say what would move the grade.
>
> Do not edit, rerun or render anything. Label every claim verified
> (executed and observed), inspected (read in source or diff) or
> inferred. Return, per paper: short title, target journal, evidence
> base (current or stale, and why), the findings table with status,
> other blockers, proposed grade, effort, and one paragraph of
> assessment in the voice of docs/66.

If `--referee` was given, run the `review-manuscript` procedure on
each stale paper first (one agent per paper; it writes a dated
`referee-report-YYYY-MM-DD.md` in the paper's folder), then assess
that paper against the new report.

When the agents return, check their grades against the scale and
against one another. Reconcile any grade that the agent's own
findings do not support, and say so in the text.

## Grade scale

- **A** submittable after copy-edit.
- **B** minor revision.
- **C** major revision, no new simulation.
- **D** needs new or rerun simulation.
- **E** incomplete draft (results section empty or pending).

A grade in parentheses, as in C (D), applies if an optional or
contested analysis turns out to be required. Grades are the
assessor's judgment and must be stated as such.

## Step 4. Cross-paper synthesis

- **Common to every paper**: blockers shared by all (authors,
  declarations, code DOI, and so on).
- **Cross-cutting findings**: inconsistent settings across papers,
  repeated arguments, companion citations to incomplete papers, stale
  side files. Mark each as checked or inferred.
- **Recommended consolidation**: merges or folds suggested by overlap
  in design, engine or argument. State that these are proposals.
- **Suggested order of work**: a numbered list, starting with the
  papers closest to submission.

## Step 5. Write the document

Write `docs/NN-manuscript-readiness-assessment-YYYY-MM-DD.md`, where
`NN` is the next free number in `docs/`. Follow docs/66:

- the YAML header of docs/66 (geometry, 10pt, line spread);
- title `Submission readiness of the N manuscripts`, then the
  timestamp line, then `**Author.** pmsimstats team`;
- **Purpose** (name the assessment this one supersedes);
- **Method and evidence** (which papers had a current referee report,
  which were stale, whether `--referee` ran, what was not done);
- the grade scale;
- a summary table:
  `| Paper | Short title | Target | Grade | Previous | Effort to submit |`;
- the paper-by-paper assessment, then the step 4 sections;
- **Limitations of this assessment** last.

Style: `analysis/report/STYLE.md` and the project CLAUDE.md (plain
declarative sentences, US English, no em dashes, 78-character lines;
the studies are exploratory, never pre-registered). Run
`perl tools/notation-lint.pl` on the file, then render it from inside
`docs/` with `bash ../tools/render.sh <file>`.

## Step 6. Report

Give the user the summary table, the grades that changed since the
previous assessment and why, the papers flagged as needing a
re-referee, and the path of the new document. List what was not
checked. Do not commit unless asked.
