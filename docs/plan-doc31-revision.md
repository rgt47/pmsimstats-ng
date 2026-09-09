# A Plan for the Revision of Document 31

*2026-09-09 09:35 PDT*

Author: pmsimstats team

## 1. Purpose

Document 31 was assembled on 2026-09-09 by merging the component
decomposition pedagogy of document 24 with the hybrid design material
of document 28, and by adding a new Part IV reporting papers 13 and
14. The assembly was mechanical. We concatenated sections, demoted
their headings by one level, and wrote the new part.

Unfortunately, a mechanical merge produces a mechanical result. The
document is complete in the sense that nothing was lost, but it is
not yet coherent as a single piece of writing. We review it here
against four criteria, namely logical consistency, redundancy,
notation, and readability, and we propose a revision in five phases.

We should state one conclusion at the outset, since it governs
everything that follows. The decision to retain document 24 as a
standalone and maintained text, taken after document 31 was
assembled, has made the larger part of document 31 redundant. Any
revision that does not address this will be repairing the trim on a
building whose foundation has shifted.

## 2. The state of the document

Document 31 runs to 2,018 lines and 14,035 words. Of these, roughly
8,200 words are document 24 reproduced verbatim and roughly 3,600 are
document 28. The new material of Part IV accounts for something under
2,300 words, or about one sixth of the whole.

## 3. Logical consistency

### 3.1 Superseded guidance is left standing

The document gives advice in three places. Document 28 contributed a
section headed 'The conclusion drawn from the 06 analysis' and
another headed 'Advice: how and when to use the hybrid design'. Part
IV then contributes 'Revised guidance', which revises them.

Neither of the earlier sections carries any indication that it has
been revised. A reader working through the document in order will
encounter the older advice first, in a section that reads as
authoritative, and will meet the revision several hundred lines
later. That is problematic. The revision is not a refinement of the
earlier advice; on the question of which components a simulation must
carry, it reaches a different conclusion.

### 3.2 The conflict with Part I is acknowledged in one direction only

Part IV, section 14, states plainly that paper 14 does what document
24 warns against. We regard that as correctly handled. The reciprocal
signpost is missing, however. Document 24's section 'Can we assume no
TV effect?' now sits in Part I, some 1,100 lines before the passage
that qualifies it, with no forward reference.

### 3.3 The frequently asked questions are stale

Part V reproduces document 24's question and answer material without
revision. Two of the questions bear directly on what Part IV
establishes. The question 'Why is the biomarker only correlated with
BR, not with PB?' is answered in terms of modeling intent. Part IV
shows that the answer has a structural component as well, since the
implementation provides no biomarker to natural-history parameter at
all. The question 'How many timepoints do I need to identify all
three components?' is answered without reference to the
positive-definiteness ceiling that Part IV derives.

### 3.4 Cross-references inherited from the sources are dangling

We count eight references of the form 'Section 3.4' or 'Section 4.1'
which were written against document 28's own numbering. That
numbering did not survive the merge. One of these, at line 1582, was
written by us in Part IV and refers to a heading that the merge had
already renumbered.

## 4. Redundancy

### 4.1 The document duplicates a living text

This is the central difficulty. Parts I through III reproduce
document 24 in full, and document 24 is now to be maintained as the
canonical pedagogical treatment. Any future correction to the
components material must therefore be applied twice, and the two
copies will drift.

We see two ways forward, and the choice between them is a genuine
trade-off rather than a matter of taste.

On the one hand, document 31 could be made self-contained, on the
grounds that a reader deciding which components to carry should not
have to consult a second file. The cost is a permanent maintenance
obligation and the near certainty that the copies will diverge.

On the other hand, document 31 could cite document 24 and reproduce
only what Part IV actually depends upon. The cost is that the reader
must have document 24 to hand. The benefit is a single source of
truth and a document of roughly 5,800 words rather than 14,000.

In our view the second is clearly preferable. The audiences differ.
Document 24 is read by someone learning the decomposition, and
document 31 by someone deciding what a particular study must carry.
The second reader is not the same person as the first and is unlikely
to want 8,000 words of introduction before reaching the guidance.

### 4.2 Two bibliographies and two limitations sections

Document 28's references section is stranded in the middle of Part
III, where it reads as though the document has ended. Document 24's
further reading sits in Part V. Similarly, document 28's 'Limitations
and open items' appears at line 1456 and Part IV's 'Open problems' at
line 1779, and the two overlap without either acknowledging the
other.

## 5. Notation

The two source documents adopted different conventions and the merge
preserved both. We find the following inconsistencies.

The moderation parameter appears as `c_bm,PB` nine times, as
`c.bm.pb` three times, and as `c_bm` elsewhere. The first is
mathematical notation, the second is the name of an R argument, and
the two are not the same object. The distinction is worth preserving,
but it needs to be made deliberately rather than by accident of
provenance.

The expectancy weight appears sixteen times as `eta` in a code font,
inherited from document 28, and elsewhere as a LaTeX symbol. The
natural history maximum appears nine times as `m_{TV}` and twice as
`m_TV`.

More generally, document 24 sets mathematics in LaTeX delimiters and
document 28 sets it in code fences. Part IV uses both. We count
ninety-one LaTeX expressions against three indented display blocks.

### 5.1 A proposed convention

We propose the following, and we note that it is a convention rather
than a discovery, so the particular choices matter less than their
consistent application.

a) Mathematical quantities take LaTeX delimiters. The moderation
   parameter is `$c_{bm}$`, the biomarker to placebo coupling is
   `$c_{bm,PB}$`, the expectancy weight is `$\eta$`, and the
   component maxima are `$m_{TV}$` and `$m_{PB}$`.

b) Software identifiers take a code font. The R argument is
   `c.bm.pb`, the components argument is `components`, and the
   function is `buildSigma()`.

c) Displayed mathematics uses LaTeX display delimiters throughout.
   The three indented plain text blocks in Part IV are converted.

d) Component names are given in full at first use in each part, and
   by initials thereafter. We adopt 'natural history (TV)' rather
   than the several variants now present.

## 6. Readability

The heading structure does not survive inspection. Part IV is
numbered from twelve to seventeen, and there are no sections one
through eleven. Document 28's original numbers, two through seven,
persist as orphans under Parts II and III, where they refer to a
scheme that no longer exists. Parts I through III contain no second
level headings at all, so the reader descends from a part title
directly to a third level heading, while Part IV uses both levels.
There are forty three fourth level headings.

The effect is that a reader cannot tell depth from typography. We
recommend a single scheme, described in phase 3 below.

## 7. The plan

We propose five phases. The first is a decision rather than work, and
it governs the scope of the remainder.

### Phase 0. Decide the scope question

The question of section 4.1 must be settled first. Is document 31 to
be self-contained, or is it to cite document 24?

Our recommendation is that it cite. Under that decision, Parts I
through III are replaced by a substantially shorter Part I which
states only what Part IV depends upon, namely the three components in
brief, the identifying contrasts, the analysis model actually fitted,
and the expectancy identification limit. We estimate this at
1,200 words against the present 11,800.

The remaining phases are written on that assumption. Should the
opposite decision be taken, phases 1 and 3 grow considerably and
phase 2 becomes more important rather than less.

### Phase 1. Restructure

Replace Parts I through III with a compact Part I as above. Retain
document 28's material on the hybrid design and the expectancy
identification limit, since Part IV depends on both and document 28
is genuinely superseded. Promote Part IV to Part II and Part V to
Part III.

Fold document 28's limitations into Part IV's open problems, and
merge the two bibliographies into a single references section at the
end.

### Phase 2. Repair the logic

Add a forward reference from the retained discussion of natural
history to the passage in Part IV that qualifies it. Mark document
28's advice as revised, with a pointer, or remove it if phase 1 has
already displaced it. Revise the two stale questions in the reference
part so that they reflect what Part IV establishes. Repair or remove
the eight dangling cross-references.

### Phase 3. Impose the notation and heading conventions

Apply the convention of section 5.1 throughout. Adopt a single
heading scheme, in which parts are numbered with Roman numerals,
sections within a part are numbered continuously from one, and no
heading below the third level is used except within the reference
part.

### Phase 4. Read for prose

The merged text has three authorial registers in it, since documents
24 and 28 were written at different times for different purposes and
Part IV was written yesterday. A read for consistency of voice is
warranted once the structure is settled. This is the least urgent of
the phases and should not be attempted before phase 3 is complete.

## 8. Estimated effect

If the phases are carried out as described, document 31 falls from
14,035 words to something near 5,800. The new material of Part IV is
untouched in substance. The pedagogy remains available, in one place,
in document 24.

## 9. Conclusions

In conclusion, four points are to be emphasized.

First, the merge was performed before the decision to retain document
24 as a standalone text, and that decision has left the larger part
of document 31 redundant. The scope question of phase 0 should be
settled before any other work is undertaken.

Second, the logical difficulties are real but small. Superseded
advice stands unmarked, two questions in the reference material are
now stale, and eight cross-references dangle. None of these requires
new analysis to fix.

Third, the notational inconsistencies are entirely a consequence of
merging two documents that had each been internally consistent. They
are straightforward to repair once a convention is chosen, and we
have proposed one.

Fourth, and we state this as a caution rather than a
recommendation, the substance of Part IV rests in part on results
that are not yet final. Paper 13's simulation was still running when
this plan was written, and its reported figures are those of a
superseded run. Any editorial pass over Part IV should wait until
those numbers have settled.
