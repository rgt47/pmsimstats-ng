## build-summary.R
##
## Builds summary.Rmd from summary-template.Rmd by replacing each
## <<chunk:LABEL>> line with the chunk of that label copied verbatim
## from report.Rmd, so the summary's tables and numbers cannot drift
## from the paper. Section cross-references (\ref{sec-...} and
## \@ref(sec-...)) become "the full report", because the summary has
## no numbered sections. Run from the repository root:
##   Rscript analysis/report/02-carryover-sensitivity/build-summary.R
dir <- 'analysis/report/02-carryover-sensitivity'
rep <- readLines(file.path(dir, 'report.Rmd'))
tpl <- readLines(file.path(dir, 'summary-template.Rmd'))

starts <- grep('^```\\{r [A-Za-z0-9-]+', rep)
labels <- sub('^```\\{r ([A-Za-z0-9-]+).*', '\\1', rep[starts])
stopifnot(!anyDuplicated(labels))
fences <- grep('^```\\s*$', rep)
chunk <- function(lab) {
  i <- starts[labels == lab]
  if (length(i) != 1) stop('no chunk labelled ', lab, ' in report.Rmd')
  j <- min(fences[fences > i])
  rep[i:j]
}

out <- unlist(lapply(tpl, function(line) {
  m <- regmatches(line, regexec('^<<chunk:([A-Za-z0-9-]+)>>$', line))[[1]]
  if (length(m) == 2) chunk(m[2]) else line
}))
## Joined, because a reference can break across lines.
txt <- paste(out, collapse = '\n')
txt <- gsub('Sections?\\s+\\\\\\\\ref\\{sec-[^}]+\\}(\\s+and\\s+\\\\\\\\ref\\{sec-[^}]+\\})?',
            'the full report', txt)
txt <- gsub('Sections?\\s+\\\\@ref\\(sec-[^)]+\\)(\\s+and\\s+\\\\@ref\\(sec-[^)]+\\))?',
            'the full report', txt)
stopifnot(!grepl('ref[({]sec-', txt))
writeLines(txt, file.path(dir, 'summary.Rmd'))
cat(sprintf('summary.Rmd: %d lines, %d chunks copied from report.Rmd\n',
            length(out), sum(grepl('^<<chunk:', tpl))))
