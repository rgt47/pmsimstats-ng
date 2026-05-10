# Integration Test: Report Rendering
# Tests that research reports can be rendered successfully

library(tinytest)
library(here)

# Manuscript Rmd file exists and is valid
manuscript_path <- here("analysis", "report", "manuscript.Rmd")
expect_true(file.exists(manuscript_path))

manuscript_content <- readLines(manuscript_path)
expect_true(length(manuscript_content) > 0)

expect_true(any(grepl("^---$", manuscript_content)))

# Report dependencies are available
required_packages <- c("rmarkdown", "knitr")

for (pkg in required_packages) {
  expect_true(requireNamespace(pkg, quietly = TRUE),
              info = paste("Package", pkg, "is required for report rendering"))
}

# Bibliography files exist
bib_path <- here("analysis", "report", "nof1-pgt.bib")
expect_true(file.exists(bib_path))

alt_bib <- here("analysis", "report", "references.bib")
expect_true(file.exists(alt_bib))

csl_path <- here("analysis", "report", "statistics-in-medicine.csl")
expect_true(file.exists(csl_path))

# Vignette documents exist
vignettes <- c(
  "monte_carlo_design_comparison.Rmd",
  "carryover_only.Rmd",
  "carryover_sensitivity_analysis.Rmd"
)

for (vig in vignettes) {
  vig_path <- here("analysis", "report", vig)
  expect_true(file.exists(vig_path), info = paste("Missing vignette:", vig))
}

# Manuscript can be parsed without errors
manuscript_path <- here("analysis", "report", "manuscript.Rmd")

yaml_content <- rmarkdown::yaml_front_matter(manuscript_path)
expect_true(is.list(yaml_content))

# Note: Actual rendering test is commented out to avoid LaTeX dependencies
# in CI. Uncomment for local testing if LaTeX is available.
# manuscript_path <- here("analysis", "report", "manuscript.Rmd")
# output_dir <- here("analysis", "report")
#
# rmarkdown::render(manuscript_path,
#                   output_dir = output_dir,
#                   quiet = TRUE)
#
# pdf_path <- here("analysis", "report", "manuscript.pdf")
# expect_true(file.exists(pdf_path))
