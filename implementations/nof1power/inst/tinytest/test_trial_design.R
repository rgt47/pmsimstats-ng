library(tinytest)

#---------------------------------------------------------------------
# build_trial_design returns list(metadata, trialpaths)
#---------------------------------------------------------------------

td <- nof1power::build_trial_design(
  name_longform = "Hybrid",
  name_shortform = "H",
  timepoints = c(2, 4, 6, 8),
  expectancies = c(0.5, 0.5, 0, 0),
  ondrug = list(c(1, 1, 0, 0), c(0, 0, 1, 1))
)

expect_true(is.list(td))
expect_true(all(c("metadata", "trialpaths") %in% names(td)))

# Two paths
expect_equal(length(td$trialpaths), 2)

# Each path is a tibble with one row per timepoint and the expected
# pmsimstats short-form columns.
for (path in td$trialpaths) {
  expect_true(is.data.frame(path))
  expect_equal(nrow(path), 4)
  expect_true(all(
    c("timeptnames", "t_wk", "e", "tod", "tsd", "tpb") %in%
      names(path)
  ))
}

#---------------------------------------------------------------------
# tod and tsd accumulate correctly on path 1 (on-drug first)
#---------------------------------------------------------------------

p1 <- td$trialpaths[[1]]
# Expected tod: 2, 4, 4, 4 (cumulative while on drug; freezes
# at discontinuation).
expect_equal(p1$tod, c(2, 4, 4, 4))
# Expected tsd: 0, 0, 2, 4 (cumulative off-drug after first
# on-drug period).
expect_equal(p1$tsd, c(0, 0, 2, 4))

#---------------------------------------------------------------------
# Metadata round-trips inputs
#---------------------------------------------------------------------

expect_equal(td$metadata$name_longform, "Hybrid")
expect_equal(td$metadata$timepoints, c(2, 4, 6, 8))
expect_equal(td$metadata$expectancies, c(0.5, 0.5, 0, 0))

#---------------------------------------------------------------------
# process_trial_design augments on_drug and t_wk_cumulative
#---------------------------------------------------------------------

aug <- nof1power::process_trial_design(td)
expect_true("on_drug" %in% names(aug))
expect_true("t_wk_cumulative" %in% names(aug))
expect_equal(aug$on_drug, c(TRUE, TRUE, FALSE, FALSE))
