# package loads correctly and exports its core simulation functions
expect_true(
  exists("buildtrialdesign", mode = "function"),
  info = "buildtrialdesign() should be exported by pmsimstats"
)
expect_true(
  exists("generateData", mode = "function"),
  info = "generateData() should be exported by pmsimstats"
)
expect_true(
  exists("generateSimulatedResults", mode = "function"),
  info = "generateSimulatedResults() should be exported by pmsimstats"
)
