repo <- '.'
b <- readRDS(file.path(repo, 'analysis/data/derived_data/component-decomposition',
  'study-b-balanced/study-b-balanced-summary.rds'))
print(names(b))
cols <- intersect(c('analysis', 'N', 'c_bm_pb', 'bias', 'mc_se_bias', 'mean_estimate'),
  names(b))
print(as.data.frame(b[, cols]))
