load("D:/fudan/Projects/2024/MetaboDeconv/Progress/build_package/MetaboDeconv/test_data/QC/data_MSE.RData")
data_MSE
# We just need chrom peak firstly.
data_MSE <- xcms::dropAdjustedRtime(data_MSE)
