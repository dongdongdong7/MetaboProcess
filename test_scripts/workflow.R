devtools::document()
# massProcesser demo data
# 0. Environment parameter
# File path
file_dir <- "D:/fudan/Projects/2024/MetaboProcess/Progress/build_package/test_data/massProcesser/"
patterns <- c(".mzXML", ".mzxml", "mzML", ".mzml")
patterns <- paste0(patterns, collapse = "|")
file_path <- list.files(file_dir, pattern = patterns)
file_path <- paste0(file_dir, file_path)
# Environment parameter
os <- get_os()
max_thread <- parallel::detectCores() # 8
# Sample info
pd <- data.frame(sample_name = sub(basename(file_path), pattern = patterns,
                                   replacement = ""),
                 sample_group = c("Case", "Case", "Case",
                                  "Control", "Control", "Control",
                                  "QC", "QC", "QC"),
                 sample_type = c(rep("TEST", 6), rep("QC", 3)),
                 sample_inject = c(2, 6, 8, 4, 7, 9, 1, 3, 5),
                 sample_path = file_path,
                 stringsAsFactors = FALSE)
pd <- dplyr::arrange(pd, sample_inject)
# 1. Peak picking
# 1.1 Load data
data <- MsExperiment::readMsExperiment(pd$sample_path, sampleData = pd)
xcms::spectra(data) %>%
  Spectra::spectraData(c("isolationWindowTargetMz", "isolationWindowLowerMz",
                         "isolationWindowUpperMz", "msLevel", "rtime"))
MsExperiment::sampleData(data)
data_QC <- data[which(pd$sample_type == "QC")]
# 1.2 Optimal parameters
ppmCut <- paramounterPart1(file_path = MsExperiment::sampleData(data_QC)$sample_path, thread = 4)
optPara_xcms <- paramounterPart2(file_path = MsExperiment::sampleData(data_QC)$sample_path, ppmCut = ppmCut, thread = 4)
