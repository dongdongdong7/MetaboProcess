devtools::document()
# ****** Test for peak picking ******
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
fileh <- mzR::openMSfile(file_path[2], backend = NULL)
instrumentInfo <- mzR::instrumentInfo(fileh)
runInfo <- mzR::runInfo(fileh)
header <- mzR::header(fileh)
mzR::close(fileh);rm(fileh)
# 1. Load data
# Whether centroid
checkCentroid(file_path = file_path)
# Prepare sample Data
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
# Load data
data <- MsExperiment::readMsExperiment(pd$sample_path, sampleData = pd)
xcms::spectra(data) %>%
  Spectra::spectraData(c("isolationWindowTargetMz", "isolationWindowLowerMz",
                "isolationWindowUpperMz", "msLevel", "rtime"))
MsExperiment::sampleData(data)
# CentWave parameter optimization
# (1) Noise estimation
data_QC <- data[which(pd$sample_type == "QC")]
MsExperiment::sampleData(data_QC)
xcms::spectra(data_QC)
n <- 1
ndata <- data_QC[n]
#chrDf_bins <- generateBin(ndata = ndata, bin = 0.05, mslevel = 1, thread = 24)
load("D:/fudan/Projects/2024/MetaboProcess/Progress/build_package/test_data/massProcesser/chrDf_bins.RData")
i <- 1021 # 1021 1022 1141 5141 5142 5143 5144 5146 10002 12009
chrDf <- chrDf_bins[[i]]
chrDf$intensity <- smoothMean(chrDf$intensity, size = 3)
chrDf_noise <- noiseEstimation(chrDf)
#chrDf$intensity <- smoothMean(chrDf$intensity, size = 3)
#chrDf_noise <- noiseEstimation2(chrDf)
#chrDf$intensity <- signal::sgolayfilt(chrDf$intensity, p = 1, n = 5)
plot_chrDf(chrDf, linewidth = 1,
           noise = chrDf_noise,
           xlim = NA)
ZOIList <- FindZOI(chrDf = chrDf, noise = noiseEstimation(chrDf), etlD = 1, IETH = 10)
length(ZOIList)
j <- 3
chrDf_j <- ZOIList[[j]]
#chrDf_j$intensity <- smoothMean(chrDf_j$intensity, size = 3)
p <- plot_chrDf(chrDf_j, linewidth = 1,
           noise = NA,
           xlim = NA)
p
chrDf_j <- baselineEs(chrDf_j, loops = 8, tol_m = 30)
p <- plot_chrDf(chrDf_j, linewidth = 1,
                noise = NA,
                xlim = NA,
                baseline = TRUE)
p
cal_IE(chrDf_j$intensity)
chrDf_j_new <- edgeTrack(chrDf = chrDf_j)
plot_chrDf(chrDf_j_new, linewidth = 1,
           noise = NA,
           xlim = NA,
           baseline = TRUE)
cal_IE(chrDf_j_new$intensity)
xcms::spectra(ndata)
fileh <- mzR::openMSfile(MsExperiment::sampleData(ndata)$sample_path, backend = NULL)
runInfo <- mzR::runInfo(fileh)
header <- mzR::header(fileh)
mzR::close(fileh);rm(fileh)
seq(runInfo$lowMz, runInfo$highMz, 0.05)
