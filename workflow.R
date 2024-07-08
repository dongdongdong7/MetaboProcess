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
data_QC <- data[which(pd$sample_type == "QC")]
# CentWave parameter optimization
# (1) Noise estimation
#chrDf_bins <- generateBin(ndata = ndata, bin = 0.05, mslevel = 1, thread = 24)
load("D:/fudan/Projects/2024/MetaboProcess/Progress/build_package/test_data/massProcesser/chrDf_bins.RData")
load("D:/fudan/Projects/2024/MetaboProcess/Progress/build_package/test_data/massProcesser/chrDf_ZOIs.RData")
load("D:/fudan/Projects/2024/MetaboProcess/Progress/build_package/test_data/massProcesser/chrDf_ZOIs_2.RData")
load("D:/fudan/Projects/2024/MetaboProcess/Progress/build_package/test_data/massProcesser/chrDf_ZOIs_3.RData")

scanNumberVec <- sapply(chrDf_bins, nrow)
chrDf_bins <- chrDf_bins[which(scanNumberVec > runInfo$scanCount * 0.05)]
chrDf_ZOIs <- bins2ZOIs(chrDfList = chrDf_bins, thread = 4, sn = 0, IETH = 10)
chrDf_ZOIs_2 <- lapply(1:length(chrDf_ZOIs), function(i) {
  print(i)
  x <- chrDf_ZOIs[[i]]
  chrDf_tmp <- tryCatch({
    edgeTrack_finer(x, tol_m = 45)
  }, error = function(e){
    edgeTrack_crude(x, tol_m = 40)
  })
})
plot_chrDf(chrDf_ZOIs_2[[1131]])
massTolVec <- sapply(chrDf_ZOIs_2, cal_massTol)
plot(x = 1:length(massTolVec), y = massTolVec)
boxplot(massTolVec)
s <- boxplot.stats(massTolVec)$stats
chrDf_ZOIs_2 <- chrDf_ZOIs_2[which(massTolVec<s[5])]
massTolVec <- massTolVec[which(massTolVec<s[5])]
plot(x = 1:length(massTolVec), y = massTolVec)
plot_chrDf(chrDf_ZOIs_2[[933]])
tol_ppm <- max(massTolVec)
cal_cwtParam(chrDfList = chrDf_ZOIs_2)



i <- 2186
plot_chrDf(chrDf_ZOIs[[i]], baseline = TRUE)
chrDf_i <- edgeTrack_finer(chrDf_ZOIs[[i]], tol_m = 45)
plot_chrDf(chrDf_i, baseline = TRUE)
chrDf_i <- edgeTrack_crude(chrDf_ZOIs[[i]], tol_m = 40)

i <- 96 # 1021 1022 1141 5141 5142 5143 5144 5146 10002 12009
chrDf <- chrDf_bins[[i]]
chrDf$intensity <- smoothMean(chrDf$intensity, size = 3)
#chrDf$intensity <- smoothSg(chrDf$intensity)
chrDf_noise <- noiseEstimation(chrDf)
plot_chrDf(chrDf, linewidth = 1,
           noise = chrDf_noise,
           xlim = NA)
ZOIList <- FindZOI(chrDf = chrDf, noise = chrDf_noise, etlD = 1, IETH = 10)
length(ZOIList)
j <- 2
chrDf_j <- ZOIList[[j]]
#chrDf_j$intensity <- smoothMean(chrDf_j$intensity, size = 3)
plot_chrDf(chrDf_j, linewidth = 1,
           noise = NA,
           xlim = NA)
chrDf_j <- baselineEs(chrDf_j, loops = 8, tol_m = 30)
plot_chrDf(chrDf_j, linewidth = 1,
           noise = NA,
           xlim = NA,
           baseline = TRUE)
#cal_IE(chrDf_j$intensity)
chrDf_j_new <- edgeTrack(chrDf = chrDf_j)
plot_chrDf(chrDf_j_new, linewidth = 1,
           noise = NA,
           xlim = NA,
           baseline = TRUE)
#cal_IE(chrDf_j_new$intensity)

# Paramounter demo data
# 0. Environment demo data
# File path
file_dir <- "D:/fudan/met/Paramounter/test1/ParamounterFiles/"
patterns <- c(".mzXML", ".mzxml", "mzML", ".mzml")
patterns <- paste0(patterns, collapse = "|")
file_path <- list.files(file_dir, pattern = patterns)
file_path <- paste0(file_dir, file_path)
# Environment parameter
os <- get_os()
max_thread <- parallel::detectCores() # 8
fileh <- mzR::openMSfile(file_path[1], backend = NULL)
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
                 sample_group = c("QC", "QC", "QC", "QC", "QC"),
                 sample_type = c(rep("QC", 5)),
                 sample_inject = c(1,2,3,4,5),
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
optParam <- optParam4xcms(data_QC = data_QC, res_dir = "D:/fudan/Projects/2024/MetaboProcess/Progress/build_package/test_data/Paramounter/urine/",
                          output_bins = TRUE, output_ZOI = TRUE, output_ZOI2 = TRUE, output_ZOI3 = TRUE,
                          thread1 = 4, thread2 = 4)
load("D:/fudan/Projects/2024/MetaboProcess/Progress/build_package/test_data/Paramounter/urine/chrDf_bins.RData")
load("D:/fudan/Projects/2024/MetaboProcess/Progress/build_package/test_data/Paramounter/urine/chrDf_ZOIs.RData")
load("D:/fudan/Projects/2024/MetaboProcess/Progress/build_package/test_data/Paramounter/urine/chrDf_ZOIs_2.RData")
load("D:/fudan/Projects/2024/MetaboProcess/Progress/build_package/test_data/Paramounter/urine/chrDf_ZOIs_3.RData")
i <- 9486
plot_chrDf(chrDf_ZOIs_3[[i]], baseline = TRUE)
cal_massTol(chrDf_ZOIs_3[[i]])
massTolVec <- sapply(chrDf_ZOIs_2, function(x) {cal_massTol(chrDf = x, factor = 1, range = 1)})
peakWidthVec <- sapply(1:length(chrDf_ZOIs_2), function(i) {
  chrDf <- chrDf_ZOIs_2[[i]]
  peakWidth <- max(chrDf$rt) - min(chrDf$rt)
  return(peakWidth)
})
massTol_s <- boxplot.stats(massTolVec[massTolVec < 100])$stats[5]
peakWidth_s_max <- boxplot.stats(peakWidthVec[peakWidthVec < 100 & peakWidthVec > 1])$stats[5]
peakWidth_s_min <- boxplot.stats(peakWidthVec[peakWidthVec < 100 & peakWidthVec > 1])$stats[1]
chrDf_ZOIs_3 <- chrDf_ZOIs_2[which(massTolVec<massTol_s & peakWidthVec<100 & peakWidthVec>1)]
massTolVec <- massTolVec[which(massTolVec < massTol_s &
                                 peakWidthVec < 100 & peakWidthVec > 1)]
tol_ppm <- sort(massTolVec)[round(length(massTolVec) * 1)]
cwtParam <- cal_cwtParam(chrDfList = chrDf_ZOIs_3, widthRange = 1)
massTolVec <- sapply(chrDf_ZOIs_3, function(x) {cal_massTol(chrDf = x, factor = factor, range = range)})
i <- 1
plot_chrDf(chrDf_ZOIs_3[[which(peakWidthVec3 > 90)[i]]], baseline = TRUE)
cal_massTol(chrDf_ZOIs_3[[which(massTolVec > 50)[i]]], factor = 2, range = 1)
