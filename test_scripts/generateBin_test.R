# Add slide parameter
library(magrittr)
file_dir <- "D:/fudan/Projects/2024/MetaboProcess/Progress/build_package/test_data/Paramounter/urine_demo/"
patterns <- c(".mzXML", ".mzxml", "mzML", ".mzml")
patterns <- paste0(patterns, collapse = "|")
file_path <- list.files(file_dir, pattern = patterns)
file_path <- paste0(file_dir, file_path)
file_path <- file_path[5]
pd <- data.frame(sample_name = sub(basename(file_path), pattern = patterns,
                                   replacement = ""),
                 sample_group = c("QC"),
                 sample_type = c("QC"),
                 sample_inject = c(1),
                 sample_path = file_path,
                 stringsAsFactors = FALSE)
pd <- dplyr::arrange(pd, sample_inject)
# Load data
ndata <- MsExperiment::readMsExperiment(pd$sample_path, sampleData = pd)
generateBin <- function(ndata, bin = 0.1, slide = 0.05, mslevel = 1, thread = 1){
  sps <- xcms::spectra(ndata) %>% Spectra::filterMsLevel(mslevel)
  raw <- MSnbase::readMSData(MsExperiment::sampleData(ndata)$sample_path)
  # raw <- xcms::xcmsRaw(MsExperiment::sampleData(ndata)$sample_path)
  # fileh <- mzR::openMSfile(MsExperiment::sampleData(ndata)$sample_path, backend = NULL)
  # runInfo <- mzR::runInfo(fileh)
  # mzR::close(fileh);rm(fileh)
  # if(is.na(lowMz)) lowMz <- runInfo$lowMz
  # if(is.na(highMz)) highMz <- runInfo$highMz
  lowMz <- round(min(sapply(Spectra::mz(sps), min)))
  highMz <- round(max(sapply(Spectra::mz(sps), max)))
  rtSps <- Spectra::rtime(sps)
  peaksData <- Spectra::peaksData(sps)
  mzS <- seq(lowMz, highMz - bin, by = slide)
  mzE <- seq(lowMz + bin, highMz, by = slide)
  maxN <- max(c(length(mzS), length(mzE)))
  #mzRange <-seq(lowMz, highMz, by = bin)
  #mzIter <- lapply(1:(length(mzRange) - 1), function(i) {c(mzRange[i], mzRange[i + 1])})
  mzIter <- lapply(1:maxN, function(i) {c(mzS[i], mzE[i])})
  #browser()
  loop <- function(x){
    mz_range <- x
    chrDf <- purrr::list_rbind(lapply(1:length(peaksData), function(j) {
      if(is.matrix(peaksData[[j]])) peakMat <- peaksData[[j]]
      else if(is.vector(peaksData[[j]])){
        peakMat <- matrix(peaksData[[j]], ncol = 2, dimnames = list(NULL, c("mz", "intensity")))
        warnings("peaksData is a vector!")
      }else{
        stop("peaksData is wrong!")
      }
      idx <- which(peakMat[, "mz"] >= mz_range[1] & peakMat[, "mz"] <= mz_range[2])
      if(length(idx) == 0) return(NULL)
      #mz <- mean(peakMat[idx,"mz"])
      mz <- peakMat[idx, "mz"][which.max(peakMat[idx,"intensity"])]
      intensity <- max(peakMat[idx,"intensity"])
      peakDf <- dplyr::tibble(mz = mz, intensity = intensity)
      peakDf$rt <- rtSps[j]
      return(peakDf)
    }))
    return(chrDf)
  }
  pb <- utils::txtProgressBar(max = length(mzIter), style = 3)
  if(thread == 1){
    chrDfList <- lapply(1:length(mzIter), function(i) {
      utils::setTxtProgressBar(pb, i)
      x <- mzIter[[i]]
      loop(x)
    })
  }else if(thread > 1){
    chrDfList <- BiocParallel::bplapply(1:length(mzIter), function(i) {
      x <- mzIter[[i]]
      loop(x)
    }, BPPARAM = BiocParallel::SnowParam(workers = thread,
                                         progressbar = TRUE))
  }else stop("Thread is wrong!")
  return(chrDfList)
}
chrDfList_bins <- generateBin(ndata = ndata, thread = 4)
length(chrDfList_bins) # 22699
logicalVec <- sapply(chrDfList_bins, function(x) {
  if(nrow(x) == 0) return(FALSE)
  else if(is.null(x)) return(FALSE)
  else return(TRUE)
})
chrDfList_bins <- chrDfList_bins[logicalVec]
length(chrDfList_bins) # 17272
saveRDS(chrDfList_bins, file = "D:/fudan/Projects/2024/MetaboProcess/Progress/generate_bins/240627/chrDfList_bins.rds")
scanNumbers <- sapply(chrDfList_bins, nrow)
idx <- which(scanNumbers >= 100)
chrDfList_bins <- chrDfList_bins[idx]
length(chrDfList_bins) # 1491
chrDfList_bins <- lapply(chrDfList_bins, function(x) {
  x$intensity <- smoothMean(x$intensity)
  return(x)
})
i <- 18
plot_chrDf(chrDfList_bins[[i]], noise = noiseEstimation(chrDfList_bins[[i]]))
ZOI <- FindZOI(chrDfList_bins[[i]], noise = noiseEstimation(chrDfList_bins[[i]]),preNum = 3)
plot_chrDf(ZOI[[1]])
ZOI <- baselineEs(chrDf = ZOI[[1]], loops = 8, tol_m = 30, threshold = 1)
plot_chrDf(ZOI, baseline = TRUE)
# deductBaseline <- function(x){
#   x$intensity <- x$intensity - x$baseline
#   return(x)
# }
# ZOI <- deductBaseline(ZOI)
ZOI <- edgeTrack(chrDf = ZOI)
plot_chrDf(ZOI, baseline = TRUE)
ZOI <- edgeTrack_finer(ZOI, tol_m = 30)
ZOI <- edgeTrack_crude(ZOI, tol_m = 30)
plot_chrDf(ZOI, baseline = TRUE)
chrDf_ZOIs_2 <- lapply(1:length(chrDf_ZOIs), function(i) {
  x <- chrDf_ZOIs[[i]]
  chrDf_tmp <- tryCatch({
    edgeTrack_finer(x, tol_m = tol_m2)
  }, error = function(e){
    edgeTrack_crude(x, tol_m = tol_m2)
  })
})
