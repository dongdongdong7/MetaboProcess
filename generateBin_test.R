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
length(chrDfList_bins) # 56748
idx1 <- which(scanNumbers > 0 & scanNumbers < 10 & !is.na(scanNumbers))
idx2 <- which(scanNumbers > 10 & scanNumbers < 50 & !is.na(scanNumbers))
idx3 <- which(scanNumbers > 100 & scanNumbers < 200 & !is.na(scanNumbers))
idx4 <- which(scanNumbers > 200 & !is.na(scanNumbers))
i <- 904
plot_chrDf(chrDfList_bins[[idx3[i]]], noise = noiseEstimation(chrDfList_bins[[idx3[i]]]))
scanNumbers <- sapply(chrDfList_bins, nrow)
plot(x = 1:length(scanNumbers[scanNumbers > 0]), y = scanNumbers[scanNumbers > 0])
