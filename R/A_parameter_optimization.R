#' @title generateBin
#' @description
#' Divide a sample into multiple chromatograms according to the length of the bin.
#'
#' @param ndata nth sample data, it is a MsExperiment object.
#' @param bin bin mz.
#' @param mslevel mslevel, 1 or 2.
#' @param thread Number of threads in parallel.
#'
#' @return A list consisting of chrDf.
#' @export
#'
#' @examples
#' chrDf_bins <- generateBin(ndata = ndata, bin = 0.05, mslevel = 1, thread = 24)
generateBin <- function(ndata, bin = 0.05, mslevel = 1, thread = 1){
  sps <- xcms::spectra(ndata) %>% Spectra::filterMsLevel(mslevel)
  fileh <- mzR::openMSfile(MsExperiment::sampleData(ndata)$sample_path, backend = NULL)
  runInfo <- mzR::runInfo(fileh)
  mzR::close(fileh);rm(fileh)
  lowMz <- runInfo$lowMz;highMz <- runInfo$highMz
  rtSps <- rtSps <- Spectra::rtime(sps)
  peaksData <- Spectra::peaksData(sps)
  mzRange <-seq(lowMz, highMz, by = bin)
  mzIter <- lapply(1:(length(mzRange) - 1), function(i) {c(mzRange[i], mzRange[i + 1])})
  loop <- function(x){
    mz_range <- x
    chrDf <- purrr::list_rbind(lapply(1:length(peaksData), function(j) {
      peakMat <- peaksData[[j]]
      idx <- which(peakMat[, "mz"] >= mz_range[1] & peakMat[, "mz"] <= mz_range[2])
      if(length(idx) == 0) return(NULL)
      mz <- mean(peakMat[idx,"mz"])
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

#' @title noiseEstimation
#' @description
#' Estimates the noise of a chromatogram, which can be used to find the apex of a chromatographic peak.
#'
#' @param chrDf A chrDf.
#'
#' @return A numeric.
#' @export
#'
#' @examples
#' noiseEstimation(chrDf_bins[[30]])
noiseEstimation <- function(chrDf){
  intensity <- sort(chrDf$intensity[chrDf$intensity > 0])
  intLength <- length(intensity)
  if(intLength == 1) return(intensity)
  for(i in (intLength - 1):1){
    int <- intensity[i + 1]
    int_vec <- intensity[i:1]
    noiEsti <- mean(int_vec) + 3 * sd(int_vec)
    if(int < noiEsti) break
  }
  noise <- intensity[i + 1]
  return(noise)
}
#' @title FindZOI
#' @description
#' Find ZOIs from a chrDf of the bin.
#'
#' @param chrDf A chrDf.
#' @param noise noise.
#' @param preNum preNum.
#' @param etlD etlD. Extend width.
#' @param IETH IE threshold.
#'
#' @return A chrDfList.
#' @export
#'
#' @examples
#' ZOIList <- FindZOI(chrDf = chrDf, noise = noiseEstimation(chrDf), etlD = 1, IETH = 10)
FindZOI <- function(chrDf, noise = NA, preNum = 3, etlD = 1,  IETH = 1.75){
  aboveTHidx <- which(chrDf$intensity >= noise)
  if(length(aboveTHidx) == 0) return(NULL)
  candidateSegInd <- split(aboveTHidx, cumsum(c(1, diff(aboveTHidx) != 1)))
  chrDfList_ZOI_tmp <- lapply(1:length(candidateSegInd), function(i) {
    idx <- candidateSegInd[[i]]
    int <- chrDf$intensity[idx]
    intLength <- length(int)
    if(intLength > preNum){
      extend <- round(intLength / etlD)
      idx1 <- idx[1] - extend
      if(idx1 < 0) idx1 <- 0
      idx2 <- idx[intLength] + extend
      if(idx2 > length(chrDf$intensity)) idx2 <- length(chrDf$intensity)
      idx <- c(idx1:(idx[1] - 1), idx, (idx[intLength] + 1):idx2)
      chrDf_tmp <- chrDf[idx, ]
      return(chrDf_tmp)
    }else return(NULL)
  })
  chrDfList_ZOI_tmp <- chrDfList_ZOI_tmp[!sapply(chrDfList_ZOI_tmp, is.null)]
  IE_vec <- sapply(1:length(chrDfList_ZOI_tmp), function(i) {
    cal_IE(chrDfList_ZOI_tmp[[i]]$intensity)
  })
  chrDfList_ZOI_tmp <- chrDfList_ZOI_tmp[which(IE_vec < IETH)]
  if(length(chrDfList_ZOI_tmp) == 0) chrDfList_ZOI_tmp <- NULL
  return(chrDfList_ZOI_tmp)
}
#' @title edgeTrack
#' @description
#' locate the edge of a peak.
#'
#' @param chrDf A chrDf with baseline.
#'
#' @return A chrDf but start and end are edges.
#' @export
#'
#' @examples
#' chrDf_j <- ZOIList[[j]]
#' chrDf_j <- baselineEs(chrDf_j, loops = 8, tol_m = 30)
#' chrDf_j_new <- edgeTrack(chrDf = chrDf_j)
edgeTrack <- function(chrDf){
  inflect <- function(x, threshold = 1){
    up   <- sapply(1:threshold, function(n) c(x[-(seq(n))], rep(NA, n)))
    down <-  sapply(-1:-threshold, function(n) c(rep(NA,abs(n)), x[-seq(length(x), length(x) - abs(n) + 1)]))
    a    <- cbind(x,up,down)
    list(minima = which(apply(a, 1, min) == a[,1]), maxima = which(apply(a, 1, max) == a[,1]))
  }
  #browser()
  deintensity <- chrDf$intensity - chrDf$baseline
  apex_idx <- which.max(deintensity)
  aboveTHidx <- which(deintensity > 0)
  if(length(aboveTHidx) == 0) return(NULL)
  candidateSegInd <- split(aboveTHidx, cumsum(c(1, diff(aboveTHidx) != 1)))
  logical_Seg <- sapply(candidateSegInd, function(x) {
    if(apex_idx %in% x) return(TRUE)
    else return(FALSE)
  })
  SegInd <- candidateSegInd[[which(logical_Seg)]]
  if((SegInd[1] - 1) >= 1) SegInd <- c(SegInd[1] - 1, SegInd)
  if((SegInd[length(SegInd)] + 1) <= nrow(chrDf)) SegInd <- c(SegInd, SegInd[length(SegInd)] + 1)
  #SegInd <- c(SegInd[1] - 1, SegInd, SegInd[length(SegInd)] + 1)
  chrDf <- chrDf[SegInd, ]
  deintensity <- chrDf$intensity - chrDf$baseline
  #tmp <- MassSpecWavelet::peakDetectionCWT(chrDf$intensity, scales = seq(1, 10, 2))
  #xcms::peaksWithCentWave(int = chrDf$intensity, rt = chrDf$rt, peakwidth = c(20, 30))
  tmp <- tryCatch({
    MassSpecWavelet::peakDetectionCWT(chrDf$intensity, scales = seq(1, 8, 2))
  },
  error = function(e){
    NULL
  })
  if(is.null(tmp)) return(chrDf)
  peakIdx <- tmp$majorPeakInfo$allPeakIndex
  apex_idx <- which.max(chrDf$intensity)
  apex_int <- chrDf$intensity[apex_idx]
  peakIdx <- peakIdx[-which.min(abs(peakIdx - apex_idx))]
  if(length(peakIdx) == 1){
    int_tmp <- chrDf$intensity
    if(apex_idx < peakIdx) int_tmp[-(apex_idx:peakIdx)] <- apex_int
    else int_tmp[-(peakIdx:apex_idx)] <- apex_int
    bottom_idx <- which.min(int_tmp)
  }else if(length(peakIdx) > 1){
    tmp <- peakIdx - apex_idx
    if(all(tmp > 0) | all(tmp < 0)){
      peakIdx <- peakIdx[which.min(abs(peakIdx - apex_idx))]
      int_tmp <- chrDf$intensity
      if(apex_idx < peakIdx) int_tmp[-(apex_idx:peakIdx)] <- apex_int
      else int_tmp[-(peakIdx:apex_idx)] <- apex_int
      bottom_idx <- which.min(int_tmp)
    }else{
      peakIdx_a <- peakIdx[tmp < 0]
      peakIdx_a <- peakIdx_a[length(peakIdx_a)]
      peakIdx_b <- peakIdx[tmp > 0]
      peakIdx_b <- peakIdx_b[1]
      int_tmp_a <- chrDf$intensity
      int_tmp_a[-(peakIdx_a:apex_idx)] <- apex_int
      int_tmp_b <- chrDf$intensity
      int_tmp_b[-(apex_idx:peakIdx_b)] <- apex_int
      bottom_idx <- c(which.min(int_tmp_a), which.min(int_tmp_b))
    }
  }else if(length(peakIdx) == 0) bottom_idx <- c(1, nrow(chrDf))
  if(length(bottom_idx) == 1){
    if(bottom_idx > apex_idx) chrDf <- chrDf[1:bottom_idx, ]
    else chrDf <- chrDf[bottom_idx:nrow(chrDf), ]
  }else if(length(bottom_idx) == 2){
    chrDf <- chrDf[bottom_idx[1]:bottom_idx[2], ]
  }else stop("Wrong!")
  return(chrDf)
}
