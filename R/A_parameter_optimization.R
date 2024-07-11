#' @title generateBin
#' @description
#' Divide a sample into multiple chromatograms according to the length of the bin.
#'
#' @param ndata nth sample data, it is a MsExperiment object.
#' @param bin bin mz.
#' @param slide slide size.
#' @param mslevel mslevel.
#' @param thread thread.
#'
#' @return A chrDf list.
#' @export
#'
#' @examples
#' chrDfList_bins <- generateBin(ndata = ndata, thread = 4)
generateBin <- function(ndata, bin = 0.05, slide = 0.05, mslevel = 1, thread = 1){
  start_time <- Sys.time()
  sps <- xcms::spectra(ndata) %>% Spectra::filterMsLevel(mslevel)
  lowMz <- round(min(sapply(Spectra::mz(sps), min)))
  highMz <- round(max(sapply(Spectra::mz(sps), max)))
  mzData <- Spectra::mz(sps)
  intData <- Spectra::intensity(sps)
  rtime <- Spectra::rtime(sps)
  mzS <- seq(lowMz, highMz - bin, by = slide)
  mzE <- seq(lowMz + bin, highMz, by = slide)
  maxN <- max(c(length(mzS), length(mzE)))
  mzIter <- lapply(1:maxN, function(i) {c(mzS[i], mzE[i])})
  loop <- function(i){
    currmzRange <- mzIter[[i]]
    chrDf <- purrr::list_rbind(lapply(1:length(mzData), function(j) {
      index <- which(mzData[[j]] >= currmzRange[1] & mzData[[j]] < currmzRange[2])
      if(length(index) == 0){mz <- NA;intensity <- 0;rt <- rtime[j]}
      else{mz <- mzData[[j]][index][which.max(intData[[j]][index])];intensity <- max(intData[[j]][index]);rt <- rtime[j]}
      return(dplyr::tibble(mz = mz, intensity = intensity, rt = rt))
    }))
    return(chrDf)
  }
  pb <- utils::txtProgressBar(max = length(mzIter), style = 3)
  if(thread == 1){
    chrDfList <- lapply(1:length(mzIter), function(i) {
      utils::setTxtProgressBar(pb, i)
      loop(i)
    })
  }else if(thread > 1){
    chrDfList <- BiocParallel::bplapply(1:length(mzIter), function(i) {
      loop(i)
    }, BPPARAM = BiocParallel::SnowParam(workers = thread,
                                         progressbar = TRUE))
  }else stop("Thread is wrong!")
  end_time <- Sys.time()
  print(end_time - start_time)
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
  #browser()
  intensity <- sort(chrDf$intensity[chrDf$intensity > 0])
  intLength <- length(intensity)
  if(intLength == 1) return(intensity)
  for(i in (intLength - 1):1){
    if(i == 1) break
    int <- intensity[i + 1]
    int_vec <- intensity[i:1]
    noiEsti <- mean(int_vec) + 3 * sd(int_vec)
    if(int < noiEsti) break
  }
  noise <- intensity[i + 1]
  return(noise)
}
#' @title pickZOI
#' @description
#' Direct peak(ZOI) finding on one chromatogram.
#'
#' @param chrDf A chrDf.
#' @param smoothPara smoothPara.
#' @param baselinePara baselinePara.
#' @param sn For filtering peaks with low signal-to-noise ratio.
#' @param preNum A peak needs to have preNum of data points greater than noise0.
#' @param tol_m Angle tolerance for baseline calculation and edge finding.
#' @param snthresh
#' The signal-to-noise ratio used for the MatchedFilter,
#' which hides the peak when the signal-to-noise ratio is lower than this value.
#'
#' @return A ZOI list contains chrDf.
#' @export
#'
#' @examples
#' data("chrDf_mrm_test")
#' plot_chrDf(chrDf_mrm_test)
#' ZOIList <- pickZOI(chrDf = chrDf_mrm_test, sn = 3)
#' length(ZOIList) # 6
#' plot_chrDf(ZOIList[[1]], baseline = TRUE)
pickZOI <- function(chrDf, smoothPara = get_smoothPara(), baselinePara = get_baselinePara(),
                    sn = 3, preNum = 3, tol_m = 10, snthresh = 0.5){ # chrDf 通常需要是一个bin窗口
  #Sample number
  sampleIdx <- attributes(chrDf)$sample
  #Smooth
  if(smoothPara$smooth){
    if(smoothPara$method == "mean") chrDf$intensity <- smoothMean(chrDf$intensity, size = smoothPara$size)
    else if(smoothPara$method == "sg") chrDf$intensity <- smoothSg(chrDf$intensity, p = smoothPara$p, n = smoothPara$n, m = smoothPara$m, ts = smoothPara$ts)
  }
  # Nise0
  noise0 <- noiseEstimation(chrDf)
  #attributes(chrDf)$noise0 <- noise0
  chrDf <- baselineEs(chrDf = chrDf, threshold = baselinePara$threshold, tol_m = baselinePara$tol_m, loops = baselinePara$loops)
  aboveTHidx <- which(chrDf$intensity > chrDf$baseline)
  if(length(aboveTHidx) == 0) return(NULL)
  candidateSegInd <- split(aboveTHidx, cumsum(c(1, diff(aboveTHidx) != 1)))
  candidateSegInd <- candidateSegInd[sapply(candidateSegInd, function(i) {
    if(length(which(chrDf$intensity[i] > noise0)) > preNum) return(TRUE)
    else return(FALSE)
  })]
  if(length(candidateSegInd) == 0) return(NULL)
  #browser()
  #deltaTime <- mean(diff(chrDf$rt))
  ZOIList <- lapply(1:length(candidateSegInd), function(i) {
    #print(i)
    idx <- candidateSegInd[[i]]
    idxLength <- length(idx)
    idx1 <- idx[1] - 1
    if(idx1 < 1) idx1 <- 1
    idx2 <- idx[idxLength] + 1
    if(idx2 > length(chrDf$intensity)) idx2 <- length(chrDf$intensity)
    idx <- unique(c(idx1:idx[1], idx, idx[idxLength]:idx2))
    ZOI <- chrDf[idx, ]
    ZOIWidth <- (max(ZOI$rt) - min(ZOI$rt)) / 2 # half
    fwhm <- round(ZOIWidth / 2)
    tmp <- tryCatch({
      xcms::peaksWithMatchedFilter(ZOI$intensity, ZOI$rt, fwhm = fwhm, snthresh = snthresh)
    },
    error = function(e){
      matrix(NA, nrow = 0, ncol = 8)
    })
    if(nrow(tmp) > 1){ # Multi peaks
      top_idx <- sort(sapply(1:nrow(tmp), function(j) {
        which(dplyr::near(ZOI$rt, tmp[j, "rt"], tol = 0.005))
      }))
      bottom_idx <- sapply(1:(length(top_idx) - 1), function(j) {
        a_idx <- top_idx[j];b_idx <- top_idx[j + 1]
        int <- ZOI$intensity
        int[1:a_idx] <- NA;int[b_idx:length(int)] <- NA
        which.min(int)
      })
      bottom_idx <- c(1, bottom_idx, length(ZOI$intensity))
      multiZOIList <- lapply(top_idx, function(t) {
        if(t %in% bottom_idx) return(NULL)
        a_idx <- bottom_idx[bottom_idx < t];a_idx <- a_idx[length(a_idx)]
        b_idx <- bottom_idx[bottom_idx > t];b_idx <- b_idx[1]
        ZOI_t <- ZOI[a_idx:b_idx, ]
        ZOI_t <- edgeTrack_crude(ZOI_t, preNum = preNum, tol_m = tol_m)
      })
      return(multiZOIList)
    }else{ # Single peak
      ZOI <- edgeTrack_crude(ZOI, preNum = preNum, tol_m = tol_m)
    }
  })
  #browser()
  ZOIList <- purrr::list_flatten(ZOIList)
  ZOIList <- ZOIList[sapply(ZOIList, function(x) {
    if(is.null(x)) return(FALSE)
    sn_tmp <- max(x$intensity) / max(x$baseline)
    if(sn_tmp > sn) return(TRUE)
    else return(FALSE)
  })]
  if(length(ZOIList) == 0) return(NULL)
  for(i in 1:length(ZOIList)){
    attributes(ZOIList[[i]])$sample <- sampleIdx
  }
  return(ZOIList)
}

# tightenZOI(chrDf_test)
tightenZOI <- function(chrDf, dight = 3){
  #browser()
  mzVec <- round(chrDf$mz, digits = dight)
  pointNum <- nrow(chrDf)
  topIdx <- which.max(chrDf$intensity)
  startIdx <- topIdx - 2;endIdx <- topIdx + 2
  if(startIdx < 1) startIdx <- 1
  if(endIdx > nrow(chrDf)) endIdx <- nrow(chrDf)
  tmp <- boxplot.stats(mzVec[startIdx:endIdx])$stats
  LB <- tmp[1];RB <- tmp[5]
  if(all(mzVec[startIdx:endIdx] > LB | abs(mzVec[startIdx:endIdx] - LB)<exp(-30)) & all(mzVec[startIdx:endIdx] < RB | abs(mzVec[startIdx:endIdx] - RB)<exp(-30))){
    startIdx <- startIdx - 1
    while(startIdx >= 1){
      tmp <- boxplot.stats(mzVec[startIdx:endIdx])$stats
      LB <- tmp[1];RB <- tmp[5]
      if((mzVec[startIdx] > LB | abs(mzVec[startIdx] - LB) < exp(-30)) & (mzVec[startIdx] < RB | abs(mzVec[startIdx] - RB) < exp(-30))){
        startIdx <- startIdx - 1
      }else{
        startIdx <- startIdx + 1
        break
      }
    }
    if(startIdx < 1) startIdx <- 1
    endIdx <- endIdx + 1
    while(endIdx <= pointNum){
      tmp <- boxplot.stats(mzVec[startIdx:endIdx])$stats
      LB <- tmp[1];RB <- tmp[5]
      if((mzVec[endIdx] > LB | abs(mzVec[endIdx] - LB) < exp(-30)) & (mzVec[endIdx] < RB | abs(mzVec[endIdx] - RB) < exp(-30))){
        endIdx <- endIdx + 1
      }else{
        endIdx <- endIdx - 1
        break
      }
    }
    if(endIdx > pointNum) endIdx <- pointNum
  }else{
    chrDf <- NULL
  }
  chrDf <- chrDf[startIdx:endIdx, ]
  return(chrDf)
}
checkShapeZOI <- function(chrDf, fwhm = NA,snthresh = 1){
  #browser()
  if(is.na(fwhm)) fwhm <- round((max(chrDf$rt) - min(chrDf$rt)))
  tmp <- xcms::peaksWithMatchedFilter(int = chrDf$intensity, rt = chrDf$rt, fwhm = fwhm, snthresh = snthresh)
  if(nrow(tmp) == 0) return(NULL)
  else{
    # if(nrow(tmp) > 1){
    #   tmp <- tmp[which.max(tmp[, "maxo"]), ]
    # }
    # rtmin <- tmp[, "rtmin"];rtmax <- tmp[, "rtmax"]
    # start_idx <- which.min(abs(chrDf$rt - rtmin))
    # end_idx <- which.min(abs(chrDf$rt - rtmax))
    # chrDf <- chrDf[start_idx:end_idx, ]
    return(chrDf)
  }
}

#' @title ZOIList2ZOITable
#' @description
#' Make ZZOIList to ZOITable.
#'
#' @param ZOIList ZOIList
#'
#' @return ZOITable.
#'
#' @examples ZOITable.
#' ZOIList <- readRDS("D:/fudan/Projects/2024/MetaboProcess/Progress/generate_bins/240705/ZOIList.rds")
#' ZOITable <- ZOIList2ZOITable(ZOIList)
ZOIList2ZOITable <- function(ZOIList){
  #browser()
  pb <- utils::txtProgressBar(max = length(ZOIList), style = 3)
  tbList <- lapply(1:length(ZOIList), function(i) {
    utils::setTxtProgressBar(pb, i)
    x <- ZOIList[[i]]
    topIdx <- which.max(x$intensity)
    mz <- x$mz[topIdx]
    mzmin <- min(x$mz)
    mzmax <- max(x$mz)
    rt <- x$rt[topIdx]
    rtmin <- min(x$rt)
    rtmax <- max(x$rt)
    deltaRt <- (rtmax - rtmin) / nrow(x)
    into <- sum(x$intensity * deltaRt)
    intb <- sum((x$intensity - x$baseline) * deltaRt)
    maxo <- x$intensity[topIdx]
    sn <- maxo / max(x$baseline)
    sample <- attributes(x)$sample
    tb <- dplyr::tibble(cpid = paste0("CP", i),mz = mz, mzmin = mzmin, mzmax = mzmax, rt = rt, rtmin = rtmin, rtmax = rtmax, into = into, intb = intb, maxo = maxo, sn = sn, sample = sample, ms_level = 1, is_filled = FALSE)
    return(tb)
  })
  ZOITable <- purrr::list_rbind(tbList)
  return(ZOITable)
}

#' @title FindIso
#' @description
#' Find isotopes from ZOI.
#'
#' @param ZOITable ZOITable.
#' @param data_QC data_QC.
#' @param thread thread.
#'
#' @return A ISOTable.
#' @export
#'
#' @examples
#' isoTable <- FindIso(ZOITable, data_QC, thread = 3)
FindIso <- function(ZOITable, data_QC, thread = 3){
  #browser()
  pb <- utils::txtProgressBar(max = length(data_QC), style = 3)
  if(thread == 1){
    isoList <- lapply(1:length(data_QC), function(i) {
      utils::setTxtProgressBar(pb, i)
      data <- data_QC[i]
      pd <- MsExperiment::sampleData(data)
      set <- xcms::xcmsSet(files = pd$sample_path, phenoData = as.data.frame(pd), mslevel = 1L, BPPARAM = BiocParallel::SerialParam())
      chromPeaks_new <- as.matrix(ZOITable[which(ZOITable$sample == i), 2:12])
      rownames(chromPeaks_new) <- ZOITable[which(ZOITable$sample == i), ]$cpid
      set@peaks <- chromPeaks_new
      set.seed(2)
      ex.cliqueGroups <- cliqueMS::getCliques(set, filter = TRUE)
      ex.Isotopes <- cliqueMS::getIsotopes(ex.cliqueGroups, ppm = 10)
      isoTable <- ex.Isotopes@peaklist[stringr::str_detect(ex.Isotopes@peaklist$isotope, "\\[*\\]"), ]
      return(isoTable)
    })
  }else if(thread > 1){
    cl <- snow::makeCluster(thread)
    doSNOW::registerDoSNOW(cl)
    opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
                                                                 n))
    isoList <- foreach::`%dopar%`(foreach::foreach(i = 1:length(data_QC),
                                                   .packages = c("MsExperiment", "cliqueMS", "xcms", "stringr", "BiocParallel"),
                                                   .options.snow = opts), {
                                                     data <- data_QC[i]
                                                     pd <- MsExperiment::sampleData(data)
                                                     set <- xcms::xcmsSet(files = pd$sample_path, phenoData = as.data.frame(pd), mslevel = 1L, BPPARAM = BiocParallel::SerialParam())
                                                     chromPeaks_new <- as.matrix(ZOITable[which(ZOITable$sample == i), 2:12])
                                                     rownames(chromPeaks_new) <- ZOITable[which(ZOITable$sample == i), ]$cpid
                                                     set@peaks <- chromPeaks_new
                                                     set.seed(2)
                                                     ex.cliqueGroups <- cliqueMS::getCliques(set, filter = TRUE)
                                                     ex.Isotopes <- cliqueMS::getIsotopes(ex.cliqueGroups, ppm = 10)
                                                     isoTable <- ex.Isotopes@peaklist[stringr::str_detect(ex.Isotopes@peaklist$isotope, "\\[*\\]"), ]
                                                     return(isoTable)
                                                   })
    snow::stopCluster(cl)
    gc()
  }
  isoTable <- purrr::list_rbind(isoList)
  return(isoTable)
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
FindZOI <- function(chrDf, noise = NA, preNum = 3, etlD = 1,  IETH = 1.75, sn = 3){
  #browser()
  aboveTHidx <- which(chrDf$intensity >= noise)
  if(length(aboveTHidx) == 0) return(NULL)
  candidateSegInd <- split(aboveTHidx, cumsum(c(1, diff(aboveTHidx) != 1)))
  chrDfList_ZOI_tmp <- lapply(1:length(candidateSegInd), function(i) {
    idx <- candidateSegInd[[i]]
    int <- chrDf$intensity[idx]
    sn_tmp <- max(int) / noise
    if(sn_tmp < sn) return(NULL)
    intLength <- length(int)
    if(intLength > preNum){
      extend <- round(intLength / etlD)
      idx1 <- idx[1] - extend
      if(idx1 < 0) idx1 <- 1
      idx2 <- idx[intLength] + extend
      if(idx2 > length(chrDf$intensity)) idx2 <- length(chrDf$intensity)
      idx <- unique(c(idx1:idx[1], idx, idx[intLength]:idx2))
      chrDf_tmp <- chrDf[idx, ]
      return(chrDf_tmp)
    }else return(NULL)
  })
  chrDfList_ZOI_tmp <- chrDfList_ZOI_tmp[!sapply(chrDfList_ZOI_tmp, is.null)]
  if(length(chrDfList_ZOI_tmp) == 0) return(NULL)
  IE_vec <- sapply(1:length(chrDfList_ZOI_tmp), function(i) {
    cal_IE(chrDfList_ZOI_tmp[[i]]$intensity)
  })
  chrDfList_ZOI_tmp <- chrDfList_ZOI_tmp[which(IE_vec < IETH)]
  if(length(chrDfList_ZOI_tmp) == 0) return(NULL)
  return(chrDfList_ZOI_tmp)
}
#' @title FindZOI_baseline
#' @description
#' Find ZOIs from a big ZOI by baseline.
#'
#' @param chrDf chrDf
#' @param noise noise
#' @param sn sn
#' @param preNum preNum
#'
#' @return A chrDf list
#' @export
#'
#' @examples
#' FindZOI_baseline(chrDf, noise = 1000)
FindZOI_baseline <- function(chrDf, noise, sn = 0,preNum = 3){
  #browser()
  aboveTHidx <- which(chrDf$intensity > chrDf$baseline)
  if(length(aboveTHidx) == 0) return(NULL)
  candidateSegInd <- split(aboveTHidx, cumsum(c(1, diff(aboveTHidx) != 1)))
  candidateSegInd <- candidateSegInd[which(sapply(candidateSegInd, length) > preNum)]
  if(length(candidateSegInd) == 0) return(NULL)
  chrDfList_ZOI_tmp <- lapply(1:length(candidateSegInd), function(i) {
    idx <- candidateSegInd[[i]]
    int <- chrDf$intensity[idx]
    if(max(int) < noise) return(NULL)
    sn_tmp <- max(int) / noise
    if(sn_tmp < sn) return(NULL)
    intLength <- length(int)
    if(intLength > preNum){
      idx1 <- idx[1] - 1
      if(idx1 < 0) idx1 <- 1
      idx2 <- idx[intLength] + 1
      if(idx2 > length(chrDf$intensity)) idx2 <- length(chrDf$intensity)
      idx <- unique(c(idx1:idx[1], idx, idx[intLength]:idx2))
      chrDf_tmp <- chrDf[idx, ]
      return(chrDf_tmp)
    }else return(NULL)
  })
  chrDfList_ZOI_tmp <- chrDfList_ZOI_tmp[!sapply(chrDfList_ZOI_tmp, is.null)]
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
  deintensity <- round(chrDf$intensity - chrDf$baseline, 4)
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
  if(is.null(tmp)) return(chrDf) # return(NULL)
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
#' @title bins2ZOIs
#' @description
#' Bins 's chrDf to ZOIs's chrDf.
#'
#' @param chrDfList A chrDfList, one chrDf means a bin.
#' @param smooth smooth type.
#' @param size smoothMean size.
#' @param p smoothSg p.
#'
#' @return A chrDfList.
#' @export
#'
#' @examples
#' bins2ZOIs(chrDfList = chrDf_bins, thread = 4)
bins2ZOIs <- function(chrDfList,
                      smooth = "mean", size = 3, p = 3,
                      etlD = 1, IETH = 1.75, preNum = 3, sn = 3,
                      loops = 8, tol_m = 30, threshold = 1,
                      thread = 1){
  loop <- function(i){
    chrDf <- chrDfList[[i]]
    if(nrow(chrDf) == 0) return(NULL)
    if(smooth == "mean") chrDf$intensity <- smoothMean(chrDf$intensity, size = size)
    else if(smooth == "sg") chrDf$intensity <- smoothSg(chrDf$intensity, p = p)
    else stop("wrong for smooth")
    chrDf_noise <- noiseEstimation(chrDf)
    plot_chrDf(chrDf, noise = chrDf_noise)
    ZOIList <- FindZOI(chrDf = chrDf, noise = chrDf_noise, etlD = etlD, IETH = IETH, preNum = preNum,sn = sn)
    if(is.null(ZOIList)) return(NULL)
    ZOIList <- lapply(1:length(ZOIList), function(j) {
      chrDf_j <- ZOIList[[j]]
      #cal_IE(chrDf_j$intensity)
      chrDf_j <- baselineEs(chrDf = chrDf_j, loops = loops, tol_m = tol_m, threshold = threshold)
      #plot_chrDf(chrDf_j, noise = chrDf_noise, baseline = TRUE)
      chrDf_j_new <- edgeTrack(chrDf = chrDf_j)
      #plot_chrDf(chrDf_j_new, noise = chrDf_noise, baseline = TRUE)
      # chrDf_j_new <- tryCatch({
      #   edgeTrack_finer(chrDf_j, tol_m = tol_m)
      # },error = function(e){
      #   edgeTrack_crude(chrDf_j, tol_m = tol_m)
      # })
      return(chrDf_j_new)
    })
    return(ZOIList)
    #plot_chrDf(chrDf, noise = chrDf_noise)
  }
  pb <- utils::txtProgressBar(max = length(chrDfList), style = 3)
  #browser()
  if(thread == 1){
    t <- lapply(1:length(chrDfList), function(i) {
      utils::setTxtProgressBar(pb, i)
      loop(i)
    })
  }else if(thread > 1){
    cl <- snow::makeCluster(thread)
    doSNOW::registerDoSNOW(cl)
    opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
                                                                 n))
    envir <- environment(loop)
    parallel::clusterExport(cl, varlist = ls(envir), envir = envir)
    envir <- environment(FindZOI)
    parallel::clusterExport(cl, varlist = ls(envir), envir = envir)
    t <- foreach::`%dopar%`(foreach::foreach(i = 1:length(chrDfList),
                                             #.packages = c("dplyr"),
                                             .export = c("smoothMean", "smoothSg", "noiseEstimation","FindZOI","baselineEs","edgeTrack"),
                                             .options.snow = opts),
                            {
                              loop(i)
                            })
    snow::stopCluster(cl)
    gc()
  }
  t <- purrr::list_flatten(t)
  t <- t[!sapply(t, is.null)]
  return(t)
}
#edgeTrack_finer(chrDf_ZOIs[[274]])
edgeTrack_finer <- function(chrDf, preNum  = 3, tol_m = 30){
  getLine <- function(A, B){
    slope <- (B[2] - A[2]) / (B[1] - A[1])
    intercept <- A[2] - slope * A[1]
    return(c(slope = slope, intercept = intercept))
  }
  #browser()
  #chrDf
  #plot_chrDf(chrDf, baseline = TRUE)
  apex_idx <- which.max(chrDf$intensity)
  int <- chrDf$intensity;rt <- chrDf$rt
  int_a <- int[1:apex_idx];int_b <- int[apex_idx:length(int)]
  while(!all(diff(int_a) >0) | !all(diff(int_b) < 0)){
    int <- smoothMean(int, size = 3)
    #plot(x = 1:length(int), y = int)
    apex_idx <- which.max(int)
    int_a <- int[1:apex_idx];int_b <- int[apex_idx:length(int)]
  }
  int_norm <- (int - min(int)) / (max(int) - min(int))
  rt_norm <- (rt - min(rt)) / (max(rt) - min(rt))
  # left
  a_start <- apex_idx - ((preNum - 1) / 2) - 1
  if(a_start <= 1) a_start <- 1
  b_start <- apex_idx + ((preNum - 1) / 2) + 1
  if(b_start >= nrow(chrDf)) b_start <- nrow(chrDf)
  a <- a_start;b <- b_start
  #browser()
  while(a!=1 | b!=nrow(chrDf)){
    a0 <- a - 1;b0 <- b + 1
    if(a == 1) a0 <- a  + 1
    if(b == nrow(chrDf)) b0 <- b - 1
    A <- c(rt_norm[a], int_norm[a]);A0 <- c(rt_norm[a0], int_norm[a0])
    B <- c(rt_norm[b], int_norm[b]);B0 <- c(rt_norm[b0], int_norm[b0])
    tmpA <- getLine(A, A0);slopeA <- tmpA["slope"];interceptA <- tmpA["intercept"]
    tmpB <- getLine(B, B0);slopeB <- tmpB["slope"];interceptB <- tmpB["intercept"]
    Cx <- (interceptB - interceptA) / (slopeA - slopeB);Cy <- slopeA * Cx + interceptA;C <- c(Cx, Cy)
    # print(a)
    # print(b)
    if(C[1] > B[1] | C[1] < A[1]){
      angleA0 <- atan(abs(slopeA)) * (180 / pi)
      angleB0 <- atan(abs(slopeB)) * (180 / pi)
      if(angleA0 < tol_m & angleB0 < tol_m){break}
      else if(angleA0 < angleB0){b <- b + 1;next}
      else if(angleA0 >= angleB0){a <- a - 1;next}
      # if(A[2] < B[2]){b <- b + 1;next}
      # else if(A[2] >= B[2]){a <- a - 1;next}
      # if(C[1] > B[1]){b <- b + 1;next}
      # else if(C[1] < A[1]){a <- a - 1;next}
    }
    # plot(x = rt_norm, y = int_norm)
    # points(x = A[1], y = A[2], col = "red")
    # points(x = B[1], y = B[2], col = "red")
    # points(x = A0[1], y = A0[2], col = "blue")
    # points(x = B0[1], y = B0[2], col = "blue")
    # points(x = C[1], y = C[2], col = "green")
    Alength <- sqrt((B[1] - C[1])^2 + (B[2] - C[2])^2)
    Blength <- sqrt((A[1] - C[1])^2 + (A[2] - C[2])^2)
    Clength <- sqrt((B[1] - A[1])^2 + (B[2] - A[2])^2)
    COSA <- (Blength^2 + Clength^2 - Alength^2) / (2*Blength*Clength)
    COSB <- (Alength^2 + Clength^2 - Blength^2) / (2*Alength*Clength)
    COSC <- (Blength^2 + Alength^2 - Clength^2) / (2*Blength*Alength)
    angleA <- acos(COSA) * (180 / pi)
    angleB <- acos(COSB) * (180 / pi)
    angleC <- acos(COSC) * (180 / pi)
    #print(a);print(b)
    # if(a == 1) browser()
    if(angleA > tol_m & angleB > tol_m){
      if(angleA > angleB) a <- a - 1
      else b <- b + 1
    }else if(angleA <= tol_m & angleB > tol_m){
      b <- b + 1
      if(b == nrow(chrDf)) break
    }else if(angleA > tol_m & angleB <= tol_m){
      a <- a - 1
      if(a == 1) break
    }else if(angleA <= tol_m & angleB <= tol_m){
      break
    }
  }
  a_end <- a;b_end <- b
  chrDf <- chrDf[a_end:b_end, ]
  #plot_chrDf(chrDf, baseline = TRUE)
  return(chrDf)
}
edgeTrack_crude <- function(chrDf, preNum = 3, tol_m = 30){
  #browser()
  getLine <- function(A, B){
    slope <- (B[2] - A[2]) / (B[1] - A[1])
    intercept <- A[2] - slope * A[1]
    return(c(slope = slope, intercept = intercept))
  }
  apex_idx <- which.max(chrDf$intensity)
  int <- chrDf$intensity;rt <- chrDf$rt
  int_a <- int[1:apex_idx];int_b <- int[apex_idx:length(int)]
  while(!all(diff(int_a) >0) | !all(diff(int_b) < 0)){
    int <- smoothMean(int, size = 3)
    #plot(x = 1:length(int), y = int)
    apex_idx <- which.max(int)
    int_a <- int[1:apex_idx];int_b <- int[apex_idx:length(int)]
  }
  int_norm <- (int - min(int)) / (max(int) - min(int))
  rt_norm <- (rt - min(rt)) / (max(rt) - min(rt))
  # left
  a_start <- apex_idx - ((preNum - 1) / 2) - 1
  if(a_start <= 1) a_start <- 1
  b_start <- apex_idx + ((preNum - 1) / 2) + 1
  if(b_start >= nrow(chrDf)) b_start <- nrow(chrDf)
  a <- a_start;b <- b_start
  while(a!=1){
    a0 <- a - 1
    if(a == 1) a0 <- a  + 1
    A <- c(rt_norm[a], int_norm[a]);A0 <- c(rt_norm[a0], int_norm[a0])
    tmpA <- getLine(A, A0);slopeA <- tmpA["slope"];interceptA <- tmpA["intercept"]
    angleA0 <- atan(abs(slopeA)) * (180 / pi)
    if(angleA0 <= tol_m){break}
    else if(angleA0 > tol_m){a <- a - 1;if(a < 1){a <- 1}}
  }
  while(b!=nrow(chrDf)){
    b0 <- b + 1
    if(b == nrow(chrDf)) b0 <- b - 1
    B <- c(rt_norm[b], int_norm[b]);B0 <- c(rt_norm[b0], int_norm[b0])
    tmpB <- getLine(B, B0);slopeB <- tmpB["slope"];interceptB <- tmpB["intercept"]
    angleB0 <- atan(abs(slopeB)) * (180 / pi)
    if(angleB0 <= tol_m){break}
    else if(angleB0 > tol_m){b <- b + 1;if(b > nrow(chrDf)) b <- nrow(chrDf)}
  }
  a_end <- a;b_end <- b
  chrDf <- chrDf[a_end:b_end, ]
  return(chrDf)
}
cal_massTol <- function(chrDf, factor = 1, range = 0.95){
  delete <- round(nrow(chrDf) * (1-range))
  mz_ref <- chrDf$mz[which.max(chrDf$intensity)]
  delta <- abs(chrDf$mz - mz_ref) / (mz_ref)
  delta <- sort(delta);delta <- delta[1:(length(delta) - delete)]
  SD <- sd(delta)
  mass_tol <- factor * SD * 10^6
  return(mass_tol)
}
cal_cwtParam <- function(chrDfList, widthRange = 0.95){
  noiseVec <- sapply(1:length(chrDfList), function(i) {
    chrDf <- chrDfList[[i]]
    noise <- chrDf$baseline[which.max(chrDf$intensity)]
    return(noise)
  })
  noise <- min(noiseVec)
  intVec <- sapply(1:length(chrDfList), function(i) {
    chrDf <- chrDfList[[i]]
    int <- max(chrDf$intensity)
    return(int)
  })
  snVec <- intVec / noiseVec
  sn <- min(snVec)
  preNumVec <- sapply(1:length(chrDfList), function(i) {
    chrDf <- chrDfList[[i]]
    return(nrow(chrDf))
  })
  preNum <- min(preNumVec)
  peakWidthVec <- sapply(1:length(chrDfList), function(i) {
    chrDf <- chrDfList[[i]]
    peakWidth <- max(chrDf$rt) - min(chrDf$rt)
    return(peakWidth)
  })
  #browser()
  s <- sort(peakWidthVec)[round(length(peakWidthVec) * widthRange)]
  #plot(x = 1:length(peakWidthVec), y = peakWidthVec)
  peakWidthVec <- peakWidthVec[peakWidthVec <= s]
  peakWidth_min <- min(peakWidthVec)
  peakWidth_max <- max(peakWidthVec)
  return(c(noise = noise, sn = sn, preNum = preNum, peakWidth_min = peakWidth_min, peakWidth_max = peakWidth_max))
}
cal_shift <- function(chrDfList, rt_tol = 5, mz_tol = 0.02, tol_nf = 0.5, method = "mz", range = 0.95){
  #browser()
  tmp <- lapply(chrDfList, function(x) {
    sample <- attributes(x)$sample
    mz <- x$mz[which.max(x$intensity)]
    rt <- x$rt[which.max(x$intensity)]
    df <- dplyr::tibble(mz = mz, rt = rt, sample = sample)
    return(df)
  })
  tmp <- purrr::list_rbind(tmp)
  tmp <- tmp %>% dplyr::arrange(sample, mz, rt)
  sample_number <- length(unique(tmp$sample))
  sample_peakNum <- sapply(1:sample_number, function(i) {
    df <- tmp %>% dplyr::filter(sample == i)
    return(nrow(df))
  })
  ref_sample_idx <- which.max(sample_peakNum)
  align_sample_idx <- setdiff(1:sample_number, ref_sample_idx)
  ref_tmp <- tmp %>% dplyr::filter(sample == ref_sample_idx)
  align_tmp <- tmp %>% dplyr::filter(sample %in% align_sample_idx)
  pb <- utils::txtProgressBar(max = nrow(ref_tmp), style = 3)
  delta <- lapply(1:nrow(ref_tmp), function(i) {
    utils::setTxtProgressBar(pb, i)
    x <- ref_tmp[i, ]
    mz_ref <- x$mz;rt_ref <- x$rt
    align_df <- align_tmp %>% dplyr::filter(dplyr::near(mz, mz_ref, tol = mz_tol) & dplyr::near(rt, rt_ref, tol = rt_tol))
    if(nrow(align_df) < round((sample_number - 1) * tol_nf)) return(NULL)
    align_df_new <- lapply(unique(align_df$sample), function(j) {
      align_df_tmp <- align_df %>% dplyr::filter(sample == j)
      if(nrow(align_df_tmp) > 1){
        if(method == "mz") align_df_tmp <- align_df_tmp[which.min(align_df_tmp$mz - mz_ref),]
        else if(method == "rt") align_df_tmp <- align_df_tmp[which.min(align_df_tmp$rt - rt_ref),]
        return(align_df_tmp)
      }else{
        return(align_df_tmp)
      }
    })
    align_df_new <- purrr::list_rbind(align_df_new)
    df_all <- rbind(align_df_new, x)
    delta_mz <- max(df_all$mz) - min(df_all$mz)
    delta_rt <- max(df_all$rt) - min(df_all$rt)
    return(list(delta_mz = delta_mz, delta_rt = delta_rt))
  })
  delta <- delta[!sapply(delta, is.null)]
  delta_mz <- sapply(delta, function(x) {x$delta_mz})
  delta_rt <- sapply(delta, function(x) {x$delta_rt})
  mzwid <- sort(delta_mz)[round(length(delta_mz) * range)]
  bw <- sort(delta_rt)[round(length(delta_rt) * range)]
  return(c(mzwid = mzwid, bw = bw))
  # plot(x = 1:length(delta_mz), delta_mz)
  # plot(x = 1:length(delta_rt), delta_rt)
}

#' @title optParam4xcms
#' @description
#' optParam4xcms.
#'
#' @param data_QC A MsExperiment object which contains QC sample.
#' @param res_dir Results path.
#' @param bin bin size.
#' @param slide slide size.
#' @param mslevel mslevel.
#' @param thread1 thread1.
#' @param output_bins Whether to output bins.
#' @param bin_scanNum bin's scanNum threshold.
#' @param smoothPara smoothPrar.
#' @param baselinePara baselinePara
#' @param sn sn
#' @param preNum preNum.
#' @param tol_m tol_m.
#' @param snthresh snthresh.
#' @param thread2 thread2.
#'
#' @return A parameter vector.
#' @export
#'
#' @examples
#' optParam4xcms(data_QC = data_QC, thread2 = 4)
optParam4xcms <- function(data_QC, res_dir = "./",
                          bin = 0.05, slide = 0.05, mslevel = 1, thread1 = 1, output_bins = FALSE, bin_scanNum = 100,
                          smoothPara = get_smoothPara(), baselinePara = get_baselinePara(), sn = 3, preNum = 3, tol_m = 10, snthresh = 0.5, thread2 = 1,
                          maxMassTol = 100, minPeakWidth = 3, maxPeakWidth = 100,thread3 = 1){
  start_time <- Sys.time()
  message("You are using optParam4xcms function for ms1!")
  message("Generate mass bin...")
  chrDf_bins_QC <- lapply(1:length(data_QC), function(n) {
    message(paste0("QC ", n, "/", length(data_QC), "\n"))
    ndata <- data_QC[n]
    tmp <- generateBin(ndata = ndata, bin = bin, slide = slide, mslevel = mslevel, thread = thread1)
    return(tmp)
  })
  #chrDf_bins_QC <- readRDS("D:/fudan/Projects/2024/MetaboProcess/Progress/generate_bins/240705/chrDf_bins_QC.rds")
  for(n in 1:length(chrDf_bins_QC)){
    for(m in 1:length(chrDf_bins_QC[[n]])){
      attributes(chrDf_bins_QC[[n]][[m]])$sample <- n
    }
  }
  chrDf_bins <- purrr::list_flatten(chrDf_bins_QC)
  chrDf_bins <- chrDf_bins[which(sapply(chrDf_bins, function(x){
    if(nrow(x) <= bin_scanNum) return(FALSE)
    else return(TRUE)
  }))]
  message(paste0("You get ", length(chrDf_bins), " bins\n"))
  if(output_bins){
    filePath <- paste0(res_dir, "chrDf_bins.rds")
    saveRDS(chrDf_bins, file = filePath)
  }
  message("Generate ZOIs...") # i 6804
  pb <- utils::txtProgressBar(max = length(chrDf_bins), style = 3)
  if(thread2 == 1){
    ZOIList <- lapply(1:length(chrDf_bins), function(i) {
      utils::setTxtProgressBar(pb, i)
      #print(i)
      pickZOI(chrDf = chrDf_bins[[i]], smoothPara = smoothPara, baselinePara = baselinePara, sn = sn, preNum = preNum, tol_m = tol_m, snthresh = snthresh)
    })
  }else if(thread2 > 1){
    cl <- snow::makeCluster(thread2)
    doSNOW::registerDoSNOW(cl)
    opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
                                                                 n))
    envir <- environment(pickZOI)
    parallel::clusterExport(cl, varlist = ls(envir), envir = envir)
    ZOIList <- foreach::`%dopar%`(foreach::foreach(i = 1:length(chrDf_bins),
                                                   .options.snow = opts),
                                  {
                                    pickZOI(chrDf = chrDf_bins[[i]], smoothPara = smoothPara, baselinePara = baselinePara, sn = sn, preNum = preNum, tol_m = tol_m, snthresh = snthresh)
                                  })
    snow::stopCluster(cl)
    gc()
  }else{
    stop("wrong thread2!")
  }
  ZOIList <- purrr::list_flatten(ZOIList)
  ZOIList <- ZOIList[sapply(ZOIList, function(x) {
    if(is.null(x)) return(FALSE)
    else return(TRUE)
  })]
  #browser()
  #saveRDS(ZOIList, file = "D:/fudan/Projects/2024/MetaboProcess/Progress/optParam/240708/ZOIList.rds")
  #ZOIList <- readRDS("D:/fudan/Projects/2024/MetaboProcess/Progress/optParam/240708/ZOIList.rds")
  ZOIList <- lapply(ZOIList, function(x) tightenZOI(x, dight = 3))
  ZOIList <- ZOIList[sapply(ZOIList, function(x) !is.null(x))]
  # ZOIList <- lapply(ZOIList, function(x) checkShapeZOI(x, snthresh = 3))
  # ZOIList <- ZOIList[sapply(ZOIList, function(x) !is.null(x))]
  massTolVec <- sapply(ZOIList, function(x) {cal_massTol(chrDf = x, factor = 2, range = 1)})
  plot(1:length(massTolVec), massTolVec)
  peakWidthVec <- sapply(1:length(ZOIList), function(i) {
    chrDf <- ZOIList[[i]]
    peakWidth <- max(chrDf$rt) - min(chrDf$rt)
    return(peakWidth)
  })
  plot(1:length(peakWidthVec), peakWidthVec)
  ZOIList <- ZOIList[massTolVec <= maxMassTol & peakWidthVec >= minPeakWidth & peakWidthVec <= maxPeakWidth]
  message("You are generating a ZOITable...")
  ZOITable <- ZOIList2ZOITable(ZOIList)
  #saveRDS(ZOITable, file = "D:/fudan/Projects/2024/MetaboProcess/Progress/generate_bins/240705/ZOITable.rds")
  #ZOITable <-readRDS("D:/fudan/Projects/2024/MetaboProcess/Progress/generate_bins/240705/ZOITable.rds")
  message("You are finding isotopic peak...")
  isoTable <- FindIso(ZOITable, data_QC, thread = 5)
  #saveRDS(isoTable, file = "D:/fudan/Projects/2024/MetaboProcess/Progress/generate_bins/240705/isoTable.rds")
  #isoTable <- readRDS("D:/fudan/Projects/2024/MetaboProcess/Progress/generate_bins/240705/isoTable.rds")
  ZOIList_new <- ZOIList[as.integer(stringr::str_replace(rownames(isoTable), "CP", ""))]
  # 需要一个指标来判断峰是不是可以纳入参数计算
  message("Calculate parameters...")
  massTolVec <- sapply(ZOIList_new, function(x) {cal_massTol(chrDf = x, factor = 2, range = 1)})
  # plot(x = 1:length(massTolVec), y = sort(massTolVec))
  # plot(x = 1:length(diff(sort(massTolVec))), y = diff(sort(massTolVec)))
  tmp <- diff(sort(massTolVec))
  #sqrt(sum((massTolVec - median(massTolVec))^2) / (length(massTolVec) - 1))
  idx <- which(tmp > 2)[1]
  massTol_s <- mean(sort(massTolVec)[c(idx - 1, idx, idx + 1)])
  peakWidthVec <- sapply(1:length(ZOIList_new), function(i) {
    chrDf <- ZOIList_new[[i]]
    peakWidth <- max(chrDf$rt) - min(chrDf$rt)
    return(peakWidth)
  })
  tmp <- hist(peakWidthVec, 100, plot = FALSE)
  for(i in 1:length(tmp$density)){
    if(sum(tmp$density[order(tmp$density, decreasing = TRUE)[1:i]]) > 0.95) break
  }
  delta <- (tmp$mids[2] - tmp$mids[1]) / 2
  peakWidth_s_min <- min(tmp$mids[order(tmp$density, decreasing = TRUE)[1:i]]) - delta
  peakWidth_s_max <- max(tmp$mids[order(tmp$density, decreasing = TRUE)[1:i]]) + delta
  ZOIList_new <- ZOIList_new[massTolVec <= massTol_s & peakWidthVec >= peakWidth_s_min & peakWidthVec <= peakWidth_s_max]
  noiseVec <- sapply(1:length(ZOIList_new), function(i) {
    chrDf <- ZOIList_new[[i]]
    noise <- min(chrDf$baseline)
    return(noise)
  })
  noise <- min(noiseVec)
  intVec <- sapply(1:length(ZOIList_new), function(i) {
    chrDf <- ZOIList_new[[i]]
    int <- max(chrDf$intensity)
    return(int)
  })
  snVec <- intVec / noiseVec
  plot(1:length(snVec), sort(snVec))
  sn <- min(snVec)
  preNumVec <- sapply(1:length(ZOIList_new), function(i) {
    chrDf <- ZOIList_new[[i]]
    return(nrow(chrDf))
  })
  preNum <- min(preNumVec)
  shiftParam <- cal_shift(chrDfList = ZOIList_new, rt_tol = 5, mz_tol = 0.02, tol_nf = 0.5, method = "mz", range = 1)
  end_time <- Sys.time()
  print(end_time - start_time)
  i <- 24
  plot_chrDf(ZOIList_new[[which(peakWidthVec > 2)[i]]])
  xcms::peaksWithMatchedFilter(int = ZOIList_new[[which(peakWidthVec > 2)[i]]]$intensity,
                               rt = ZOIList_new[[which(peakWidthVec > 2)[i]]]$rt,
                               fwhm = 3,
                               snthresh = 1)
  return(c(ppm = tol_ppm, noise = noise, sn = sn, preNum = preNum, peakWidth_min = peakWidth_min, peakWidth_max = peakWidth_max,
           shiftParam))
}

#' @title optParam4xcms_old
#' @description
#' Find the optimal xcms parameters
#'
#' @param data_QC A MsExperiment object which contains QC sample.
#' @param res_dir Results path.
#' @param bin bin size.
#' @param mslevel mslevel.
#' @param thread1 thread1.
#' @param output_bins Whether to output bins.
#' @param bin_scanNum bin's scanNum threshold.
#' @param smooth smooth methods, c("mean", "sg").
#' @param size mean size.
#' @param p p.
#' @param etlD etlD.
#' @param IETH IETH.
#' @param preNum preNum.
#' @param sn sn.
#' @param loops loops.
#' @param tol_m1 tol_m1.
#' @param thread2 thread2.
#' @param threshold threshold.
#' @param tol_m2 tol_m2
#' @param output_ZOI Whether to output ZOIs.
#' @param output_ZOI2 Whether to output ZOI2s.
#' @param factor factor.
#' @param range range.
#' @param massRange massRange.
#' @param output_ZOI3 Whether to output ZOI3s.
#' @param rt_tol rt_tol.
#' @param mz_tol mz_tol.
#' @param tol_nf tol_nf.
#' @param method method.
#' @param shiftRange shiftRange.
#' @param widthRange widthRange.
#' @param maxMassTol maxMassTol.
#' @param minPeakWidth minPeakWidth.
#' @param maxPeakWidth maxPeakWidth.
#'
#' @return A list contains optimal parameters and ggplot2 object.
#' @export
#'
#' @examples
#' optParam4xcms(data_QC)
optParam4xcms_old <- function(data_QC, res_dir = "./",
                          bin = 0.05, mslevel = 1, thread1 = 1, output_bins = FALSE, bin_scanNum = 10,
                          smooth = "mean", size = 3, p = 3, etlD = 1, IETH = 10, preNum = 3, sn = 0, loops = 8, tol_m1 = 30, threshold = 1, thread2 = 1, tol_m2 = 40, output_ZOI = FALSE, output_ZOI2 = FALSE,
                          factor = 1, range = 1,massRange = 1, output_ZOI3 = FALSE,
                          rt_tol = 5, mz_tol = 0.02, tol_nf = 0.5, method = "mz", widthRange = 1, shiftRange = 0.95,
                          maxMassTol = 100, minPeakWidth = 1, maxPeakWidth = 100){
  #browser()
  start_time <- Sys.time()
  message("You are using optParam4xcms function for ms1!")
  message("Generate mass bin...")
  chrDf_bins_QC <- lapply(1:length(data_QC), function(n) {
    ndata <- data_QC[n]
    tmp <- generateBin(ndata = ndata, bin = bin, mslevel = mslevel, thread = thread1)
    return(tmp)
  })
  for(n in 1:length(chrDf_bins_QC)){
    for(m in 1:length(chrDf_bins_QC[[n]])){
      attributes(chrDf_bins_QC[[n]][[m]])$sample <- n
    }
  }
  chrDf_bins <- purrr::list_flatten(chrDf_bins_QC)
  chrDf_bins <- chrDf_bins[which(sapply(chrDf_bins, function(x){
    if(nrow(x) <= bin_scanNum) return(FALSE)
    else return(TRUE)
  }))]
  if(output_bins){
    filePath <- paste0(res_dir, "chrDf_bins.RData")
    save(chrDf_bins, file = filePath)
  }
  message("Generate ZOIs...")
  chrDf_ZOIs <- bins2ZOIs(chrDfList = chrDf_bins, smooth = smooth, size = size, p = p, etlD = etlD, loops = loops, tol_m = tol_m1, threshold = threshold, thread = thread2, sn = sn, IETH = IETH, preNum = preNum)
  chrDf_ZOIs_2 <- lapply(1:length(chrDf_ZOIs), function(i) {
    x <- chrDf_ZOIs[[i]]
    chrDf_tmp <- tryCatch({
      edgeTrack_finer(x, tol_m = tol_m2)
    }, error = function(e){
      edgeTrack_crude(x, tol_m = tol_m2)
    })
  })
  if(output_ZOI){
    filePath <- paste0(res_dir, "chrDf_ZOIs.RData")
    save(chrDf_ZOIs, file = filePath)
  }
  if(output_ZOI2){
    filePath <- paste0(res_dir, "chrDf_ZOIs_2.RData")
    save(chrDf_ZOIs_2, file = filePath)
  }
  message("Calculate parameters...")
  massTolVec <- sapply(chrDf_ZOIs_2, function(x) {cal_massTol(chrDf = x, factor = factor, range = range)})
  peakWidthVec <- sapply(1:length(chrDf_ZOIs_2), function(i) {
    chrDf <- chrDf_ZOIs_2[[i]]
    peakWidth <- max(chrDf$rt) - min(chrDf$rt)
    return(peakWidth)
  })
  massTol_s <- boxplot.stats(massTolVec[massTolVec < maxMassTol])$stats[5]
  peakWidth_s_min <- boxplot.stats(peakWidthVec[peakWidthVec<maxPeakWidth & peakWidthVec > minPeakWidth])$stats[1]
  peakWidth_s_max <- boxplot.stats(peakWidthVec[peakWidthVec<maxPeakWidth & peakWidthVec > minPeakWidth])$stats[5]
  #peakWidth_s <- boxplot.stats(peakWidthVec)$stats[5]
  #s <- sort(massTolVec)[round(length(massTolVec) * massRange)]
  chrDf_ZOIs_3 <- chrDf_ZOIs_2[which(massTolVec<=massTol_s & peakWidthVec<=peakWidth_s_max & peakWidthVec>=peakWidth_s_min)]
  if(output_ZOI3){
    filePath <- paste0(res_dir, "chrDf_ZOIs_3.RData")
    save(chrDf_ZOIs_3, file = filePath)
  }
  massTolVec <- massTolVec[which(massTolVec<=massTol_s & peakWidthVec<=peakWidth_s_max & peakWidthVec>=peakWidth_s_min)]
  tol_ppm <- sort(massTolVec)[round(length(massTolVec) * massRange)]
  cwtParam <- cal_cwtParam(chrDfList = chrDf_ZOIs_3, widthRange = widthRange)
  shiftParam <- cal_shift(chrDfList = chrDf_ZOIs_3, rt_tol = rt_tol, mz_tol = mz_tol, tol_nf = tol_nf, method = method, range = shiftRange)
  end_time <- Sys.time()
  print(end_time - start_time)
  return(c(c(ppm = tol_ppm), cwtParam, shiftParam))
}

