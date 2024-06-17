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
#' @param widthRangewidthRange.
#'
#' @return A list contains optimal parameters and ggplot2 object.
#' @export
#'
#' @examples
#' optParam4xcms(data_QC)
optParam4xcms <- function(data_QC, res_dir = "./",
                          bin = 0.05, mslevel = 1, thread1 = 1, output_bins = FALSE, bin_scanNum = 10,
                          smooth = "mean", size = 3, p = 3, etlD = 1, IETH = 10, preNum = 3, sn = 0, loops = 8, tol_m1 = 30, threshold = 1, thread2 = 1, tol_m2 = 40, output_ZOI = FALSE, output_ZOI2 = FALSE,
                          factor = 1, range = 0.95,massRange = 0.95, output_ZOI3 = FALSE,
                          rt_tol = 5, mz_tol = 0.02, tol_nf = 0.5, method = "mz", widthRange = 0.95, shiftRange = 0.95){
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
  massTol_s <- boxplot.stats(massTolVec)$stats[5]
  peakWidth_s <- boxplot.stats(peakWidthVec)$stats[5]
  #s <- sort(massTolVec)[round(length(massTolVec) * massRange)]
  chrDf_ZOIs_3 <- chrDf_ZOIs_2[which(massTolVec<massTol_s & peakWidthVec<peakWidth_s)]
  if(output_ZOI3){
    filePath <- paste0(res_dir, "chrDf_ZOIs_3.RData")
    save(chrDf_ZOIs_3, file = filePath)
  }
  massTolVec <- massTolVec[which(massTolVec<massTol_s & peakWidthVec<peakWidth_s)]
  tol_ppm <- sort(massTolVec)[round(length(massTolVec) * massRange)]
  cwtParam <- cal_cwtParam(chrDfList = chrDf_ZOIs_3, widthRange = widthRange)
  shiftParam <- cal_shift(chrDfList = chrDf_ZOIs_3, rt_tol = rt_tol, mz_tol = mz_tol, tol_nf = tol_nf, method = method, range = shiftRange)
  end_time <- Sys.time()
  print(end_time - start_time)
  return(c(c(ppm = tol_ppm), cwtParam, shiftParam))
}

