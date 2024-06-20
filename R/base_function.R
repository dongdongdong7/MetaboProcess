#' @title checkCentroid
#' @description
#' Check if the input file is centroid mode.
#'
#' @param file_path A vector of file's path.
#'
#' @return A logical vector.
#' @export
#'
#' @examples
#' file_dir <- "D:/fudan/Projects/2024/MetaboProcess/Progress/build_package/test_data/massProcesser/"
#' file_path <- list.files(file_dir, pattern = ".mzXML")
#' file_path <- paste0(file_dir, file_path)
#' checkCentroid(file_path = file_path)
checkCentroid <- function(file_path){
  sapply(file_path, function(x) {
    ifile <- x
    ifileh <- mzR::openMSfile(ifile, backend = NULL)
    logical <- all(unique(mzR::header(ifileh)$centroided))
    mzR::close(ifileh)
    rm(ifileh)
    return(logical)
  })
}
#' @title smoothMean
#' @description
#' Smooth function using mean.
#'
#' @param int intensity.
#' @param size mean size.
#'
#' @return intensity vector.
#' @export
#'
#' @examples
#' intensity <- c(1701.196, 2121.159, 2401.849, 2043.653, 1542.964,  723.803)
#' smoothMean(int = intensity)
smoothMean <- function(int, size = 3){
  if(size %% 2 == 0){
    stop("size should be a singular!")
  }
  half_size <- floor(size / 2)
  N <- length(int)
  special_idx <- c(seq_len(half_size), rev((N:(N - half_size + 1))))
  smooth_int <- sapply(1:N, function(i) {
    if(i %in% special_idx){
      return(int[i])
    }
    y <- sapply((i - half_size):(i + half_size), function(j){
      int[j]
    })
    smooth_int <- mean(y)
  })
  return(smooth_int)
}
#' @title smoothSg
#' @description
#' Smooth data with a Savitzky-Golay smoothing filter.
#'
#'
#' @param int intensity.
#' @param p filter order.
#' @param n filter length (must be odd).
#' @param m return the m-th derivative of the filter coefficients.
#' @param ts time scaling factor.
#'
#' @return intensity vector
#' @export
#'
#' @examples
#' intensity <- c(1701.196, 2121.159, 2401.849, 2043.653, 1542.964,  723.803)
#' smoothSg(int = intensity)
smoothSg <- function(int, p = 3, n = p + 3 - p%%2, m = 0, ts = 1){
  smooth_int <- signal::sgolayfilt(int, p = p, n = n, m = m, ts = ts)
  return(smooth_int)
}
#' @title baselineEs
#' @description
#' Estimating the baseline for a chrDf.
#'
#' @param chrDf A chrDf.
#' @param threshold threshold for inflect.
#' @param tol_m Tolerable angles.
#' @param loops loop number.
#'
#' @return A new chrDf with baseline column.
#' @export
#'
#' @examples
#' chrDf_j <- ZOIList[[j]]
#' chrDf_j <- baselineEs(chrDf_j, loops = 8, tol_m = 30)
baselineEs <- function(chrDf, threshold = 1, tol_m = 30, loops = 6){
  inflect <- function(x, threshold = 1){
    up   <- sapply(1:threshold, function(n) c(x[-(seq(n))], rep(NA, n)))
    down <-  sapply(-1:-threshold, function(n) c(rep(NA,abs(n)), x[-seq(length(x), length(x) - abs(n) + 1)]))
    a    <- cbind(x,up,down)
    list(minima = which(apply(a, 1, min) == a[,1]), maxima = which(apply(a, 1, max) == a[,1]))
  }
  fineEstimation <- function(xs, ys, int, tol_m = 0.5){
    top_idx <- inflect(ys, 1)$maxima
    bottom_idx <- inflect(int, 1)$minima
    bottom_idx <- setdiff(bottom_idx, top_idx)
    idx_seq_list <- lapply(1:(length(bottom_idx) - 1), function(x){
      idx_seq <- bottom_idx[x]:bottom_idx[x + 1]
    })
    idx_seq_logical <- rep(FALSE, length(idx_seq_list))
    for(i in 1:length(top_idx)){
      apex_idx <- top_idx[i]
      bottom_idx_tmp <- bottom_idx[which(int[bottom_idx] < ys[apex_idx])]
      idx_logical <- which(sapply(idx_seq_list, function(x){
        apex_idx %in% x
      }))
      if(length(idx_logical) == 0){
        next
      }
      if(!idx_seq_logical[idx_logical]){
        idx_seq_logical[idx_logical] <- TRUE
      }else{
        next
      }
      start_idx <- bottom_idx_tmp[which(apex_idx > bottom_idx_tmp)[length(which(apex_idx > bottom_idx_tmp))]]
      end_idx <- bottom_idx_tmp[which(apex_idx < bottom_idx_tmp)[1]]
      if(length(start_idx) == 0 | length(end_idx) == 0){
        next
      }
      if(is.na(start_idx) | is.na(end_idx)){
        next
      }
      xs_norm <- (xs - min(xs)) / (max(xs) - min(xs))
      ys_norm <- (ys - min(int)) / (max(int) - min(int))
      A <- c(xs_norm[apex_idx], ys_norm[apex_idx])
      B <- c(xs_norm[start_idx], ys_norm[start_idx])
      C <- c(xs_norm[end_idx], ys_norm[end_idx])
      Alength <- sqrt((B[1] - C[1])^2 + (B[2] - C[2])^2)
      Blength <- sqrt((A[1] - C[1])^2 + (A[2] - C[2])^2)
      Clength <- sqrt((B[1] - A[1])^2 + (B[2] - A[2])^2)
      COSA <- (Blength^2 + Clength^2 - Alength^2) / (2*Blength*Clength)
      COSB <- (Alength^2 + Clength^2 - Blength^2) / (2*Alength*Clength)
      COSC <- (Blength^2 + Alength^2 - Clength^2) / (2*Blength*Alength)
      angleA <- acos(COSA) * (180 / pi)
      angleB <- acos(COSB) * (180 / pi)
      angleC <- acos(COSC) * (180 / pi)
      if(angleB > tol_m | angleC > tol_m){
        a <- (ys[end_idx] - ys[start_idx]) / (xs[end_idx] - xs[start_idx])
        b <- (ys[start_idx] - a * xs[start_idx])
        ys[(start_idx + 1):(end_idx - 1)] <- a * xs[(start_idx + 1):(end_idx - 1)] + b
      }else{
        next
      }
    }
    return(ys)
  }
  rt <- chrDf$rt
  int <- chrDf$intensity
  #browser()
  bottom_idx <- inflect(int, threshold)$minima
  bottom_idx <- unique(c(1, bottom_idx, length(int)))
  if(length(bottom_idx) < 3){
    bottom_idx <- sort(unique(c(bottom_idx, order(int)[1:3])))
  }
  xs <- rt
  ys <- pracma::pchip(rt[bottom_idx], int[bottom_idx], xs)
  for(loop in 1:loops){
    #if(loop == 10) browser()
    tryCatch({ys <- fineEstimation(xs, ys, int, tol_m = tol_m)},
             error = function(e){
               return(ys <- ys)
             })
  }
  chrDf$baseline <- sapply(1:length(ys), function(i) {
    if(ys[i] > chrDf$intensity[i]) return(chrDf$intensity[i])
    else return(ys[i])
  })
  return(chrDf)
}
#' @title cal_IE
#' @description
#' Calculate IE.
#'
#'
#' @param int int.
#'
#' @return IE.
#' @export
#'
#' @examples
#' intensity <- c(1701.196, 2121.159, 2401.849, 2043.653, 1542.964,  723.803)
#' cal_IE(int = intensity)
cal_IE = function(int){
  len = length(int)
  int_max = max(int)
  int_max_idx = which(int == int_max)[1] # 因为可能有相同高度的int_max,选择第一个

  int_idx_a = 1:(int_max_idx - 1)
  int_idx_b = (int_max_idx + 1):len

  int_a = int[int_idx_a]
  int_b = int[int_idx_b]

  diff_a = diff(int_a)
  h_a = abs(diff_a[which(diff_a < 0)])
  diff_b = diff(int_b)
  h_b = abs(diff_b[which(diff_b > 0)])
  h = c(h_a,int_max,h_b)
  p = h / sum(h)
  S = -sum(p * log(p))

  return(S)
}
#' @title getChromPeakTable
#' @description
#' Get a chromPeakTable tibble.
#'
#' @param data A xcmsExperiment object after peak picking.
#' @param style c("xcms", "neatms")
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' getChromPeakTable(data)
getChromPeakTable <- function(data, style = "xcms"){
  chromPeakTable <- dplyr::as_tibble(cbind(xcms::chromPeaks(data),
                                           xcms::chromPeakData(data)),
                                     rownames = "cpid")
  if(style == "xcms") return(chromPeakTable)
  else if(style == "neatms"){
    sampleData <- MsExperiment::sampleData(data)
    sampleData$sample_order <- 1:nrow(sampleData)
    sampleData <- dplyr::as_tibble(sampleData) %>%
      dplyr::select(sample_order, sample_name)
    colnames(sampleData) <- c("sample", "sample_name")
    chromPeakTable <- dplyr::left_join(chromPeakTable, sampleData, by = c("sample" = "sample")) %>%
      dplyr::select(mz, mzmin, mzmax, rt, rtmin, rtmax, into, intb, maxo, sn, sample, sample_name)
    chromPeakTable$sample_name <- paste0(chromPeakTable$sample_name, ".mzML")
    return(chromPeakTable)
  }else stop("style is wrong!")
}
#' @title getChromPeak
#' @description
#' Get a chromPeak from XcmsExperiment.
#'
#' @param data A XcmsExperiment object.
#' @param cpid cpid.
#' @param expandRt expandRt.
#' @param smooth Whether or not smooth.
#' @param method Smooth method.
#' @param size size.
#' @param p p.
#' @param n n.
#' @param m m.
#' @param ts ts.
#'
#' @return A chrDf.
#' @export
#'
#' @examples
#' getChromPeak(data, cpid = cpid)
getChromPeak <- function(data, cpid, expandRt = 0,
                         smooth = FALSE, method = "mean",size = 3, p = 3, n = p + 3 - p%%2, m = 0, ts = 1){
  #browser()
  #chrs <- xcms::chromPeakChromatograms(data)
  #chrs <- xcms::chromPeakChromatograms(data, peaks = cpid)
  pk_mat <- xcms::chromPeaks(data)[cpid, ]
  pk_meta <- xcms::chromPeakData(data)[cpid, ]
  dataOrigin <- MsExperiment::sampleData(data)$spectraOrigin[pk_mat[["sample"]]]
  sps <- xcms::spectra(data) %>%
    Spectra::filterDataOrigin(dataOrigin) %>%
    Spectra::filterMsLevel(pk_meta[1, "ms_level"]) %>%
    Spectra::filterRt(c(pk_mat[["rtmin"]] - expandRt, pk_mat[["rtmax"]] + expandRt)) %>%
    Spectra::filterMzRange(c(pk_mat[["mzmin"]], pk_mat[["mzmax"]]))
  rtsps <- Spectra::rtime(sps)
  peaksData <- Spectra::peaksData(sps)
  chrDf <- purrr::list_rbind(lapply(1:length(sps), function(j){
    rt_tmp <- rtsps[j]
    peakMat <- peaksData[[j]]
    if(nrow(peakMat) == 0) return(NULL)
    else{
      peakDf <- as.data.frame(peakMat)
      peakDf$rt <- rt_tmp
      return(peakDf)
    }
  }))
  if(smooth){
    if(method == "mean") chrDf$intensity <- smoothMean(chrDf$intensity, size = size)
    else if(method == "sg") chrDf$intensity <- smoothSg(chrDf$intensity, p = p, n = n, m = m, ts = ts)
    else stop("method is wrong!")
  }
  return(chrDf)
  #plot_chrDf(chrDf)
}
#' @title getFeatureTable
#' @description
#' Get feature table.
#'
#' @param data A XcmsExperiment object.
#' @param method See xcms::featureValues function.
#' @param value See xcms::featureValues function.
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' fearureTable <- getFeatureTable(data)
getFeatureTable <- function(data, method = "medret", value = "intb"){
  #browser()
  chromPeaks <- xcms::chromPeaks(data)
  cpidVec <- rownames(chromPeaks)
  featureTable <- dplyr::as_tibble(xcms::featureDefinitions(data), rownames = "ftid")
  peakidList <- lapply(1:nrow(featureTable), function(i) {
    peakidx <- featureTable[i, ]$peakidx[[1]]
    return(cpidVec[peakidx])
  })
  featureTable$peakid <- peakidList
  featureValues <- dplyr::as_tibble(xcms::featureValues(data, method = method, value = value))
  colnames(featureValues) <- stringr::str_replace(colnames(featureValues),
                                                  pattern = ".mzXML|.mzxml|.mzML|.mzml",
                                                  replacement = "")
  featureTable <- dplyr::as_tibble(cbind(featureTable, featureValues))
  return(featureTable)
}
