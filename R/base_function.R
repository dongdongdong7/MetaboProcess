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
    bottom_idx <- unique(c(1, bottom_idx, length(int)))
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

#' @title getChromPeaksDf
#' @description
#' Get chromPeak date frame.
#'
#' @param data A XcmsExperiment object.
#' @param cpid cpid vector.
#' @param noise1 noise for ms1.
#' @param noise2 noise for ms2.
#' @param smoothPara smoothPara.
#' @param expandRt expandRt.
#' @param expandMz expandMz
#'
#' @return A chrDfList.
#' @export
#'
#' @examples
#' getChromPeaksDf(data)
getChromPeaksDf <- function(data, cpid = NA, noise1 = 0, noise2 = 0, smoothPara = get_smoothPara(), expandRt = 0, expandMz = 0){
  tmp <- xcms::chromPeakChromatograms(data, peaks = cpid, expandRt = expandRt, expandMz = expandMz)
  chrDfList <- lapply(1:nrow(tmp), function(i) {
    XChr <- tmp[i, 1]
    msLevel <- XChr@msLevel
    if(msLevel == 1) noise <- noise1
    else if(msLevel == 2) noise <- noise2
    else stop("msLevel wrong!")
    intensity <- XChr@intensity
    intensity[which(intensity <= noise)] <- NA
    rt <- XChr@rtime
    chrDf <- dplyr::tibble(intensity = intensity, rt = rt) %>%
      dplyr::filter(!is.na(intensity))
    chrDf$mz <- mean(XChr@mz)
    if(smoothPara$smooth){
      if(smoothPara$method == "mean") chrDf$intensity <- smoothMean(chrDf$intensity, size = smoothPara$size)
      else if(smoothPara$method == "sg") chrDf$intensity <- smoothSg(chrDf$intensity, p = smoothPara$p, n = smoothPara$n, m = smoothPara$m, ts = smoothPara$ts)
    }
    return(chrDf)
  })
  names(chrDfList) <- cpid
  return(chrDfList)
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

#' @title peakAnnotation.
#'
#' @param data A XcmsExperiment object.
#' @param polarity positive or negative
#' @param adinfo adinfo <- data(positive.adinfo, package = "cliqueMS")
#' @param thread thread.
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' data(positive.adinfo, package = "cliqueMS")
#' positive.adinfo <- positive.adinfo[positive.adinfo$adduct %in%
#'                                      c("[M+H]+", "[M+H-H2O]+", "[M+Na]+", "[M+H-NH3]+", "[M+K]+", "[M+NH4]+"),]
#' resTable <- peakAnnotation(data = data_filter, polarity = "positive", adinfo = positive.adinfo, thread = length(data_filter))
peakAnnotation <- function(data, polarity = "positive", adinfo, thread = 1){
  #browser()
  iterList <- lapply(1:length(data), function(i) {
    data_n <- data[i]
    attributes(data_n)$order <- i
    return(data_n)
  })
  loop <- function(data_n){
    data_xs <- xcms:::.XCMSnExp2xcmsSet(data_n)
    set.seed(2)
    ex.cliqueGroups <- cliqueMS::getCliques(data_xs, filter = TRUE, silent = FALSE)
    ex.Isotopes <- cliqueMS::getIsotopes(ex.cliqueGroups, ppm = 10)
    ex.Adducts <- cliqueMS::getAnnotation(ex.Isotopes, ppm = 10,
                                          adinfo = adinfo, polarity = polarity,
                                          normalizeScore = TRUE)
    resTable <- ex.Adducts@peaklist[, c("mz", "rt", "maxo", "sample", "cliqueGroup", "isotope", "mass1", "an1")]
    return(resTable)
  }
  pb <- utils::txtProgressBar(max = length(data), style = 3)
  if(thread == 1){
    resList <- lapply(1:length(iterList), function(i) {
      utils::setTxtProgressBar(pb, i)
      data_n <- iterList[[i]]
      loop(i)
    })
  }else if(thread > 1){
    cl <- snow::makeCluster(thread)
    doSNOW::registerDoSNOW(cl)
    opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
                                                                 n))
    resList <- foreach::`%dopar%`(foreach::foreach(i = iterList,
                                                   .packages = c("cliqueMS", "xcms"),
                                                   .options.snow = opts),
                                  {
                                    loop(i)
                                  })
    snow::stopCluster(cl)
    gc()
  }else stop("Thread wrong!")
  resTable <- purrr::list_rbind(resList)
  resTable <- dplyr::as_tibble(resTable, rownames = "cpid")
  return(resTable)
}
