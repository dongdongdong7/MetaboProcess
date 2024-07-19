devtools::document()
# Amino
file_dir <- "D:/fudan/Projects/2024/MetaboProcess/Progress/build_package/test_data/mrm_amino/"
patterns <- c(".mzXML", ".mzxml", ".mzML", ".mzml")
patterns <- paste0(patterns, collapse = "|")
file_path <- list.files(file_dir, pattern = patterns)
file_path <- paste0(file_dir, file_path)
raw_data <- MSnbase::readSRMData(files = file_path)

i <- 75;j <- 3
chr <- raw_data[i, j]
chrDf_bin <- dplyr::tibble(mz = chr@productMz[1], intensity = chr@intensity, rt = chr@rtime * 60)
plot_chrDf(chrDf_bin)
pickZOI <- function(chrDf, smoothPara = get_smoothPara(), baselinePara = get_baselinePara(),
                    sn = 1, preNum = 3, tol_m = 10, snthresh = 0.5){ # chrDf 通常需要是一个bin窗口
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
    sn_tmp <- max(x$intensity) / mean(x$baseline)
    if(sn_tmp > sn) return(TRUE)
    else return(FALSE)
  })]
  if(length(ZOIList) == 0) return(NULL)
  return(ZOIList)
}
ZOIList <- pickZOI(chrDf = chrDf_bin, sn = 3)
plot_chrDf(ZOIList[[1]], baseline = TRUE)

pdf("D:/fudan/Projects/2024/MetaboProcess/Progress/MRM_test/240628/test.pdf")
for(i in 1:nrow(raw_data)){
  for(j in 1:ncol(raw_data)){
    print(i)
    print(j)
    chr <- raw_data[i, j]
    chrDf_bin <- dplyr::tibble(mz = chr@productMz[1], intensity = chr@intensity, rt = chr@rtime * 60)
    ZOIList <- pickZOI(chrDf = chrDf_bin, sn = 3)
    noise0 <- noiseEstimation(chrDf_bin)
    baselinePara = get_baselinePara()
    chrDf <- baselineEs(chrDf = chrDf_bin, threshold = baselinePara$threshold, tol_m = baselinePara$tol_m, loops = baselinePara$loops)
    p0 <- plot_chrDf(chrDf = chrDf, baseline = TRUE, noise = noise0)
    edgeRt <- purrr::list_c(lapply(ZOIList, function(x) {
      c(x$rt[1], x$rt[nrow(x)])
    }))
    if(is.null(edgeRt)){
      p <- p0 +
        ggplot2::labs(title = paste0("i: ", i, "-", "j: ", j))
    }else{
      p <- p0 +
        ggplot2::geom_vline(xintercept = edgeRt) +
        ggplot2::labs(title = paste0("i: ", i, "-", "j: ", j))
    }
    print(p)
  }
}
dev.off()
