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
