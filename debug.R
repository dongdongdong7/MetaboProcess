chrDf_test <- readRDS("D:/fudan/Projects/2024/MetaboProcess/Progress/build_package/debug/240627/test_chrDf.rds")
plot_chrDf(chrDf_test)
chrDf_test <- baselineEs(chrDf = chrDf_test, loops = 8, tol_m = 30, threshold = 1)
plot_chrDf(chrDf_test, baseline = TRUE)
#MetaboProcess::plot_chrDf(chrDf_j, baseline = TRUE)
chrDf_j_new <- MetaboProcess::edgeTrack(chrDf = chrDf_j)

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
FindZOI_baseline(chrDf = chrDf_test, noise = 2378)
