#' @title paramounterPart1
#' @description
#' Paramounter part1 to calculate ppmCut.
#'
#' @param file_path Files path of QC sample.
#' @param massSDrange massSDrange.
#' @param smooth smooth.
#' @param cutoff cutoff.
#' @param thread thread.
#' @param msLevel msLevel.
#'
#' @return A number - ppmCut.
#' @export
#'
#' @examples
#' ppmCut <- paramounterPart1(file_path = file_path, thread = 24)
paramounterPart1 <- function(file_path, massSDrange = 2, smooth = 0, cutoff = 0.95, thread = 1, msLevel = 1L){
  start_time <- Sys.time()
  ppm2DList <- lapply(1:length(file_path), function(q) {
    message(paste0(q, "/", length(file_path), "\n"))
    ms1data <- MSnbase::readMSData(files = file_path[q], mode = "onDisk", msLevel. = msLevel)
    mzTmp <- MSnbase::mz(ms1data)
    mzRange <- c(min(unlist(mzTmp)), max(unlist(mzTmp)))
    ROI <- seq(mzRange[1], mzRange[2], 0.05)
    mzData <- mzTmp
    intData <- MSnbase::intensity(ms1data)
    rtime <- MSnbase::rtime(ms1data)

    peak_smooth <- function(x,level=smooth){
      n <- level
      if(length(x) < 2*n){
        return(x)
      } else if(length(unique(x))==1){
        return(x)
      } else{
        y <- vector(length=length(x))
        for(i in 1:n){
          y[i] <- sum(c((n-i+2):(n+1),n:1)*x[1:(i+n)])/sum(c((n-i+2):(n+1),n:1))
        }
        for(i in (n+1):(length(y)-n)){
          y[i] <-  sum(c(1:(n+1),n:1)*x[(i-n):(i+n)])/sum(c(1:(n+1),n:1))
        }
        for(i in (length(y)-n+1):length(y)){
          y[i] <- sum(c(1:n,(n+1):(n+i-length(x)+1))*x[(i-n):length(x)])/sum(c(1:n,(n+1):(n+i-length(x)+1)))
        }
        return(y)
      }
    }
    calROI <- function(i){
      # Obtain data lists in each m/z bin
      #browser()
      ppm2Ddist <- as.data.frame(matrix(ncol = 3, nrow = 1))
      colnames(ppm2Ddist) <- c("mz", "rt", "ppm")
      currmzRange <- c(ROI[i], ROI[i+1])
      tmpMZdata <- mzData
      tmpINTdata <- intData
      for(j in 1:length(mzData)){
        index <- which(tmpMZdata[[j]] >= currmzRange[1] & tmpMZdata[[j]] < currmzRange[2])
        tmpMZdata[[j]] <- tmpMZdata[[j]][index]
        tmpINTdata[[j]] <- tmpINTdata[[j]][index]
      }
      # Extract the intensity vectors from each m/z bin
      eicINTraw <- c()
      eicINT <- c()
      eicRT <- c()
      for(k in 1:length(mzData)){
        if(length(tmpINTdata[[k]]) > 0){
          eicINTraw[k] <- mean(tmpINTdata[[k]])
        }else{
          eicINTraw[k] <- 0
        }
        eicRT[k] <- rtime[k]
      }
      if(sum(eicINTraw != 0) == 0) return(ppm2Ddist)
      # Sort the intensity vectors from each m/z bin, estimate the noise cut off and average
      eicINT <- peak_smooth(eicINTraw)
      eicNon0 <- sort(eicINT[eicINT > 0])
      if(length(eicNon0) > 10){
        for(x in seq(10,length(eicNon0), 10)){
          sd <- sd(eicNon0[1:x])
          blk <- sum(eicNon0[1:x])/x
          thres <- blk + 3*sd
          if(x+1 <= length(eicNon0)){
            if(eicNon0[x+1] >= thres) break()
          }
        }
        cutOFF <- eicNon0[x]
      }else{
        cutOFF <- max(eicNon0)
      }

      aboveTHindex <- which(eicINT > cutOFF)
      if(length(aboveTHindex) == 0) return(ppm2Ddist)
      candidateSegInd <- split(aboveTHindex, cumsum(c(1, diff(aboveTHindex) != 1)))
      peakInd <- c()
      for(x in 1:length(candidateSegInd)){
        peakInd[x] <- which(eicINT[candidateSegInd[[x]]] == max(eicINT[candidateSegInd[[x]]]))[1] + min(candidateSegInd[[x]]) - 1
      }
      refMZvec <- c()
      for(y in 1:length(peakInd)){
        highestINT <- which(tmpINTdata[[peakInd[y]]] == max(tmpINTdata[[peakInd[y]]]))[1]
        refMZvec[y] <- tmpMZdata[[peakInd[y]]][highestINT]
      }

      # Estimate the universal parameters (mass tolerance, peak height, and peak width) for each m/z bin
      ppmDiff <- c()
      for(z in 1:length(peakInd)){
        currPeakInd <- peakInd[z]
        currRefMz <- refMZvec[z]
        currSamePeakMass <- c()
        currSamePeakMass <- c(currSamePeakMass, currRefMz)
        leftInd <- currPeakInd-1
        rightInd <- currPeakInd+1
        if(leftInd > 0){
          while (length(tmpMZdata[[leftInd]]) > 0 & mean(tmpINTdata[[leftInd]]) >= cutOFF) {
            if (length(tmpMZdata[[leftInd]]) == 1){
              currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[leftInd]])
              if(eicINT[leftInd] > eicINT[leftInd+1] & length(currSamePeakMass) > 5){
                Q1 <- as.numeric(summary(currSamePeakMass)[2])
                Q3 <- as.numeric(summary(currSamePeakMass)[5])
                LB <- Q1 - 1.5 *(Q3 - Q1)
                RB <- Q3 + 1.5 *(Q3 - Q1)
                if (currSamePeakMass[length(currSamePeakMass)] < LB || currSamePeakMass[length(currSamePeakMass)] > RB) break()
              }
            } else {
              abvector <- abs(tmpMZdata[[leftInd]] - currRefMz)
              NearInd <- which(abvector == min(abvector))[1]
              currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[leftInd]][NearInd])
              if(eicINT[leftInd] > eicINT[leftInd+1] & length(currSamePeakMass) > 5){
                Q1 <- as.numeric(summary(currSamePeakMass)[2])
                Q3 <- as.numeric(summary(currSamePeakMass)[5])
                LB <- Q1 - 1.5 *(Q3 - Q1)
                RB <- Q3 + 1.5 *(Q3 - Q1)
                if (currSamePeakMass[length(currSamePeakMass)] < LB || currSamePeakMass[length(currSamePeakMass)] > RB) break()
              }
            }
            leftInd <- leftInd-1
            if(leftInd <= 0) break()
          }
        }
        if(rightInd <= length(tmpMZdata)){
          while (length(tmpMZdata[[rightInd]]) > 0 & mean(tmpINTdata[[rightInd]]) >= cutOFF) {
            if (length(tmpMZdata[[rightInd]]) == 1){
              currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[rightInd]])
              if(eicINT[rightInd] > eicINT[rightInd-1] & length(currSamePeakMass) > 5){
                Q1 <- as.numeric(summary(currSamePeakMass)[2])
                Q3 <- as.numeric(summary(currSamePeakMass)[5])
                LB <- Q1 - 1.5 *(Q3 - Q1)
                RB <- Q3 + 1.5 *(Q3 - Q1)
                if (currSamePeakMass[length(currSamePeakMass)] < LB || currSamePeakMass[length(currSamePeakMass)] > RB) break()
              }
            } else {
              abvector <- abs(tmpMZdata[[rightInd]] - currRefMz)
              NearInd <- which(abvector == min(abvector))[1]
              currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[rightInd]][NearInd])
              if(eicINT[rightInd] > eicINT[rightInd-1] & length(currSamePeakMass) > 5){
                Q1 <- as.numeric(summary(currSamePeakMass)[2])
                Q3 <- as.numeric(summary(currSamePeakMass)[5])
                LB <- Q1 - 1.5 *(Q3 - Q1)
                RB <- Q3 + 1.5 *(Q3 - Q1)
                if (currSamePeakMass[length(currSamePeakMass)] < LB || currSamePeakMass[length(currSamePeakMass)] > RB) break()
              }
            }
            rightInd <- rightInd+1
            if(rightInd > length(tmpMZdata)) break()
          }
        }

        if(length(currSamePeakMass) > 1){
          ppmDiff[z] <- (massSDrange*sd(currSamePeakMass))/currRefMz * 1e6
          ppm2Ddist <- rbind(ppm2Ddist, c(currRefMz, rtime[[peakInd[z]]], ppmDiff[z]))
        }
      }
      return(ppm2Ddist)
    }

    pb <- utils::txtProgressBar(max = (length(ROI) - 1), style = 3)
    if(thread == 1){
      ppm2Ddist <- purrr::list_rbind(lapply(1:(length(ROI) - 1), function(i) {
        #print(i)
        utils::setTxtProgressBar(pb, i)
        calROI(i)
      }))
    }else if(thread > 1){
      ppmList <- BiocParallel::bplapply(1:(length(ROI) - 1), function(i) {
        calROI(i)
      }, BPPARAM = BiocParallel::SnowParam(workers = thread,
                                           progressbar = TRUE))
      ppm2Ddist <- purrr::list_rbind(ppmList)
    }else{
      stop("thread is wrong!")
    }
    return(ppm2Ddist)
  })
  ppm2D <- purrr::list_rbind(ppm2DList)
  ppm2D <- ppm2D[complete.cases(ppm2D),]
  ppm2D <- ppm2D[order(ppm2D[,3]),]
  ppm2D <- ppm2D[1:round(nrow(ppm2D)*0.97),]
  #plot(ppm2D$mz, ppm2D$ppm, ylab = "ppm", xlab = "m/z", pch=1, cex.main=4, cex.lab=1.7, cex.axis=2)
  ppm2Ddash <- ppm2D[1:round(nrow(ppm2D)*cutoff),]
  dashline <- max(ppm2Ddash[,3])
  long <- length(ppm2D[,3])
  cutoffvalue <- rep(dashline,long)
  #lines(ppm2D$mz, cutoffvalue, lty = "dashed", lwd = "3", col = "red")
  print(Sys.time() - start_time)
  return(dashline)
}

#' @title paramounterPart2
#' @description
#' Paramounter function part2 for xcms parameters.
#'
#' @param file_path Files path of QC sample.
#' @param massSDrange massSDrange.
#' @param smooth smooth.
#' @param ppmCut ppmCut calculated by part1.
#' @param thread thread.
#' @param msLevel msLevel.
#'
#' @return A vector of XCMS parameters.
#' @export
#'
#' @examples
#' parameters <- paramounterPart2(file_path = file_path, ppmCut = ppmCut, thread = 16)
paramounterPart2 <- function(file_path, massSDrange = 2, smooth = 0, ppmCut = NA, thread = 1, msLevel = 1L){
  if(is.na(ppmCut)) stop("You should run paramounterPart1 first!")

  mzDiff <- c()
  ppm <- c()
  noiselevel <- c()
  peakWidth <- c()
  peakScans <- c()
  SNRatio <- c()
  peakHeight <- c()
  massShiftALL <- c()
  rtShiftALL <- c()
  mzDiff2D <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(mzDiff2D) <- c("mz", "rt", "mzdiff")
  ppm2D <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(ppm2D) <- c("mz", "rt", "ppm")
  Shift <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(Shift) <- c("mz", "rt", "Sample")
  for(q in 1:length(file_path)){
    message(paste0(q, "/", length(file_path), "\n"))
    ms1data <- MSnbase::readMSData(files = file_path[q], mode = "onDisk", msLevel. = msLevel)
    mzTmp <- MSnbase::mz(ms1data)
    mzRange <- c(min(unlist(mzTmp)), max(unlist(mzTmp)))
    ROI <- seq(mzRange[1], mzRange[2], 0.05)
    mzData <- mzTmp
    intData <- MSnbase::intensity(ms1data)
    rtime <- MSnbase::rtime(ms1data)

    peak_smooth <- function(x,level=smooth){
      n <- level
      if(length(x) < 2*n){
        return(x)
      } else if(length(unique(x))==1){
        return(x)
      } else{
        y <- vector(length=length(x))
        for(i in 1:n){
          y[i] <- sum(c((n-i+2):(n+1),n:1)*x[1:(i+n)])/sum(c((n-i+2):(n+1),n:1))
        }
        for(i in (n+1):(length(y)-n)){
          y[i] <-  sum(c(1:(n+1),n:1)*x[(i-n):(i+n)])/sum(c(1:(n+1),n:1))
        }
        for(i in (length(y)-n+1):length(y)){
          y[i] <- sum(c(1:n,(n+1):(n+i-length(x)+1))*x[(i-n):length(x)])/sum(c(1:n,(n+1):(n+i-length(x)+1)))
        }
        return(y)
      }
    }
    calROI <- function(i){
      snALL <- c()
      sALL <- c()
      peakWidthALL <- c()
      peakScansALL <- c()
      ShiftTable <- as.data.frame(matrix(ncol = 3, nrow = 1))
      colnames(ShiftTable) <- c("mz", "rt", "Sample")
      ppm2Ddist <- as.data.frame(matrix(ncol = 3, nrow = 1))
      colnames(ppm2Ddist) <- c("mz", "rt", "ppm")
      mzdiff2Ddist <- as.data.frame(matrix(ncol = 3, nrow = 1))
      colnames(mzdiff2Ddist) <- c("mz", "rt", "mzdiff")
      noiseALL <- c()
      blankList <- list(snALL = snALL, sALL = sALL, peakWidthALL = peakWidthALL, peakScansALL = peakScansALL,
                        ShiftTable = ShiftTable, ppm2Ddist = ppm2Ddist, mzdiff2Ddist = mzdiff2Ddist, noiseALL = noiseALL)
      # Obtain data lists in each m/z bin
      currmzRange <- c(ROI[i], ROI[i+1])
      tmpMZdata <- mzData
      tmpINTdata <- intData
      for(j in 1:length(mzData)){
        index <- which(tmpMZdata[[j]] >= currmzRange[1] & tmpMZdata[[j]] < currmzRange[2])
        tmpMZdata[[j]] <- tmpMZdata[[j]][index]
        tmpINTdata[[j]] <- tmpINTdata[[j]][index]
      }
      # Extract the intensity vectors from each m/z bin
      eicINTraw <- c()
      eicINT <- c()
      eicRT <- c()
      for(k in 1:length(mzData)){
        if(length(tmpINTdata[[k]]) > 0){
          eicINTraw[k] <- mean(tmpINTdata[[k]])
        }else{
          eicINTraw[k] <- 0
        }
        eicRT[k] <- rtime[k]
      }
      if(sum(eicINTraw != 0) == 0) return(blankList)
      # Sort the intensity vectors from each m/z bin, estimate the noise cut off and average
      eicINT <- peak_smooth(eicINTraw)
      eicNon0 <- sort(eicINT[eicINT > 0])
      if(length(eicNon0) > 10){
        for(x in seq(10,length(eicNon0), 10)){
          sd <- sd(eicNon0[1:x])
          blk <- sum(eicNon0[1:x])/x
          thres <- blk + 3*sd
          if(x+1 <= length(eicNon0)){
            if(eicNon0[x+1] >= thres) break()
          }
        }
        cutOFF <- eicNon0[x]
      }else{
        cutOFF <- max(eicNon0)
        sd <- 0
        blk <- 0
      }
      noiseALL[i] <- cutOFF


      # Find the Reference m/z in each ROI from each m/z bin
      aboveTHindex <- which(eicINT > cutOFF)
      if(length(aboveTHindex) == 0) return(blankList)
      candidateSegInd <- split(aboveTHindex, cumsum(c(1, diff(aboveTHindex) != 1)))
      peakInd <- c()
      for(x in 1:length(candidateSegInd)){
        peakInd[x] <- which(eicINT[candidateSegInd[[x]]] == max(eicINT[candidateSegInd[[x]]]))[1] + min(candidateSegInd[[x]]) - 1
      }
      refMZvec <- c()
      for(y in 1:length(peakInd)){
        highestINT <- which(tmpINTdata[[peakInd[y]]] == max(tmpINTdata[[peakInd[y]]]))[1]
        refMZvec[y] <- tmpMZdata[[peakInd[y]]][highestINT]
      }

      # Estimate the universal parameters (mass tolerance, peak height, and peak width) for each m/z bin
      for(z in 1:length(peakInd)){
        currPeakInd <- peakInd[z]
        currRefMz <- refMZvec[z]
        currSamePeakMass <- c()
        currSamePeakMass <- c(currSamePeakMass, currRefMz)
        leftInd <- currPeakInd-1
        rightInd <- currPeakInd+1
        if(leftInd > 0){
          while (length(tmpMZdata[[leftInd]]) > 0 & mean(tmpINTdata[[leftInd]]) >= cutOFF) {
            if (length(tmpMZdata[[leftInd]]) == 1){
              currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[leftInd]])
              if(eicINT[leftInd] > eicINT[leftInd+1] & length(currSamePeakMass) > 5){
                Q1 <- as.numeric(summary(currSamePeakMass)[2])
                Q3 <- as.numeric(summary(currSamePeakMass)[5])
                LB <- Q1 - 1.5 *(Q3 - Q1)
                RB <- Q3 + 1.5 *(Q3 - Q1)
                if (currSamePeakMass[length(currSamePeakMass)] < LB || currSamePeakMass[length(currSamePeakMass)] > RB) break()
              }
            } else {
              abvector <- abs(tmpMZdata[[leftInd]] - currRefMz)
              NearInd <- which(abvector == min(abvector))[1]
              currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[leftInd]][NearInd])
              if(eicINT[leftInd] > eicINT[leftInd+1] & length(currSamePeakMass) > 5){
                Q1 <- as.numeric(summary(currSamePeakMass)[2])
                Q3 <- as.numeric(summary(currSamePeakMass)[5])
                LB <- Q1 - 1.5 *(Q3 - Q1)
                RB <- Q3 + 1.5 *(Q3 - Q1)
                if (currSamePeakMass[length(currSamePeakMass)] < LB || currSamePeakMass[length(currSamePeakMass)] > RB) break()
              }
            }
            leftInd <- leftInd-1
            if(leftInd <= 0) break()
          }
        }
        if(rightInd <= length(tmpMZdata)){
          while (length(tmpMZdata[[rightInd]]) > 0 & mean(tmpINTdata[[rightInd]]) >= cutOFF) {
            if (length(tmpMZdata[[rightInd]]) == 1){
              currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[rightInd]])
              if(eicINT[rightInd] > eicINT[rightInd-1] & length(currSamePeakMass) > 5){
                Q1 <- as.numeric(summary(currSamePeakMass)[2])
                Q3 <- as.numeric(summary(currSamePeakMass)[5])
                LB <- Q1 - 1.5 *(Q3 - Q1)
                RB <- Q3 + 1.5 *(Q3 - Q1)
                if (currSamePeakMass[length(currSamePeakMass)] < LB || currSamePeakMass[length(currSamePeakMass)] > RB) break()
              }
            } else {
              abvector <- abs(tmpMZdata[[rightInd]] - currRefMz)
              NearInd <- which(abvector == min(abvector))[1]
              currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[rightInd]][NearInd])
              if(eicINT[rightInd] > eicINT[rightInd-1] & length(currSamePeakMass) > 5){
                Q1 <- as.numeric(summary(currSamePeakMass)[2])
                Q3 <- as.numeric(summary(currSamePeakMass)[5])
                LB <- Q1 - 1.5 *(Q3 - Q1)
                RB <- Q3 + 1.5 *(Q3 - Q1)
                if (currSamePeakMass[length(currSamePeakMass)] < LB || currSamePeakMass[length(currSamePeakMass)] > RB) break()
              }
            }
            rightInd <- rightInd+1
            if(rightInd > length(tmpMZdata)) break()
          }
        }
        if(sum(is.na(currSamePeakMass)) > 0) next()
        if(length(currSamePeakMass) > 1 && length(currSamePeakMass) < 200){
          ppmCheck <- (massSDrange*sd(currSamePeakMass))/currRefMz * 1e6
          if(ppmCheck < ppmCut){
            daDiff <- massSDrange*sd(currSamePeakMass)
            ppmDiff <- (massSDrange*sd(currSamePeakMass))/currRefMz * 1e6
            currPeakWidth <- rtime[[rightInd - 1]] - rtime[[leftInd + 1]]
            currPeakScans <- rightInd - leftInd
            ppm2Ddist <- rbind(ppm2Ddist, c(currRefMz, rtime[[peakInd[z]]], ppmDiff))
            mzdiff2Ddist <- rbind(mzdiff2Ddist, c(currRefMz, rtime[[peakInd[z]]], daDiff))
            peakWidthALL <- c(peakWidthALL, currPeakWidth)
            peakScansALL <- c(peakScansALL, currPeakScans)
            if (cutOFF > 0){
              snALL <- c(snALL, (eicINT[currPeakInd]-blk)/sd)
            } else {
              snALL <- c(snALL, eicINT[currPeakInd])
            }
            sALL <- c(sALL, eicINT[currPeakInd])
            if (z == 1) {
              if(is.na(peakInd[z+1])){
                ShiftTable <- rbind(ShiftTable, c(currRefMz, rtime[[currPeakInd]], q))
              }
              if(!is.na(peakInd[z+1])) {
                if(rtime[[currPeakInd]] < (rtime[[peakInd[z+1]]] - 300)) {
                  ShiftTable <- rbind(ShiftTable, c(currRefMz, rtime[[currPeakInd]], q))
                }
              }
            }
            if (z > 1) {
              if(is.na(peakInd[z+1]) & rtime[[currPeakInd]] > (rtime[[peakInd[z-1]]] + 300)){
                ShiftTable <- rbind(ShiftTable, c(currRefMz, rtime[[currPeakInd]], q))
              }
              if(!is.na(peakInd[z+1])) {
                if(rtime[[currPeakInd]] < (rtime[[peakInd[z+1]]] - 300)) {
                  ShiftTable <- rbind(ShiftTable, c(currRefMz, rtime[[currPeakInd]], q))
                }
              }
            }
          }
        }
      }
      return(list(snALL = snALL, sALL = sALL, peakWidthALL = peakWidthALL, peakScansALL = peakScansALL,
                  ShiftTable = ShiftTable, ppm2Ddist = ppm2Ddist, mzdiff2Ddist = mzdiff2Ddist, noiseALL = noiseALL))
    }

    pb <- utils::txtProgressBar(max = (length(ROI) - 1), style = 3)
    if(thread == 1){
      tmpList <- lapply(1:(length(ROI) - 1), function(i) { # (length(ROI) - 1)
        utils::setTxtProgressBar(pb, i)
        calROI(i)
      })
    }else if(thread > 1){
      cl <- snow::makeCluster(thread)
      doSNOW::registerDoSNOW(cl)
      opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
                                                                   n))
      tmpList <- foreach::`%dopar%`(foreach::foreach(i = 1:(length(ROI) - 1),
                                                     .options.snow = opts), {
                                                       calROI(i)
                                                     })
      snow::stopCluster(cl)
      gc()
    }else{
      stop("thread is wrong!")
    }
    snALLList <- lapply(tmpList, function(x) {
      purrr::pluck(x, "snALL")
    })
    snALL <- purrr::list_c(snALLList)
    sALLList <- lapply(tmpList, function(x) {
      purrr::pluck(x, "sALL")
    })
    sALL <- purrr::list_c(sALLList)
    peakWidthALLList <- lapply(tmpList, function(x) {
      purrr::pluck(x, "peakWidthALL")
    })
    peakWidthALL <- purrr::list_c(peakWidthALLList)
    peakScansALLList <- lapply(tmpList, function(x) {
      purrr::pluck(x, "peakScansALL")
    })
    peakScansALL <- purrr::list_c(peakScansALLList)
    ShiftTableList <- lapply(tmpList, function(x) {
      purrr::pluck(x, "ShiftTable")
    })
    ShiftTable <- purrr::list_rbind(ShiftTableList)
    ppm2DdistList <- lapply(tmpList, function(x) {
      purrr::pluck(x, "ppm2Ddist")
    })
    ppm2Ddist <- purrr::list_rbind(ppm2DdistList)
    mzdiff2DdistList <- lapply(tmpList, function(x) {
      purrr::pluck(x, "mzdiff2Ddist")
    })
    mzdiff2Ddist <- purrr::list_rbind(mzdiff2DdistList)
    noiseALLList <- lapply(tmpList, function(x) {
      purrr::pluck(x, "noiseALL")
    })
    noiseALL <- purrr::list_c(noiseALLList)
    mzDiff2D <- rbind(mzDiff2D, mzdiff2Ddist)
    ppm2D <- rbind(ppm2D, ppm2Ddist)
    peakWidth <- c(peakWidth, peakWidthALL)
    peakScans <- c(peakScans, peakScansALL)
    SNRatio <- c(SNRatio, snALL)
    peakHeight <- c(peakHeight, sALL)
    Shift <- rbind(Shift, ShiftTable)
    noiselevel <- c(noiselevel, noiseALL)
  }

  # Estimation of the instrumental shift parameters
  if (length(file_path) > 1) {
    Shift <- Shift[complete.cases(Shift),]
    Shiftlist <- split(Shift, Shift[,3])
    if (length(Shiftlist) < 3){
      AlignmentT1 <- data.frame(matrix(nrow = 0, ncol = (ncol(Shiftlist[[1]]) + ncol(Shiftlist[[2]]))))
      for(i in 1:nrow(Shiftlist[[1]])){
        mass.lower.limit <- as.numeric(Shiftlist[[1]]$mz[i]) - 0.015
        mass.upper.limit <- as.numeric(Shiftlist[[1]]$mz[i]) + 0.015
        rt.lower.limit <- as.numeric(Shiftlist[[1]]$rt[i]) - 30
        rt.upper.limit <- as.numeric(Shiftlist[[1]]$rt[i]) + 30
        temp <- Shiftlist[[2]][(as.numeric(Shiftlist[[2]]$mz) >= as.numeric(mass.lower.limit) &
                                  as.numeric(Shiftlist[[2]]$mz) <= as.numeric(mass.upper.limit) &
                                  as.numeric(Shiftlist[[2]]$rt) >= as.numeric(rt.lower.limit) &
                                  as.numeric(Shiftlist[[2]]$rt) <= as.numeric(rt.upper.limit)),]
        temp <- temp[complete.cases(temp),]
        if(nrow(temp) > 0) {
          AlignmentT1 <- rbind(AlignmentT1, c(Shiftlist[[1]][i,], temp[1,]))
        }
      }
      AlignmentT1 <- AlignmentT1[complete.cases(AlignmentT1),]
      colnames(AlignmentT1) <- c(colnames(Shiftlist[[1]]), colnames(Shiftlist[[2]]))
    } else {
      AlignmentT1 <- data.frame(matrix(nrow = 0, ncol = (ncol(Shiftlist[[1]]) + ncol(Shiftlist[[2]]))))
      for(i in 1:nrow(Shiftlist[[1]])){
        mass.lower.limit <- as.numeric(Shiftlist[[1]]$mz[i]) - 0.015
        mass.upper.limit <- as.numeric(Shiftlist[[1]]$mz[i]) + 0.015
        rt.lower.limit <- as.numeric(Shiftlist[[1]]$rt[i]) - 30
        rt.upper.limit <- as.numeric(Shiftlist[[1]]$rt[i]) + 30
        temp <- Shiftlist[[2]][(as.numeric(Shiftlist[[2]]$mz) >= as.numeric(mass.lower.limit) &
                                  as.numeric(Shiftlist[[2]]$mz) <= as.numeric(mass.upper.limit) &
                                  as.numeric(Shiftlist[[2]]$rt) >= as.numeric(rt.lower.limit) &
                                  as.numeric(Shiftlist[[2]]$rt) <= as.numeric(rt.upper.limit)),]
        temp <- temp[complete.cases(temp),]
        if(nrow(temp) > 0) {
          AlignmentT1 <- rbind(AlignmentT1, c(Shiftlist[[1]][i,], temp[1,]))
        }
      }
      AlignmentT1 <- AlignmentT1[complete.cases(AlignmentT1),]
      colnames(AlignmentT1) <- c(colnames(Shiftlist[[1]]), colnames(Shiftlist[[2]]))

      for(k in 3:length(Shiftlist)){
        output <- data.frame(matrix(nrow = 0, ncol = ncol(AlignmentT1) + ncol(Shiftlist[[k]])))
        for(j in 1:nrow(AlignmentT1)){
          mass.lower.limit <- as.numeric(AlignmentT1[j,1]) - 0.015
          mass.upper.limit <- as.numeric(AlignmentT1[j,1]) + 0.015
          rt.lower.limit <- as.numeric(AlignmentT1[j,2]) - 30
          rt.upper.limit <- as.numeric(AlignmentT1[j,2]) + 30
          temp <- Shiftlist[[k]][(as.numeric(Shiftlist[[k]]$mz) >= as.numeric(mass.lower.limit) &
                                    as.numeric(Shiftlist[[k]]$mz) <= as.numeric(mass.upper.limit) &
                                    as.numeric(Shiftlist[[k]]$rt) >= as.numeric(rt.lower.limit) &
                                    as.numeric(Shiftlist[[k]]$rt) <= as.numeric(rt.upper.limit)),]
          temp <- temp[complete.cases(temp),]
          if(nrow(temp) > 0) {
            tmprow <- cbind(AlignmentT1[j,], temp[1,])
            colnames(tmprow) <- colnames(output)
            output <- rbind(output, tmprow)
          }
        }
        colnames(output) <- c(colnames(AlignmentT1), colnames(Shiftlist[[k]]))
        AlignmentT1 <- data.frame(matrix(nrow = 0, ncol = ncol(output)))
        AlignmentT1 <- output[complete.cases(output),]
      }
    }
    massShift <- as.data.frame(matrix(ncol = length(Shiftlist), nrow = 1))
    rtShift <- as.data.frame(matrix(ncol = length(Shiftlist), nrow = 1))
    massShift <- AlignmentT1[,seq(1, 3*length(Shiftlist)-2 ,3)]
    for (i in 1:nrow(massShift)){
      massShiftALL[i] <- max(massShift[i,]) - min(massShift[i,])
    }
    rtShift <- AlignmentT1[,seq(2, 3*length(Shiftlist)-1 ,3)]
    for (i in 1:nrow(rtShift)){
      rtShiftALL[i] <- max(rtShift[i,]) - min(rtShift[i,])
    }
  }
  noiselevel <- noiselevel[!is.na(noiselevel)]
  noiselevel <- noiselevel[order(noiselevel)]
  noiselevel <- noiselevel[1:round(length(noiselevel)*0.97)]
  ppm2D <- ppm2D[complete.cases(ppm2D),]
  ppm2D <- ppm2D[order(ppm2D[,3]),]
  mzDiff2D <- mzDiff2D[complete.cases(mzDiff2D),]
  mzDiff2D <- mzDiff2D[order(mzDiff2D[,3]),]
  mzDiff2D <- mzDiff2D[1:round(nrow(mzDiff2D)*0.97),]
  peakWidth <- peakWidth[!is.na(peakWidth)]
  peakWidth <- peakWidth[order(peakWidth)]
  peakWidth <- peakWidth[1:round(length(peakWidth)*0.97)]
  peakScans <- peakScans[!is.na(peakScans)]
  peakScans <- peakScans[order(peakScans)]
  peakScans <- peakScans[1:round(length(peakScans)*0.97)]
  SNRatio <- SNRatio[!is.na(SNRatio)]
  SNRatio <- SNRatio[order(-SNRatio)]
  SNRatio <- SNRatio[1:round(length(SNRatio)*0.97)]
  peakHeight <- peakHeight[!is.na(peakHeight)]
  peakHeight <- peakHeight[order(-peakHeight)]
  peakHeight <- peakHeight[1:round(length(peakHeight)*0.97)]
  if (length(file_path) > 1) {
    massShiftALL <- massShiftALL[!is.na(massShiftALL)]
    massShiftALL <- massShiftALL[order(massShiftALL)]
    massShiftALL <- massShiftALL[1:round(length(massShiftALL)*0.97)]
    rtShiftALL <- rtShiftALL[!is.na(rtShiftALL)]
    rtShiftALL <- rtShiftALL[order(rtShiftALL)]
    rtShiftALL <- rtShiftALL[1:round(length(rtShiftALL)*0.97)]
    maxmassshift <- max(massShiftALL)
    maxrtshift <- max(rtShiftALL)
  }
  maxppm <- max(ppm2D[,3])
  maxppm <- ceiling(maxppm)
  maxmzdiff <- max(mzDiff2D[,3])
  maxmzdiff <- (ceiling(maxmzdiff*100))/100
  minnoise <- min(noiselevel)
  minnoise <- floor(minnoise)
  W <- mean(peakWidth, trim=0.05, na.rm = TRUE)
  H <- mean(peakHeight, trim=0.05, na.rm = TRUE)
  ratio <- H/W
  minpeakwidth <- min(peakWidth)
  maxpeakwidth <- max(peakWidth)
  if (maxpeakwidth > 35 & ratio > 515) {
    minpeakwidth <- 0
    maxpeakwidth <- (ceiling(maxpeakwidth)+7)/2
  } else {
    minpeakwidth <- (ceiling(minpeakwidth))+4
    maxpeakwidth <- ceiling(maxpeakwidth)+5
  }
  minpeakscan <- min(peakScans)
  minpeakscan <- floor(minpeakscan)
  maxpeakscan <- max(peakScans)
  maxpeakscan <- floor(maxpeakscan)
  minSN <- min(SNRatio)
  if (minSN < 3){
    minSN <- 3
  }
  minpeakheight <- min(peakHeight)
  minpeakheight <- floor(minpeakheight)
  P <- c("ppm", "minimum peakwidth", "maximum peakwidth", "signal/noise threshold", "mzdiff(Please refer to the user manual for more detailed explanation)", "Integration method",
         "prefilter peaks", "prefilter intensity", "noise filter", "bw", "minfrac", "mzwid", "minsamp", "max")
  V <- c(maxppm, minpeakwidth, maxpeakwidth, minSN, -0.01, 1, minpeakscan, minnoise, minnoise, 5, 0.5, maxmassshift, 1, 100)
  parameters <- round(V, 3)
  names(parameters) <- P
  return(parameters)
}
