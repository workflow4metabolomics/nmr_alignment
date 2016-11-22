################################################################################################
# NMR SPECTRA ALIGNEMENT USING THE CLUPA ALGORITHM (Vu et al. 2013)                            #
# Included in the speaq package                                                                #
# User : Galaxy                                                                                #
# Original data : --                                                                           #
# Starting date : 05-10-2016                                                                   #
# Version 1 : --                                                                               #
#                                                                                              #
# Input files : datra.frame containing spectra to align                                        #
################################################################################################


################################################################################################
# Bruker files reading                                                                         #
# Input parameters                                                                             #
#   - directory: name of your folder containing all experimental samples as sub-directories    #
#   - leftBorder: upper boundary: values greater than this value are not used;                 #
#                 Default value: 10.0 ppm                                                      #
#   - rightBorder: lower boundary: values lower than this value are not used;                  #
#                 Default value: 10.0 ppm                                                      #
#   - exclusionZones: spectral regions to exclude (water, solvent or contaminant resonance);   #
#                   boolean                                                                    #
#   - exclusionZonesBorders: upper and lower boudaries of exclusions zones                     #
# Output parameters                                                                            #
#   - truncatedSpectrum_matrice: n x p datamatrix                                              #
################################################################################################
NmrRead <- function(directory, leftBorder= 10.0, rightBorder= 0.5, exclusionZones=FALSE, exclusionZonesBorders=NULL)
{
  ## Option
  ##---------------
  strAsFacL <- options()$stringsAsFactors
  options(stingsAsFactors=FALSE)
  options(warn=-1)
  
  
  ## Constants
  ##---------------
  topEnvC <- environment()
  flgC <- "\n"
  
  ## Log file (in case of integration into Galaxy)
  ##----------------------------------------------
  #  if(!is.null(savLog.txtC))
  #    sink(savLog.txtC, append=TRUE)
  
  ## Functions definition
  ##---------------------  
  ## RAW BRUKER FILE READING FUNCTION
  NmRBrucker_read <- function(DataDir,SampleSpectrum)
  {
    
    bruker.get_param <- function (ACQ,paramStr)
    {
      regexpStr <- paste("^...",paramStr,"=",sep="")
      as.numeric(gsub("^[^=]+= ","" ,ACQ[which(simplify2array(regexpr(regexpStr,ACQ))>0)]))
    }
    
    ACQFILE <- "acqus"
    SPECFILE <- paste(DataDir,"/1r",sep="")
    PROCFILE <- paste(DataDir,"/procs",sep="")
    
    ACQ <- readLines(ACQFILE)
    TD      <- bruker.get_param(ACQ,"TD")
    SW      <- bruker.get_param(ACQ,"SW")
    SWH     <- bruker.get_param(ACQ,"SW_h")
    DTYPA   <- bruker.get_param(ACQ,"DTYPA")
    BYTORDA <- bruker.get_param(ACQ,"BYTORDA")
    #ENDIAN=ifelse( BYTORDA==0, "little", "big")
    ENDIAN <- "little"
    SIZE=ifelse( DTYPA==0, 4, 8)
    
    PROC <- readLines(PROCFILE)
    OFFSET <- bruker.get_param(PROC,"OFFSET")
    SI <- bruker.get_param(PROC,"SI")
    
    to.read=file(SPECFILE,"rb")
    maxTDSI=max(TD,SI)
    #  signal<-rev(readBin(to.read, what="int",size=SIZE, n=TD, signed=TRUE, endian=ENDIAN))
    signal<-rev(readBin(to.read, what="int",size=SIZE, n=maxTDSI, signed=TRUE, endian=ENDIAN))
    close(to.read)
    
    td <- length(signal)
    
    #  dppm <- SW/(TD-1)
    dppm <- SW/(td-1)
    pmax <- OFFSET
    pmin <- OFFSET - SW
    ppmseq <- seq(from=pmin, to=pmax, by=dppm)
    signal <- 100*signal/max(signal)
    
    SampleSpectrum <- cbind(ppmseq,signal)
    return(SampleSpectrum)
  }
  
  
  # File names
  FileNames <- list.files(directory)
  n <- length(FileNames)

  # Reading and Bucketing
  directory <- paste(directory,"/",sep="")
  
  i <- 1
  while (i <= n)
  {
    # File reading
    SampleDir <- paste(directory,FileNames[i],"/1/",sep="")
    setwd(SampleDir)
    DataDir <- "pdata/1"
    rawSpectrum <- NmRBrucker_read(DataDir,rawSpectrum)
    
    orderedSpectrum <- rawSpectrum[order(rawSpectrum[,1],decreasing=T), ]
    
    # Removal of chemical shifts > leftBorder or < rightBorder boundaries
    truncatedSpectrum <- orderedSpectrum[orderedSpectrum[,1] < leftBorder & orderedSpectrum[,1] > rightBorder, ]
    truncatedSpectrum[,1] <- round(truncatedSpectrum[,1],4)

    # Exclusion zones
      if (!is.null(exclusionZonesBorders))
      {
        truncatedSpectrum[truncatedSpectrum[,1] < exclusionZonesBorders[[1]][1] & truncatedSpectrum[,1] > exclusionZonesBorders[[1]][2],2] <- 0
        
        if (length(exclusionZonesBorders) > 1)
          for (k in 2:length(exclusionZonesBorders))
            truncatedSpectrum[truncatedSpectrum[,1] < exclusionZonesBorders[[k]][1] & truncatedSpectrum[,1] > exclusionZonesBorders[[k]][2],2] <- 0
      }

    # spectrum Concatenation
    if (i==1)
      truncatedSpectrum_matrice <- truncatedSpectrum
    if (i > 1)
      truncatedSpectrum_matrice <- cbind(truncatedSpectrum_matrice,truncatedSpectrum[,2])
    colnames(truncatedSpectrum_matrice)[i+1] <- FileNames[i]
    
    # Next sample
    rm(spectrum.bucket)
    i <- i +1
  }
  
  identifiants <- gsub("([- , * { } | \\[ ])","_",colnames(truncatedSpectrum_matrice)[-1])
  colnames(truncatedSpectrum_matrice) <- c(colnames(truncatedSpectrum_matrice)[1],identifiants)
  
  # Directory
  setwd(directory)
  return(truncatedSpectrum_matrice)
}



################################################################################################
# Peak detection for spectra                                                                   #
# Input parameters                                                                             #
#   - X: spectral dataset in matrix format in which each row contains a single sample          #
#   - nDivRange: size of a single small segment after division of spectra                      #
#                 Default value: 64                                                            #
#   - baselineThresh: removal of all the peaks with intensity lower than this threshold        #
#                 Default value: 50000                                                         #
# Output parameters                                                                            #
#   - peak lists of the spectra                                                                #
################################################################################################
detectSpecPeaks <- function (X, nDivRange, scales=seq(1, 16, 2), baselineThresh,SNR.Th=-1, verbose) 
{
  nFea = ncol(X)
  nSamp = nrow(X)
  noiseEsp = 0.005
  if (SNR.Th < 0) 
    SNR.Th = max(scales) * 0.05
  pList = NULL
  for (i in 1:nSamp) {
    myPeakRes = NULL
    mySpec = X[i, ]
    for (k in 1:length(nDivRange)) {
      divR = nDivRange[k]
      for (j in 1:(trunc(nFea/divR) - 3)) {
        startR = (j - 1) * divR + 1
        if (startR >= nFea) 
          startR = nFea
        endR = (j + 3) * divR
        if (endR > nFea) 
          endR = nFea
        xRange = mySpec[startR:endR]
        xMean = mean(xRange)
        xMedian = median(xRange)
        if ((xMean == xMedian) || abs(xMean - xMedian)/((xMean + 
                                                         xMedian) * 2) < noiseEsp) {
          next
        }
        else {
          peakInfo = peakDetectionCWT(mySpec[startR:endR], 
                                      scales = scales, SNR.Th = SNR.Th)
          majorPeakInfo = peakInfo$majorPeakInfo
          if (length(majorPeakInfo$peakIndex) > 0) {
            myPeakRes = c(myPeakRes, majorPeakInfo$peakIndex + 
                            startR - 1)
          }
        }
      }
    }
    pList[i] = list(myPeakRes)
    pList[[i]] = sort(unique(pList[[i]]))
    pList[[i]] = pList[[i]][which(X[i, pList[[i]]] > baselineThresh)]
    pList[[i]] = sort(pList[[i]])
    if (verbose) 
      cat("\n Spectrum ", i, " has ", length(pList[[i]]), 
          " peaks")
  }
  return(pList)
}


################################################################################################
# Heuristical detection of a reference spectrum                                                #
# Input parameters                                                                             #
#   - peakList: peak lists of the spectra                                                      #
# Output parameters: list containing                                                           #
#   - refInd: index of the reference spectrum found by the algorithm                           #
#   - orderSpec: sorted array of the spectra by their goodness values                          #
################################################################################################
findRef <- function (peakList) 
{
  disS = matrix(data = NA, ncol = length(peakList), nrow = length(peakList))
  sumDis = double(length(peakList))
  for (refInd in 1:length(peakList)) {
    for (tarInd in 1:length(peakList)) if (refInd != tarInd) {
      disS[refInd, tarInd] = 0
      for (i in 1:length(peakList[[tarInd]])) disS[refInd, 
                                                   tarInd] = disS[refInd, tarInd] + min(abs(peakList[[tarInd]][i] - 
                                                                                              peakList[[refInd]]))
    }
  }
  for (refInd in 1:length(peakList)) {
    disS[refInd, refInd] = 0
    sumDis[refInd] = sum(disS[refInd, ])
  }
  orderSumdis = order(sumDis)
  
  return(list(refInd = orderSumdis[1], orderSpec = orderSumdis))
}


################################################################################################
# CluPA alignment for multiple spectra                                                         #
# Input parameters                                                                             #
#   - X: spectral dataset in the matrix format in which each row contains a single sample      #
#   - peakList: peak lists of the spectra                                                      #
#   - refInd: index of the reference spectrum                                                  # 
#   - maxShift:  maximum number of the points for a shift step                                 #
# Output parameters                                                                            #
#   - aligned spectra: dataframe?                                                              #
################################################################################################
dohCluster <- function (X, peakList, refInd = 0, maxShift = maxshift, acceptLostPeak = TRUE, verbose) 
{
  Y = X
  peakListNew = peakList
  if (verbose) 
    startTime = proc.time()
  refSpec = Y[refInd, ]
  for (tarInd in 1:nrow(X)) if (tarInd != refInd) {
    if (verbose) 
      cat("\n aligning spectrum ", tarInd)
    targetSpec = Y[tarInd, ]
    myPeakList = c(peakList[[refInd]], peakList[[tarInd]])
    myPeakLabel = double(length(myPeakList))
    for (i in 1:length(peakList[[refInd]])) myPeakLabel[i] = 1
    startP = 1
    endP = length(targetSpec)
    res = hClustAlign(refSpec, targetSpec, myPeakList, myPeakLabel, 
                      startP, endP, maxShift = 50, acceptLostPeak = TRUE)
    Y[tarInd, ] = res$tarSpec
    peakListNew[[tarInd]] = res$PeakList[(length(peakList[[refInd]]) + 
                                            1):length(myPeakList)]
  }
  peakList = peakListNew
  if (verbose) {
    endTime = proc.time()
    cat("\n Alignment time: ", (endTime[3] - startTime[3])/60, 
        " minutes")
  }
  return(Y)
}


################################################################################################
# Spectra display                                                                              #
# Input parameters                                                                             #
#   - X: datamatrix                                                                            #
#   - ...                                                                                      #
# Output parameters                                                                            #
#   - truncatedSpectrum_matrice: n x p datamatrix                                              #
################################################################################################
drawSpec <- function (X, startP = -1, endP = -1, groupLabel = NULL, useLog = -1, highBound = -1, 
                      lowBound = -1, xlab = NULL, ylab = NULL, main = NULL, nAxisPos = 4, offside = 0) 
{
  groupLabel_name = groupLabel
  X = as.data.frame(X)
  #  colnames(X) = c(1:ncol(X))
  X = as.matrix(X)
  if (highBound != -1) {
    for (i in 1:nrow(X)) {
      myIndex = which(X[i, ] > highBound)
      X[i, myIndex] = highBound
    }
  }
  if (lowBound != -1) {
    for (i in 1:nrow(X)) {
      myIndex = which(X[i, ] < lowBound)
      X[i, myIndex] = lowBound
    }
  }
  if (is.null(groupLabel)) {
    groupLabel = c(1:nrow(X))
    groupLabel = as.factor(groupLabel)
  }
  else {
    levels(groupLabel) = c(1:length(levels(groupLabel)))
  }
  if (startP == -1) 
    startP = 1
  if (endP == -1) 
    endP = ncol(X)
  if (is.null(xlab)) {
    xlab = "index"
  }
  if (is.null(ylab)) {
    ylab = "intensity"
  }
  if (is.null(main)) {
    main = paste(" ", startP + offside, "-", endP + offside)
  }
  GraphRange <- c(startP:endP)
  yn <- X[, GraphRange]
  if (useLog != -1) 
    yn = log(yn)
  plot(yn[1, ], ylim = c(min(yn), max(yn)), type = "n", ylab = ylab, xlab = xlab, main = main, xaxt = "n")
  tempVal = trunc(length(GraphRange)/nAxisPos)
  xPos = c(0:nAxisPos) * tempVal
  #  axis(1, at = xPos, labels = xPos + startP + offside)
  axis(1, at = xPos, labels = colnames(X)[xPos + startP + offside])
  for (i in 1:length(levels(groupLabel))) {
    groupLabelIdx = which(groupLabel == levels(groupLabel)[i])
    color <- palette(rainbow(length(levels(groupLabel))))
    for (j in 1:length(groupLabelIdx)) {
      #      lines(yn[groupLabelIdx[j], ], col = as.integer(levels(groupLabel)[i]))
      lines(yn[groupLabelIdx[j], ], col = color[i])
    }
  }
  if (!is.null(groupLabel_name)) {
    legendPos = "topleft"
    legend(legendPos, levels(groupLabel_name), col = as.integer(levels(groupLabel)), text.col = "black", pch = c(19, 19), bg = "gray90")
  }
}


################################################################################################
# Spectra alignment                                                                            #
# Input parameters                                                                             #
#   - data: n x p datamatrix                                                                   #
#   - nDivRange: size of a single small segment after division of the whole spectrum           #
#                 Default value: 64                                                            #
#   - reference: number of the spectrum reference; if NULL, automatic detection                #
#                Default value: NULL                                                           #
#   - baselineThresh: removal of all the peaks with intensity lower than this threshold        #
#                 Default value: 50000                                                         #
# Output parameters                                                                            #
#   - Y: dataframe (?)                                                                         #
################################################################################################
cluPA.alignment <- function(data, reference=reference, nDivRange, scales = seq(1, 16, 2), baselineThresh,  SNR.Th = -1, maxshift=maxshift, 
                            verbose)
{
  ## Peak picking
  cat("\n detect peaks....")
  startTime <- proc.time()
  peakList <- detectSpecPeaks(X=data, nDivRange=nDivRange, scales=scales, baselineThresh=baselineThresh,  
                              SNR.Th = SNR.Th, verbose=verbose)
  endTime <- proc.time()
  cat("Peak detection time:",(endTime[3]-startTime[3])/60," minutes")

  ## Reference spectrum determination
  if (reference == 0)
  {
    cat("\n Find the spectrum reference...")
    resFindRef<- findRef(peakList)
    refInd <- resFindRef$refInd
    cat("\n Order of spectrum for reference \n")
    for (i in 1:length(resFindRef$orderSpec))
    {
      cat(paste(i, ":",resFindRef$orderSpec[i],sep=""), " ")
      if (i %% 10 == 0) 
        cat("\n")
    }
    cat("\n The reference is: ", refInd)
  }
  else
  {
    refInd=reference
  }
  ## Spectra alignment to the reference
  maxshift <- 50
  Y <- dohCluster(data, peakList=peakList, refInd=refInd, maxShift=maxShift, acceptLostPeak, verbose)

  ## Output  
  return(Y)
}


################################################################################################
# Principal function                                                                           #
# Input parameters                                                                             #
#   - directory: name of your folder containing all experimental samples as sub-directories    #
#   - leftBorder: upper boundary: values greater than this value are not used;                 #
#                 Default value: 10.0 ppm                                                      #
#   - rightBorder: lower boundary: values lower than this value are not used;                  #
#                 Default value: 10.0 ppm                                                      #
#   - exclusionZones: spectral regions to exclude (water, solvent or contaminant resonance);   #
#                   boolean                                                                    #
#   - exclusionZonesBorders: upper and lower boudaries of exclusions zones                     #
#   - verbose: printing out process information (boolean)                                      #
#   - reference: reference spectrum number, if exists                                          #
#                Default value: 0 (automatic detection                                         #
#   - nDivRange: size of a single small segment after division of the whole spectrum           #
#                 Default value: 64                                                            #
#   - baselineThresh: removal of all the peaks with intensity lower than this threshold        #
#                 Default value: 50000                                                         #
# Output parameters: list containing                                                           #
#   - data.read: n x p matrix                                                                  #
#   - data.aligned: n x p matrix                                                               #
################################################################################################
nmr.alignment <- function(directory, leftBorder= 10.0, rightBorder= 0.5, exclusionZones=FALSE, 
                          exclusionZonesBorders=NULL, reference=0, nDivRange=64, baselineThresh=50000, maxshift=50, verbose=FALSE)
{
  data.read <- NmrRead(directory=directory, leftBorder=leftBorder, rightBorder=rightBorder, exclusionZones=exclusionZones, exclusionZonesBorders=exclusionZonesBorders)
  rownames(data.read) <- data.read[,1]
  data.read <- data.read[,-1]
  data.read <- t(data.read)
  data.aligned <- cluPA.alignment(data=data.read, reference=reference, nDivRange=nDivRange, 
                                  baselineThresh=baselineThresh, maxshift=maxshift, verbose=verbose)

  return(list(data.read,data.aligned))
}

