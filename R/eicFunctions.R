#' This function isolates extracted ion chromatograms for a given list of m/z ratios
#' 
#' @param msData A OnDiskExp from which EICs shall be extracted
#' @param msPeaks List of peaks for which EICs shall be extracted
#' @param rtWindow Window in seconds to enlarge EICs around peak boundaries
#' @param msLevel Indicates from which MS level EICs shall be extracted
#' @param adjustedRtime TRUE or FALSE if adjusted retention times shall be used
#' @param fileIndex index of file for which EICs shall be extracted
#' 
#' @return returns a data frame containing the m/z value of the candidates and their EICs.
getEics <- function(msData, msPeaks, rtwindow, msLevel, adjustedRtime = TRUE, fileIndex = 1) {
  
  #create empty data frame for EICs
  eicDf <- data.frame()
  
  if(!is.null(nrow(msPeaks))) {
    #iterate through candidates
    for(i in 1:nrow(msPeaks)) {
      
      #create extracted ion chromatogram
      #extract a bit larger than needed to have enough data points
      bpis <- chromatogram(msData,
                           msLevel = msLevel,
                           mz=c(msPeaks[i,2], msPeaks[i,3]),
                           rt=c(msPeaks[i,5] - rtwindow, msPeaks[i,6] + rtwindow),
                           adjustedRtime = adjustedRtime)
      
      #add to data frame
      eicDf <- rbind.data.frame(eicDf, cbind.data.frame(mz = msPeaks[i,1],
                                                        rtime = bpis[[fileIndex]]@rtime,
                                                        intensity = bpis[[fileIndex]]@intensity))
      
    }
  } else {
    
    #create extracted ion chromatogram
    #extract a bit larger than needed to have enough data points
    bpis <- chromatogram(msData,
                         msLevel = msLevel,
                         mz=c(msPeaks[2], msPeaks[3]),
                         rt=c(msPeaks[5] - rtwindow, msPeaks[6] + rtwindow),
                         adjustedRtime = adjustedRtime)
    
    #add to data frame
    eicDf <- rbind.data.frame(eicDf, cbind.data.frame(mz = msPeaks[1],
                                                      rtime = bbpis[[fileIndex]]@rtime,
                                                      intensity = bpis[[fileIndex]]@intensity))
    
    
  }
  
  #return data frame
  return(eicDf)
  
}

#' This function aligns the MS1 and MS2 traces to overcome the lag between them based on the different dwell times
#' 
#' @param ms1df Data frame with the EIC for the precursor
#' @param ms2df data frome with the EICs for all fragment candidates
#' 
#' @return returns a data frame containing the m/z value of the candidates and their EICs corrected according to the MS1 EICs
alignMsLevel <- function(ms1df, ms2df) {
  
  #create data frame to store corrected EICs
  ms2dfCor <- data.frame()
  
  #iterate through all MS2 candidate EICs
  for(i in unique(ms2df$mz)) {
    
    #isolate individual EIC to work with
    candEIC <- ms2df[which(ms2df$mz == i),]
    
    #use linear approximation to get values at RTs of MS1
    candEICCor <- approx(candEIC$rtime, candEIC$intensity, ms1df$rtime)
    rtimeCor <- candEICCor$x
    intensityCor <- candEICCor$y
    
    #add to result data frame
    ms2dfCor <- rbind.data.frame(ms2dfCor, cbind.data.frame(mz = i,
                                                            rtimeCor = rtimeCor,
                                                            intensityCor = intensityCor))
  }
  
  #return new data frame
  return(ms2dfCor)
  
}

#' This function correlates the MS1 and MS2 traces and returns a data frame with the pearson correlation coefficients
#' 
#' @param ms1df Data frame with the EIC for the precursor
#' @param ms2dfCor data frome with the EICs for all fragment candidates
#' 
#' @return a data frame containing m/z ratios and correlation coefficients
correlateMsLevel <- function(ms1df, ms2dfCor) {
  
  #create data frame for correlation values
  corDf <- data.frame()
  
  #iterate through all MS2 candidate EICs
  for(i in unique(ms2dfCor$mz)) {
    
    #isolate individual corrected EIC to work with
    corEIC <- ms2dfCor[which(ms2dfCor$mz == i),]
    
    #perform correlation with Ms1 level
    corValue <- cor(ms1df$intensity, corEIC$intensityCor)
    
    #add results to data frame
    corDf <- rbind.data.frame(corDf, cbind.data.frame(mz = i,
                                                      corValue = corValue))
    
  }
  
  #return new data frame
  return(corDf)
  
}