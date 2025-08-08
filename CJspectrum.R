CJspectrum <- function(data, qMinimum, qMaximum, qStep) { 
  # data is be a column vector
  # qMinimum is minimum q that you want to test
  # qMaximum is maximum q that you want to test
  # qStep is q incrementing
  MinimumBinSize <- 4 # Set up minimum bin size
  LengthOfData <- length(data) # Determine length of data
  MaximumBinSize <- floor(LengthOfData / 4) # Set up maximum bin size
  AllBinSizes <- 2^c(2:floor(log(MaximumBinSize) / log(2))) # Set up the vector of bin sizes
  qVector <- seq(qMinimum, qMaximum, by = qStep) # Set up the q vector
  LengthOfqVector <- length(qVector) # Determine length of q vector
  FqOriginal <- mat.or.vec((qMaximum - qMinimum) / qStep + 1, 1) # Generate a 21-row vector for each q’s estimate of the Hausdorff dimension
  FCorrCoefOriginal <- mat.or.vec((qMaximum - qMinimum) / qStep + 1, 1) # Generate a 21-row vector for each q’s correlation coefficient relating negative Shannon entropy and logarithmic bin size
  AqOriginal <- mat.or.vec((qMaximum - qMinimum) / qStep + 1, 1) # Generate a 21-row vector for each q’s estimate of singularity strength
  ACorrCoefOriginal <- mat.or.vec((qMaximum - qMinimum) / qStep + 1, 1) # Generate a 21-row vector for each q’s correlation coefficient relating mass-weighted logarithmic proportion and logarithmic bin size
  FqNumer <- mat.or.vec(length(AllBinSizes), LengthOfqVector) # Create matrix to store values of Entropymass for different bin size and q
  AqNumer <- mat.or.vec(length(AllBinSizes), LengthOfqVector) # Create matrix to store values of sMWlBinProportion for different bin size and q
  for (Count in 1:length(AllBinSizes)) { # Set up the loop for bin size
    CurrentIncrement <- AllBinSizes[Count] # Set up current increment
    NumberOfIncrements <- floor(length(data) / CurrentIncrement) # Compute the number of increments
    BinnableData <- data[1:(NumberOfIncrements * CurrentIncrement)] # Cut out the binnable portion of the series
    DataMatrix <- matrix(BinnableData, nrow = CurrentIncrement, byrow = FALSE) # Reshape the data matrix
    MeasureInEachBin <- colMeans(DataMatrix) # Compute the column means
    TotalMeasureInBins <- sum(MeasureInEachBin) # Compute the total bin proportion across the entire series
    BinProportion <- MeasureInEachBin/TotalMeasureInBins # Compute bin proportions
    for (qCount in 1:LengthOfqVector) { # Set up the loop for q
      qParameter <- qVector[qCount] # Assign value of q
      qBinProportion <- BinProportion^qParameter # Compute qth-order bin proportion
      SumOfqBinProportion <- sum(qBinProportion) # Compute the total mass across the entire series
      BinMass <- (qBinProportion)/SumOfqBinProportion # Compute bin mass as bin’s proportion of total mass
      lBinMass <- log(BinMass) # Compute logarithmic bin mass
      EntropyBinMass <- BinMass * lBinMass # Compute negative Shannon entropy as mass-weighted logarithmic bin mass
      EntropyMass <- sum(EntropyBinMass) # Compute the total negative Shannon entropy across the entire series
      FqNumer[Count, qCount] <- EntropyMass # Assign to FqNumer
      lBinProportion <- log(BinProportion) # Compute logarithmic bin proportion
      MlBinProportion <- BinMass * lBinProportion # Compute mass-weighted logarithmic bin proportion
      SumOfMlBinProportion <- sum(MlBinProportion) # Compute the total mass-weighted logarithmic bin proportion across the entire series
      AqNumer[Count, qCount] <- SumOfMlBinProportion # Assign to AqNumer
    }
  }
  for (qCount in 1:LengthOfqVector) {
    FqNumerq <- FqNumer[, qCount] # Extract values of FqNumerq for specific q
    AqNumerq <- AqNumer[, qCount] # Extract values of AqNumerq for specific q
    FqNumerq <- FqNumerq[is.finite(FqNumerq)] # Remove nonfinite values from FqNumerq
    AqNumerq <- AqNumerq[is.finite(AqNumerq)] # Remove nonfinite values from AqNumerq
    FqNumerq <- FqNumerq[!is.na(FqNumerq)] # Remove NaN values from FqNumerq
    AqNumerq <- AqNumerq[!is.na(AqNumerq)] # Remove NaN values from AqNumerq
    LengthDifferenceF <- length(FqNumer[, qCount]) - length(FqNumerq) # Find the number of values removed from FqNumerq
    LengthDifferenceA <- length(AqNumer[, qCount]) - length(AqNumerq) # Find the number of values removed from AqNumerq
    AllBinSizesF <- AllBinSizes[(1 + LengthDifferenceF):length(AllBinSizes)] # Select bin sizes with valid values of FqNumerq
    AllBinSizesA <- AllBinSizes[(1 + LengthDifferenceA):length(AllBinSizes)] # Select bin sizes with valid values of AqNumerq
    lBinSizeF <- log(AllBinSizesF) # Compute logarithmic bin size for SlopeOfF
    lBinSizeA <- log(AllBinSizesA) # Compute logarithmic bin size for SlopeOfA
    SlopeOfF <- lsfit(lBinSizeF,FqNumerq) # Compute the slope of negative Shannon entropy vs. logarithmic bin size
    SlopeOfA <- lsfit(lBinSizeA,AqNumerq) # Compute the slope of mass-weighted logarithmic bin proportion vs. logarithmic bin size
    FqOriginal[qCount] <- SlopeOfF$coefficients[2] # Store slopes for negative Shannon entropy
    AqOriginal[qCount] <- SlopeOfA$coefficients[2] # Store slopes for mass-weighted logarithmic bin proportion vs. logarithmic bin size
    FCorrCoefOriginal[qCount] <- cor(lBinSizeF,FqNumerq) # Find correlation coefficient of negative Shannon entropy vs. logarithmic bin size 
    ACorrCoefOriginal[qCount] <- cor(lBinSizeA,AqNumerq) # Find correlation coefficient of mass-weighted logarithmic bin proportion vs. logarithmic bin size
  }
  CriticalAlpha <- 0 # Initiate variable to store the points on the spectrum that meet the correlation-coefficient criterion
  for (Point in 1:length(AqOriginal)) { # Set up the loop
    if (FCorrCoefOriginal[Point] > 0.995) { # Apply threshold criterion on the Hausdorff dimension
      if(ACorrCoefOriginal[Point] > 0.995) { # Apply threshold criterion on singularity strength
        CriticalAlpha <- cbind(CriticalAlpha, AqOriginal[Point]) # Extract values of AqOriginal that meet the correlation-coefficient criterion
      }
    }
  }
  RangeA <- range(CriticalAlpha[2:length(CriticalAlpha)]) # First value of is CriticalAlpha is 0 because of cbind() building of CriticalAlpha
  RangeA <- RangeA[2] - RangeA[1] # Compute multifractal spectrum width
  return(RangeA)
}