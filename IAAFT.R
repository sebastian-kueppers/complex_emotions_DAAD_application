IAAFT <- function(data, NumberOfSurrogates, MaximumNumberOfIteractions) {
  # data is time series
  # NumberOfSurrogates is the number of surrogates
  # MaximumNumberOfIteractions is the maximum number of iterations, usually 1000
  fftdata <- fft(data) # Compute fast Fourier transformation of the original series
  Amplitude <- abs(fftdata) # Remove sign of the computed fast Fourier transformation of the original series
  Surrogates <- matrix(rep(data, NumberOfSurrogates), ncol = NumberOfSurrogates) # Replicate for each to-be-generated surrogate
  Surrogates <- apply(Surrogates, 2, sample) # Randomize surrogates within samples
  SortOutput <- sort(data, index.return = TRUE) # Sort the original series in ascending order, storing the index of the sorted values from the original series in sortOutput$ix
  fftSurrogates <- mvfft(Surrogates) # Perform fast Fourier transformation concurrently on each colum of Surrogates
  fftSurrogates.re <- Re(fftSurrogates) # Save the real component of fftSurrogates
  fftSurrogates.im <- Im(fftSurrogates) # Save the imaginary component of fftSurrogates
  fftRatio <- fftSurrogates.im / fftSurrogates.re # Compute the ratio of the imaginary and real Fourier components
  PhaseSpectra <- atan(fftRatio) # Computes phase spectrum as the arctangent of fftRatio
  
  # print(PhaseSpectra)
  
  nn <- 1
  sortOutput <- sort(data, index.return = TRUE) # Sort the original series in ascending order, storing the index of the sorted values from the original series in sortOutput$ix
  SortedValues <- sortOutput$x # Store the ascending-sorted values of the series
  IndexOfSortedValues <- sortOutput$ix # Store the indices of ascending-sorted values of the series
  while (nn <= MaximumNumberOfIteractions) {
    Phase.i <- PhaseSpectra * 1i # Multiply the phase spectrum by the imaginary number i to produce an exponent for the Euler formula for the Fourier transform (note that this phase spectrum was computed from a randomization of the original series, and hence, it refers to the different sequence than that of the original series)
    # print(Phase.i)
    ExponentiatedPhase <- exp(Phase.i) # Exponentiate Phase.i
    
    Surrogates <- Amplitude * ExponentiatedPhase # Replace all columns of Surrogates with the Euler formula for the Fourier series, i.e., with the amplitude spectrum multiplied by the exponentiation of an imaginary phase spectrum
    
    # AmplitudeMatrix <- matrix(Amplitude, nrow = length(Amplitude), ncol = NumberOfSurrogates)
    # Surrogates <- AmplitudeMatrix * ExponentiatedPhase
    # print(Surrogates)
    
    ifftSurrogates <- mvfft(Surrogates, inverse = TRUE) # Compute the inverse-Fourier transform of these Fourier series
    Surrogates <- Re(ifftSurrogates) # Store the real components of the inverse-Fourier transform, i.e., of the surrogate series obtained in the current iteration
    lst <- apply(Surrogates, 2, sort, index.return = TRUE) # Sort each column of Surrogates and store the rank-order information
    SortStructure <- do.call(rbind, Map(cbind, colInd = seq_along(lst), lapply(lst, function(x) do.call(cbind, x))))
    Indices <- matrix(SortStructure[, 3], ncol = NumberOfSurrogates) # Store the indices of the values in rank order in a reshaped matrix
    lst <- apply(Indices, 2, sort, index.return = TRUE) # Sort Indices by column
    SortStructure <- do.call(rbind, Map(cbind, colInd = seq_along(lst), lapply(lst, function(x) do.call(cbind, x)))) # Overwrite the old SortStructure with a vectorization from lst
    Indices <- matrix(SortStructure[, 3], ncol = NumberOfSurrogates) # Reshape the 3rd column of SortStructure into a matrix, using these inverse-rank numbers to overwrite the values in Indices
    Surrogates <- matrix(SortedValues[Indices], ncol = NumberOfSurrogates) # Overwrites Surrogates with a reshaping of SortedValues using the matrix of inverse ranks in Indices as indices
    nn <- nn + 1
    fftSurrogates <- mvfft(Surrogates) # Perform fast Fourier transformation concurrently on each column of Surrogates
    fftSurrogates.re <- Re(fftSurrogates) # Save the real component of fftSurrogates
    fftSurrogates.im <- Im(fftSurrogates) # Save the imaginary component of fftSurrogates
    fftRatio <- fftSurrogates.im / fftSurrogates.re # Compute the ratio of the imaginary and real Fourier components
    PhaseSpectra <- atan(fftRatio) # Computes phase spectrum as the arctangent of fftRatio
  }
  return(Surrogates)
}
