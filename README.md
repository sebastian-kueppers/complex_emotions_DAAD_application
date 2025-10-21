# Complexity in Emotion ESM Data

In the R script at hand, I analyze two ESM-type datasets on emotions from [Bringmann et al. (2013)](https://doi.org/10.1371/journal.pone.0060188) and [Fried et al. (2022)](https://doi.org/10.1177/21677026211017839) for the prevalence of certain complexity markers. This analysis is inspired by [Olthof et al. (2020)]( https://doi.org/10.1186/s12916-020-01727-2) who run a similar investigation on the _n_ = 1 dataset from [Wichers & Groot (2016)](https://doi.org/10.1159/000441458).

### Complexity Markers
To each individual time series in the data, I apply the same markers [Olthof et al. (2020)](https://doi.org/10.1186/s12916-020-01727-2) employed: the Bartels rank test, the inspection of (partial) autocorrelation functions, and the application of a time-varying autoregressive model to test for the complex systems property _memory_; the KPSS test for level stationarity and a nonparametric change point analysis to test for _regime shifts_; and the forecast skill (Sugihara-May algorithm) to test for _sensitive dependence on initial conditions_. Additionally, I will include the KPSS test for trend stationarity as a test for regime shifts, as well as Keenan’s nonlinearity test, Tsay’s nonlinearity test, and the multifractality test ([Kelty-Stephen et al., 2022](https://doi.org/10.3758/s13428-022-01866-9)) to investigate the complex property _nonlinearity_. 

### Hypotheses
This research is largely exploratory, so no hypotheses are defined beforehand. However, I expect that most ESM time series will display complexity. Moreover, I expect longer time series to display more complexity, as more observations provide greater statistical power for the tests on complexity markers.

### Results
#### Bringmann et al. (2013)
43.5% of all significance tests (Bartels rank test, TV-AR model, KPPS tests, Tsay’s and Keenan’s nonlinearity tests, multifractality test) performed on participant x item time series detected a complexity marker (α = .05). On average, each time series displayed 0.63 change points, a forecast decay of −0.46, and 1.48 significant partial autocorrelations. On average, for a given time series, the significant partial correlation with the longest lag had a lag of 6.22 observations.
#### Fried et al. (2022)
27.7% of all significance tests found complexity markers (α = .05). On average, a particular time series exhibited zero change points, a forecast decay of -0.66, and 0.69 significant partial autocorrelations with the longest one having a lag of 5.13 observations.

### Conclusion
Complexity was prevalent in both datasets. As expected, the data from [Bringmann et al. (2013)](https://doi.org/10.1371/journal.pone.0060188) displayed more complexity than the data from [Fried et al. (2022)](https://doi.org/10.1177/21677026211017839) (apart from forecast decay; see chi squared and _t_ tests in R script).
