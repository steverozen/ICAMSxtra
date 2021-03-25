# ICAMSxtra 0.0.4.90xx
## Added
* Added an extra argument `row.names` in function `WriteExposure`.

* Added extra argument `cex.yaxis` in functions `PlotExposure` and `PlotExposureToPdf`.

## Fixed
* Fixed a bug in function `PlotExposure` when the number of columns in exposure
matrix is less than `samples.per.line`.

* Fixed a bug in function `PlotExposureInternal` when user provides specified y axis label.

<br/>

# ICAMSxtra 0.0.3
## Added
* Extra argument `mean.weighted` for functions `MeanOfSpectraAsSig` and `PlotSpectraAsSigsWithUncertainty`.

* New dependency packages `boot` and `simpleboot`.

* New internal function `CalculateConfidenceInterval` to calculate bootstrap confidence interval.

* Extra argument `conf.int` in function `PlotSpectraAsSigsWithUncertainty` to
plot the error bars as confidence interval for the mean.

* Extra argument `num.of.bootstrap.replicates` in *exported* function `PlotSpectraAsSigsWithUncertainty`.

<br/>

# ICAMSxtra 0.0.2
## Added
* Functions for matching sets of signatures (see `match.sigs.R`).

<br/>

# ICAMSxtra 0.0.1
## Added
* Functions for reading, writing, sorting and plotting exposures were added: `ReadExposure`, `WriteExposure`, `SortExposure`, `PlotExposure`, `PlotExposureToPdf`.

* Various functions generating, analyzing and plotting ID115 catalogs. See the reference manual for more details.
