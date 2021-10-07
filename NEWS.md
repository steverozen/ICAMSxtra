# ICAMSxtra 0.0.7
* Fixed a bug in plotting exposure legend angle and density.

# ICAMSxtra 0.0.6
* Fixed bugs related to usage of package data `catalog.row.order$ID115`.

# ICAMSxtra 0.0.5
* Updated code not to use deprecated functions in ICAMS.

<br/>

# ICAMSxtra 0.0.4
## Added
* Added an extra argument `row.names` in function `WriteExposure`.

* Added extra argument `cex.yaxis`, `cex.xaxis`, `plot.sample.names` and
`yaxis.labels` in functions `PlotExposureInternal`, `PlotExposure` and
`PlotExposureToPdf`.

* Added extra argument `width`, `height` in function `PlotExposureToPdf`.

## Fixed
* Fixed a bug in function `PlotExposure` when the number of columns in exposure
matrix is less than `samples.per.line`.

* Fixed a bug in function `PlotExposureInternal` when user provides specified y axis label.

* Fixed a bug in function `PlotExposureInternal` when user provides specified space (as a fraction of the average bar width) left before each bar.

* Fixed a bug in function `PlotExposureInternal` when user provides specified density of shading lines for the bars.

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
