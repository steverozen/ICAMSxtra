#' Return the mean of multiple spectra as a signature
#' 
#' @param spectra An \code{\link[ICAMS]{ICAMS}} spectrum catalog.
#'   Convert each spectrum to a signature and then compute the
#'   mean.
#'
#' @param mean.weighted Logical. Whether to weigh the samples according to the
#'   number of mutations in them to calculate the weighted mean as the consensus
#'   signature. Default is TRUE. If FALSE, then arithmetic mean will be
#'   calculated.
#'   
#' @param title The name of the output signature.
#'   
#' @export
MeanOfSpectraAsSig <- 
  function(spectra, mean.weighted = TRUE, title = "sig.from.spectra.mean") {
  
  if (is.null(spectra)) {
    stop("Argument spectra is NULL")
  }
  
  ctype <- attr(spectra, "catalog.type", exact = TRUE)
  
  if (is.null(ctype)) {
    stop("Cannot call MeanOfSpectraAsSig when catalog.type is NULL")
  } else if (ctype == "counts") {
    tctype <- "counts.signature"
  } else if (ctype == "density") {
    tctype <- "density.signature"
  } else {
    stop("Cannot call MeanOfSpectraAsSig when catalog.type is ", ctype)
  }
  
  target.abundance <- attr(spectra, "abundance", exact = TRUE)
  target.region    <- attr(spectra, "region",    exact = TRUE)
  
  # stopifnot(!is.null(target.abundance))
  if (is.null(target.abundance)) {
    warning("Using NULL target.abundance in MeanOfSpectraAsSig")
  }
  
  sigs <- ICAMS::TransformCatalog(spectra, 
                                  target.catalog.type = tctype,
                                  target.abundance = target.abundance)
  
  if (mean.weighted == TRUE) {
    wts <- colSums(spectra)
    mean.sig <- apply(X = sigs, MARGIN = 1, stats::weighted.mean, w = wts)
  } else {
    mean.sig <- apply(X = sigs, MARGIN = 1, mean)
  }
  
  mean.sig <- matrix(mean.sig, ncol = 1)
  
  mean.sig <- 
    ICAMS::as.catalog(mean.sig, 
                      catalog.type   = tctype, 
                      region         = target.region, 
                      abundance      = target.abundance,
                      infer.rownames = TRUE)
  
  colnames(mean.sig) <- title
  
  return(list(mean.sig = mean.sig, constituent.sigs = sigs))
}

#' Convert spectra to signatures and then plot mean with "error" bars
#' 
#' @inheritParams  MeanOfSpectraAsSig
#' 
#' @param conf.int A number specifying the required confidence intervals. The
#'   error bars will be plotted as confidence intervals for the mean. If NULL,
#'   then use the maximum and minimum value to plot error bars.
#' 
#' @export
#'
#' @return The mean of the spectra as a signature, the
#'   constituent spectra as signatures, and the y positions of the
#'   arrowheads.
PlotSpectraAsSigsWithUncertainty <- 
  function(spectra, mean.weighted = TRUE, conf.int = 0.95, 
           title = "Mean.as.signature") {
    xx <- MeanOfSpectraAsSig(spectra = spectra, mean.weighted = mean.weighted,
                             title = title)
    if (!is.null(conf.int)) {
      retval <- CalculateConfidenceInterval(xx$constituent.sigs, 
                                            conf.int = conf.int)
      arrow.tops <- sapply(retval, FUN = "[", 2)
      arrow.bottoms <- sapply(retval, FUN = "[", 1)
    } else {
      arrow.tops <- apply(xx$constituent.sigs, 1, max)
      arrow.bottoms <- apply(xx$constituent.sigs, 1, min)
    }
    
    PlotCatalogWithArrows(xx$mean.sig[ , 1, drop = FALSE],
                          arrow.tops, arrow.bottoms)
    xx$arrow.tops    <- arrow.tops
    xx$arrow.bottoms <- arrow.bottoms
    return(invisible(xx))
  }

#' @importFrom boot boot.ci
#' @importFrom  simpleboot one.boot
#' @keywords internal
CalculateConfidenceInterval <- 
  function(constituent.sigs, num.of.bootstrap.replicates = 10^4, 
           conf.int = 0.95) {
    retval <- apply(X = constituent.sigs, MARGIN = 1, FUN = function(x) {
      simpleboot::one.boot(data = x, FUN = mean, 
                           R = num.of.bootstrap.replicates)
    })
    
    retval2 <- lapply(retval, FUN = function(x){
      # Get the bootstrap confidence interval
      boot::boot.ci(boot.out = x, conf = conf.int, type = "bca")$bca[4:5]
    })
    
    return(retval2)
  }

#' @keywords internal
PlotCatalogWithArrows <- function(catalog, arrow.tops, arrow.bottoms) {
  oldopt <- getOption("digits")
  on.exit(options(digits = oldopt))
  options(digits = 2)
  bp <- ICAMS::PlotCatalog(
    catalog = catalog,
    upper   = TRUE,
    xlabels = TRUE,
    cex     = 0.8,
    ylim    = c(min(arrow.bottoms), max(arrow.tops) + 0.005)
  )
  AddArrows(bp$plot.object, arrow.tops, arrow.bottoms)
}

#' @keywords internal
AddArrows <- function(bp, tops, bottoms) {
  oldopt <- getOption("warn")
  on.exit(options(warn = oldopt))
  options(warn = -1) # Does not seem to turn off warnings
  which0 <- which((tops - bottoms) == 0)
  tops[which0] <- tops[which0] + 1e-5 
  suppressWarnings(
    # Necessary because generates warnings for 0-length arrows
    arrows(
      x0     = bp,
      y0     = tops,    # location of up arrow tips
      y1     = bottoms, # location of down arrow tips
      angle  = 90,      # use "T" arrow heads
      code   = 3,       # use arrow heads at both ends
      length = 0.025    # width of arrow heads
    ))
  # Try this to get rid of the warnings:
  assign("last.warning", NULL, envir = baseenv())
  
}
