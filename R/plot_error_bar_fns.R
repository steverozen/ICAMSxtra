#' Return the mean of multiple spectra as a signature.
#' 
#' @param spectra An \code{\link[ICAMS]{ICAMS}} spectrum catalog.
#'   Convert each spectrum to a signature and then compute the
#'   mean.
#'   
#' @param title The name of the output signature.
#'   
#' @export
#' 
MeanOfSpectraAsSig <- function(spectra, title = "sig.from.spectra.mean") {
  
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
  
  mean.sig <- apply(X = sigs, MARGIN = 1, mean)
  
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
#' @export
#'
#' @return The mean of the spectra as a signature, the
#'   constituent spectra as signatures, and the y positions of the
#'   arrowheads.
PlotSpectraAsSigsWithUncertainty <- 
  function(spectra, title = "Mean.as.signature") {
    xx <- MeanOfSpectraAsSig(spectra = spectra, title = title)
    arrow.tops <- apply(xx$constituent.sigs, 1, max)
    arrow.bottoms <- apply(xx$constituent.sigs, 1, min)
    reval <- PlotCatalogWithArrows(xx$mean.sig[ , 1, drop = FALSE],
                                   arrow.tops, arrow.bottoms)
    inivisble(retval)
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
  xx$arrow.tops    <- arrow.tops
  xx$arrow.bottoms <- arrow.bottoms
  return(invisible(xx))
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
