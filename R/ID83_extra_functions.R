###############################################################################
# Transcription-associated damage functions
###############################################################################

#' Plot indel counts on transcribed and nontranscribed strands to pdf
#'
#' @param list.of.vcfs List of in-memory ID VCFs. The list names will be
#' the sample ids in the output catalog.
#' 
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'   
#' @param names.of.vcfs list of names of vcfs
#' @param file The name of the PDF file to be produced.
#' 
#' @section Note: 
#' The strand bias statistics are Benjamini-Hochberg q-values based on two-sided
#' binomial tests of the mutation counts on the transcribed and untranscribed strands
#' relaitve to the actual abundances of C and T on the transcribed strand. On the
#' plot, asterisks indicate q-values as follows *, \eqn{Q<0.05}; **, \eqn{Q<0.01}; ***,
#'  \eqn{Q<0.001}.
#'
#' @return a list of tables of p-values for each vcf
#' @export
#'
#' @examples 
#' dirs <- c(system.file("extdata/Strelka-ID-vcf", "Strelka.ID.GRCh37.s1.vcf", package="ICAMSxtra"),
#'           system.file("extdata/Strelka-ID-vcf", "Strelka.ID.GRCh37.s2.vcf", package="ICAMSxtra"))
#' 
#' list.of.vcfs <- ICAMS::ReadStrelkaIDVCFs(dirs, names.of.VCFs = c("s1","s2"))
#' PlotTranscriptionAssociatedDamageToPdf(list.of.vcfs = list.of.vcfs, 
#'                                        ref.genome = "hg19", 
#'                                        names.of.vcfs = c("s1","s2"), 
#'                                        file = file.path(tempdir(), "test.pdf"))
PlotTranscriptionAssociatedDamageToPdf <- function(list.of.vcfs, ref.genome, names.of.vcfs, file){
  
  annotated.vcfs <- AnnotateIDVCFsWithTransRanges(list.of.vcfs, ref.genome, vcf.names = names.of.vcfs)
  
  mutation.annotated.vcfs <- lapply(FUN = function(dt){dt[, mutation := ICAMS:::CanonicalizeID(seq.context, 
                                                                                       REF, 
                                                                                       ALT, 
                                                                                       seq.context.width+1)]},
                                    annotated.vcfs)
  
  # tabulation of transcription-associated / not indels
  list.of.counts <- list()
  n <- length(mutation.annotated.vcfs)
  for (i in 1:n){
    counts <- matrix(data = 0L, nrow = 2, ncol = 83)
    rownames(counts) <- c("coding", "noncoding")
    colnames(counts) <- catalog.row.order$ID
    DT <- mutation.annotated.vcfs[[i]]
    m <- nrow(DT)
    for (j in 1:m){
      col <- DT[j, mutation]
      row <- ifelse(!is.na(DT[j, trans.gene.symbol]), "coding", "noncoding")
      counts[row, col] <- counts[row, col] + 1
    }
    list.of.counts[[names(mutation.annotated.vcfs)[i]]] <- counts
  }
  
  # define helper function
  t_col <- function(color, percent = 50, name = NULL) {
    #      color = color name
    #    percent = % transparency
    #       name = an optional name for the color
    
    ## Get RGB values for named color
    rgb.val <- col2rgb(color)
    
    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100 - percent) * 255 / 100,
                 names = name)
    
    ## Save the color
    invisible(t.col)
  }
  
  # calculate density based on counts
  list.of.densities <- list.of.counts
  for (i in 1:n){
    density <- list.of.densities[[i]]
    list.of.densities[[i]]["coding",] <- density["coding", ] * 1000000 / GRCh37.proportions$coding.bp;
    list.of.densities[[i]]["noncoding",] <- density["noncoding", ] * 1000000 / GRCh37.proportions$noncoding.bp;
  }
  
  # calculate p-values based on counts
  list.of.p.values <- list()
  for (i in 1:n){
    counts <- list.of.counts[[i]]
    p.values <- vector()
    for (category in colnames(counts)){
      htest <- binom.test(x = counts[, category], p = GRCh37.proportions$prop.coding, 
                          alternative = "two.sided")
      p.values[category] <- htest$p.value
    }
    list.of.p.values[[names(mutation.annotated.vcfs)[i]]] <- p.values
  }
  
  # plot to pdf
  grDevices::pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  par(mfrow = c(8, 1), mar = c(3, 4, 2, 2), oma = c(3, 2, 1, 1))
  
  class.col <- c("#fdbe6f",
                 "#ff8001",
                 "#b0dd8b",
                 "#36a12e",
                 "#fdcab5",
                 "#fc8a6a",
                 "#f14432",
                 "#bc141a",
                 "#d0e1f2",
                 "#94c4df",
                 "#4a98c9",
                 "#1764ab",
                 "#e2e2ef",
                 "#b6b6d8",
                 "#8683bd",
                 "#61409b")
  
  bg.col <- mapply(t_col, class.col, percent=90)
  
  strand.col <- c("#394398",
                  "#e83020")
  cols <- rep(strand.col, 83)
  cex <- par("cex")
  
  # indel specific parameters
  maj.class.names <- c("1bp deletion", "1bp insertion",
                       ">1bp deletions at repeats\n(Deletion length)",
                       ">1bp insertions at repeats\n(Insertion length)",
                       "Deletions with microhomology\n(Deletion length)")
  
  category.lab <- c(rep(c("C", "T"), 2), rep(c("2", "3", "4", "5+"), 3))
  category.col <- c(rep(c("black", "white"), 2),
                    rep(c("black", "black", "black", "white"), 3))
  
  for (i in 1:n){
    density <- list.of.densities[[i]]
    counts <- list.of.counts[[i]]
    ymax <- max(density) * 1.3
    
    bp <- barplot(unname(density), beside = TRUE, ylim = c(0, ymax),
                  axes = FALSE, lwd = 3, xaxs = "i",
                  border = NA, col = cols, xpd = NA, ylab = "mut/million",
                  cex.lab = cex * par("cex.lab") * 1.25, 
                  las = 2)
    
    # left and right x ticks
    x.left <- bp[2 * c(seq(0, 66, 6), 72, 73, 75, 78) + 1] - 0.5
    x.right <- bp[2 * c(seq(6, 72, 6), 73, 75, 78, 83)] + 0.5
    class.pos <- c((x.left[seq(1, 4, 2)] + x.right[seq(2, 5, 2)]) / 2,
                   (x.left[c(6, 10)] + x.right[c(8, 12)] - 12) / 2,
                   (x.left[13] + x.right[length(x.left)]) / 2)
    
    # Draw background color
    rect(xleft = x.left - 0.5, 0, xright = x.right + 0.5, ymax,
         col = bg.col, border = 'grey90', lwd = 1.5)
    
    # Plot again
    barplot(unname(density), beside = TRUE, ylim = c(0, ymax),
            axes = FALSE, ann = FALSE, lwd = 3,
            border = NA, col = cols, xpd = NA, add = TRUE)
    
    # Draw grid lines
    segments(x.left[1], seq(0, ymax, ymax / 4), x.right[length(x.right)],
             seq(0, ymax, ymax / 4), col = "grey60", lwd = 0.5, xpd = NA)
    
    # Draw y axis
    y.axis.values <- seq(0, ymax, length.out = 5)
    y.axis.labels <- round(y.axis.values, digits = 2)
    Axis(side = 2, at = y.axis.values, las = 1, labels = FALSE)
    text(-2, y.axis.values, labels = y.axis.labels,
         las = 1, adj = 1, xpd = NA, cex = cex)
    
    # Draw x axis labels
    mut.type <- c(rep(c("1", "2", "3", "4", "5", "6+"), 2),
                  rep(c("0", "1", "2", "3", "4", "5+"), 2),
                  rep(c("1", "2", "3", "4", "5", "6+"), 4),
                  rep(c("0", "1", "2", "3", "4", "5+"), 4),
                  "1", "1", "2", "1", "2", "3", "1", "2", "3", "4", "5+")
    bottom.pos <- c((x.left[1] + x.right[2]) / 2, (x.left[3] + x.right[4]) / 2,
                    class.pos[3 : length(class.pos)])
    bottom.lab <- c("Homopolymer length", "Homopolymer length",
                    "Number of repeat units", "Number of repeat units",
                    "Microhomology length")
    rect(xleft = x.left, -ymax * 0.09, xright = x.right, -ymax * 0.01,
         col = class.col, border = NA, xpd = NA)
    text((bp[1,]+bp[2,])/2, -ymax * 0.15, labels = mut.type, cex = 0.65, xpd = NA)
    text(bottom.pos, -ymax * 0.27, labels = bottom.lab, cex = 0.75, xpd = NA)
    
    # Draw lines above each class
    rect(xleft = x.left, ymax * 1.02, xright = x.right, ymax * 1.11,
         col = class.col, border = NA, xpd = NA)
    text((x.left + x.right) / 2, ymax * 1.06, labels = category.lab,
         cex = 0.65, col = category.col, xpd = NA)
    
    # Draw mutation class labels at the top of the figure
    text(class.pos, ymax * 1.27, labels = maj.class.names, cex = 0.75, xpd = NA)
    
    # Write the name of the sample
    text(3, ymax * 0.6, labels = names(list.of.densities)[[i]], adj = 0, cex = cex, font = 2)
    
    # Write mutation counts
    count.labels = mapply(function(cols){sum(counts[,cols])}, 
                          list(1:6, 7:12, 13:18, 19:24, 25:48, 49:72, 73:83))
    text(bp[2*c(6, 12, 18, 24, 48, 72, 83)], ymax * 0.92, labels = count.labels,
         adj = c(1, 1), xpd = NA, cex = cex)
    
    
    # Draw asterisks if significant
    p.values <- list.of.p.values[[i]]
    for (i in 1:83){
      if (p.values[[i]] < 0.05){
        x1 <- bp[1, i]
        x2 <- bp[2, i]
        y1 <- y2 <- max(density[, i]) + max(density) * 0.03
        segments(x1, y1, x2, y2)
        x3 <- mean(c(x1, x2))
        y3 <- y1 + max(density) * 0.03
        label <- ICAMS:::AssignNumberOfAsterisks(p.values[[i]])
        text(x3, y3, label)
      }
    }
  }
  grDevices::dev.off()
  retval <- list(plot.success=TRUE, list.of.p.values = list.of.p.values)
  return(retval)
}
