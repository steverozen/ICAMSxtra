#' Transcription bias of indels classified into 115 categories (pyrimidine)
#'
#' This data is designed to be used as an example in function \cr
#' \code{\link{PlotTransBiasID115}} and \code{\link{PlotTransBiasID115ToPdf}}.
#'
#' @format A \code{vector} which contains the 115 categories of indel events, standardised 
#' to pyrimidine format.
#'   
#' @name Target
NULL

#' @rdname Target
"target"

#' Transcription bias of indels classified into 115 categories (purine)
#'
#' This data is designed to be used as an example in function \cr
#' \code{\link{PlotTransBiasID115}} and \code{\link{PlotTransBiasID115ToPdf}}.
#'
#' @format A \code{vector} which contains the 115 categories of indel events, but in purine
#' format
#'   
#' @name Reverse
NULL

#' @rdname Reverse
"reverse"

#' Transcription bias of indels classified into 13 categories (pyrimidine)
#'
#' This data is designed to be used as an example in function \cr
#' \code{\link{PlotTransBiasID115}} and 
#' \code{\link{PlotTransBiasID115}} 
#'
#' @format A \code{vector} which contains the 13 categories of indel events, standardised 
#' to pyrimidine format.
#'   
#' @name Target_pooled
NULL

#' @rdname Target_pooled
"target_pooled"

#' Transcription bias of indels classified into 13 categories (purine)
#'
#' This data is designed to be used as an example in function \cr
#' \code{\link{PlotTransBiasID115}} and 
#' \code{\link{PlotTransBiasID115}} when \code{pool = TRUE}.
#'
#' @format A \code{vector} which contains the 13 categories of indel events, standardised 
#' to purine format.
#'   
#' @name Reverse_pooled
NULL

#' @rdname Reverse_pooled
"reverse_pooled"

#' Genic/intergenic size and proportions for H.sapiens BSgenome GRCh37
#' 
#' This data is designed to be used in function
#' \code{\link{PlotTranscriptionAssociatedDamageToPdf}} 
#' 
#' @format A list of 5 numbers with the names: 
#' 1. total.bp 2. coding.bp 3. noncoding.bp 4. prop.coding and 5. prop.noncoding
"GRCh37.proportions"

#' Genic/intergenic size and proportions for H.sapiens BSgenome GRCh38
#' 
#' This data is designed to be used in function
#' \code{\link{PlotTranscriptionAssociatedDamageToPdf}} 
#' 
#' @format A list of 5 numbers with the names: 
#' 1. total.bp 2. coding.bp 3. noncoding.bp 4. prop.coding and 5. prop.noncoding
"GRCh38.proportions"

#' Standard order of row names in a catalog
#'
#' This data is designed for those
#' who need to create their own catalogs from formats not
#' supported by this package. The rownames denote the mutation
#' types.  For example, for SBS96 catalogs, the rowname
#'  AGAT represents a mutation from AGA > ATA.
#'
#' @format A list of character vectors indicating the standard
#'   orders of row names in catalogs.
#'
#' @note In ID (small insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation deletion
#'   repeat sizes range from 1 to 6+. In ID83 catalogs, deletion repeat sizes
#'   range from 0 to 5.
#'
#' @name CatalogRowOrder
#' 
#' @examples 
#' catalog.row.order$ID115
#' # There are altogether 115 row names to denote the mutation types
#' # in ID115 catalog.
#' # The difference from the .$ID class in \code{\link{ICAMS::catalog.row.order}} is that 
#' # single base nonhomopolymer indels have trinucleotide context added to them in the format
#' # INS/DEL:C/T:1:0_PF where P is the preceding base and F is the following base.
"catalog.row.order"

# Quiets concerns of R CMD check about no visible binding for global variable
utils::globalVariables(c("N", "target_pooled", "target", "target_pooled",
                         "target", "trans.gene.symbol", ".", "REF", "ALT",
                         "trans.strand", "seq.context", "seq.context.width", 
                         "trans.start.pos", "trans.end.pos", "CHROM", "POS",
                         "reverse_pooled", "reverse", "mutation", 
                         "catalog.row.order", "mutation", "GRCh37.proportions",
                         "GRCh38.proportions", "CheckAndAssignAttributes", "cols"
))