###############################################################################
# Utility functions
###############################################################################

#' Given a deletion and its sequence context, categorize it
#' 
#' This function is primarily for internal use, but we export it
#' to document the underlying logic.
#' 
#' See 
#' \url{https://github.com/steverozen/ICAMS/raw/master/data-raw/PCAWG7_indel_classification_2017_12_08.xlsx}
#' for additional information on deletion 
#' mutation classification.
#' 
#' This function first handles deletions in homopolymers, then
#' handles deletions in simple repeats with
#' longer repeat units, (e.g. \code{CACACACA}, see
#' \code{\link{FindMaxRepeatDel}}),
#' and if the deletion is not in a simple repeat,
#' looks for microhomology (see \code{\link{FindDelMH}}).
#' 
#' See the code for unexported function \code{\link{CanonicalizeID}}
#' and the functions it calls for handling of insertions.
#'
#' @param context The deleted sequence plus ample surrounding
#'   sequence on each side (at least as long as \code{del.seq}).
#'
#' @param del.seq The deleted sequence in \code{context}.
#'
#' @param pos The position of \code{del.sequence} in \code{context}.
#'
#' @param trace If > 0, then generate messages tracing
#' how the computation is carried out.
#' 
#' @param strand NULL by default. But when called by PlotTransBiasInternal,
#' \code{strand} is either \code{+} or \code{-}. The return value will 
#' include :trans or :nontrans indicating whether the deletion occured on
#' the transcribed or non-transcribed strand.
#' 
#' @return A string that is the canonical representation
#'  of the given deletion type. Return \code{NA} 
#'  and raise a warning if
#'  there is an un-normalized representation of
#'  the deletion of a repeat unit.
#'  See \code{FindDelMH} for details.
#'  (This seems to be very rare.)
#'  
#'  @keywords internal
Canonicalize1Del115 <- function(context, del.seq, pos, trace = 0, strand) {
  # Is the deletion involved in a repeat?
  rep.count <- FindMaxRepeatDel(context, del.seq, pos)
  
  rep.count.string <- ifelse(rep.count >= 5, "5+", as.character(rep.count))
  deletion.size <- nchar(del.seq)
  deletion.size.string <- 
    ifelse(deletion.size >= 5, "5+", as.character(deletion.size))
  
  # Category is "1bp deletion"
  if (deletion.size == 1) {
    if (rep.count == 0){
      if (del.seq %in% c("A","G")){
        context <- revc(context)
        pos <- nchar(context) - pos + 1
        del.seq <- substr(context,pos,pos)
        preceding.base <- substr(context,pos+1,pos+1)
        proceeding.base <- substr(context,pos-1,pos-1)
        retval <- paste0("DEL:",del.seq,":1:0_",preceding.base,proceeding.base)
        if (!is.na(strand)){
          retval <- paste0(retval, ifelse(strand=="+", ":nontrans", ":trans"))
        }
      } else {
        preceding.base <- substr(context,pos-1,pos-1)
        proceeding.base <- substr(context,pos+1,pos+1)
        retval <- paste0("DEL:",del.seq,":1:0_",preceding.base,proceeding.base)
        if (!is.na(strand)){
          retval <- paste0(retval, ifelse(strand=="+", ":trans", ":nontrans"))
        }
      }
      return(retval)
      
    } else {
      if (del.seq == "A") del.seq <- "T"
      if (del.seq == "G") del.seq <- "C"
      return(paste0("DEL:", del.seq, ":1:", rep.count.string))
    }
  }
  
  # Category is ">2bp deletion"
  if (rep.count > 0) {
    return(
      paste0("DEL:repeats:", deletion.size.string, ":", rep.count.string))
  }
  
  # We have to look for microhomology
  microhomology.len <- FindDelMH(context, del.seq, pos, trace = trace)
  if (microhomology.len == -1) {
    warning("Non-normalized deleted repeat ignored:",
            "\ncontext: ", context,
            "\ndeleted sequence: ", del.seq,
            "\nposition of deleted sequence: ", pos)
    return(NA)
  }
  
  if (microhomology.len == 0) {
    stopifnot(rep.count.string == 0)
    # Categorize and return non-repeat, non-microhomology deletion
    return(paste0("DEL:repeats:", deletion.size.string, ":0"))
  }
  
  microhomology.len.str <-
    ifelse(microhomology.len >= 5, "5+", as.character(microhomology.len))
  
  return(paste0(
    "DEL:MH:", deletion.size.string, ":", microhomology.len.str))
}


#' Given an insertion and its sequence context, categorize it.
#' @param context 
#' @param ins.sequence 
#' @param pos The position of \code{ins.sequence} in \code{context}.
#' @param trace 
#' @param strand
#'
#' @return A string that is the canonical representation of 
#' the given insertion type.
#' 
#' @details same as Canonicalize1INS except when 
#' \code{insertion.size == 1} & \code{rep.count == 0}
#' in which the format will be \code{INS:Y:1:0_PF}
#' where Y is C or T, P is the preceding nucleotide, and
#' F is the following nucleotide
#' 
#' @keywords internal
Canonicalize1INS115 <- function(context, ins.sequence, pos, trace = 0, strand) {
  if (trace > 0) {
    message("Canonicalize1ID(", context, ",", ins.sequence, ",", pos, "\n")
  }
  rep.count <- FindMaxRepeatIns(context, ins.sequence, pos)
  rep.count.string <- ifelse(rep.count >= 5, "5+", as.character(rep.count))
  insertion.size <- nchar(ins.sequence)
  insertion.size.string <-
    ifelse(insertion.size >= 5, "5+", as.character(insertion.size))
  
  if (insertion.size == 1) {
    if (rep.count == 0){
      if (ins.sequence %in% c("A","G")){
        context <-  revc(context)
        ins.sequence <-  revc(ins.sequence)
        pos <- nchar(context) - pos + 1
        preceding.base <- substr(context, pos, pos)
        proceeding.base <- substr(context, pos-1, pos-1)
        retval <- paste0("INS:",ins.sequence,":1:0_",preceding.base,proceeding.base)
        if (!is.na(strand)){
          retval <- paste0(retval, ifelse(strand=="+",":nontrans", ":trans"))
        }
      } else {
        preceding.base <- substr(context, pos, pos)
        proceeding.base <- substr(context, pos+1, pos+1)
        retval <- paste0("INS:",ins.sequence,":1:0_",preceding.base,proceeding.base)
        if (!is.na(strand)){
          retval <- paste0(retval, ifelse(strand=="+", ":trans", ":nontrans"))
        }
      }
    } else {
      if (ins.sequence == "A") ins.sequence <- "T"
      if (ins.sequence == "G") ins.sequence <- "C"
      retval <-
        (paste0("INS:", ins.sequence, ":1:", rep.count.string))
    }
  } else { 
    retval <-
      paste0("INS:repeats:", insertion.size.string, ":", rep.count.string)
  }
  if (trace > 0) message(retval)
  return(retval)
}

#' @title Given a single insertion or deletion in context categorize it.
#'
#' @param context Ample surrounding
#'   sequence on each side of the insertion or deletion.
#'
#' @param ref The reference allele (vector of length 1)
#'
#' @param alt The alternative allele (vector of length 1)
#'
#' @param pos The position of \code{ins.or.del.seq} in \code{context}.
#'
#' @param trace If > 0, then generate messages tracing
#' how the computation is carried out.
#' 
#' @param strand NULL by default. But when called by \code{ID115_StrandBiasGetPlottables},
#'  \code{strand} is either \code{+} or \code{-}. The return value will 
#'  include :trans or :nontrans indicating whether the deletion occured on
#'  the transcribed or non-transcribed strand.
#
#' @return A string that is the canonical representation
#'  of the type of the given
#'  insertion or deletion.
#'  Return \code{NA} 
#'  and raise a warning if
#'  there is an un-normalized representation of
#'  the deletion of a repeat unit.
#'  See \code{FindDelMH} for details.
#'  (This seems to be very rare.)
#'
#' @keywords internal
Canonicalize1ID115 <- function(context, ref, alt, pos, trace = 0, strand) {
  if (trace > 0) {
    message("Canonicalize1ID115(", context, ",", ref, ",", alt, ",", pos, "\n")
  }
  if (nchar(alt) < nchar(ref)) {
    # A deletion
    return(Canonicalize1Del115(context, ref, pos + 1, trace, strand))
  } else if (nchar(alt) > nchar(ref)) {
    # An insertion
    return(Canonicalize1INS115(context, alt, pos, trace, strand))
  } else {
    stop("Non-insertion / non-deletion found: ", ref, " ", alt, " ", context)
  }
}

#' @title Determine the mutation types of insertions and deletions.
#'
#' @param context A vector of ample surrounding
#'   sequence on each side the variants
#'
#' @param ref Vector of reference alleles
#'
#' @param alt Vector of alternative alleles
#'
#' @param pos Vector of the positions of the insertions and deletions in
#'  \code{context}.
#'  
#' @param strand NULL by default. But when called by \code{ID115_StrandBiasGetPlottables},
#'  \code{strand} is either \code{+} or \code{-}. The return value will 
#'  include :trans or :nontrans indicating whether the deletion occured on
#'  the transcribed or non-transcribed strand.
#' 
#' @return A vector of strings that are the canonical representations
#'  of the given insertions and deletions.
#'
#' @importFrom utils head
#'
#' @keywords internal
CanonicalizeID115 <- function(context, ref, alt, pos, strand) {
  
  if (all(substr(ref, 1, 1) == substr(alt, 1, 1))) {
    ref <- substr(ref, 2, nchar(ref))
    alt <- substr(alt, 2, nchar(alt))
  } else {
    stopifnot(ref != "" | alt != "")
  }
  ret <- mapply(Canonicalize1ID115, context, ref, alt, pos, 0, strand)  
  return(ret)
}

#' @title Create one column of the matrix for an indel catalog from *one* in-memory VCF.
#'
#' @param ID.vcf An in-memory VCF as a data.frame annotated by the
#'   \code{\link{AnnotateIDVCF}} function. It must only
#'   contain indels and must \strong{not} contain SBSs
#'   (single base substitutions), DBSs, or triplet
#'   base substitutions, etc.
#'
#'   One design decision for variant callers is the representation of "complex
#'   indels", e.g. mutations e.g. CAT > GC. Some callers represent this as C>G,
#'   A>C, and T>_. Others might represent it as CAT > CG. Multiple issues can
#'   arise. In PCAWG, overlapping indel/SBS calls from different callers were
#'   included in the indel VCFs.
#'
#' @param SBS.vcf This argument defaults to \code{NULL} and
#'   is not used. Ideally this should be an in-memory SBS VCF 
#'   as a data frame. The rational is that for some data,
#'   complex indels might be represented as an indel with adjoining
#'   SBSs. 
#'
#' @return A list of a 1-column matrix containing the mutation catalog
#'   information and the annotated VCF with ID categories information added.
#'
#' @keywords internal
CreateOneColID115Matrix <- function(ID.vcf, SBS.vcf = NULL) {
  if (nrow(ID.vcf) == 0) {
    # Create 1-column matrix with all values being 0 and the correct row labels.
    catID <- matrix(0, nrow = length(ICAMS::catalog.row.order$ID115), ncol = 1)
    rownames(catID) <- ICAMS::catalog.row.order$ID115
    return(catID)
  }
  
  if (!is.null(SBS.vcf)) 
    warning("Argument SBS.vcf in CreateOneColIDMatrix is always ignored")
  
  canon.ID <- CanonicalizeID115(ID.vcf$seq.context,
                                ID.vcf$REF,
                                ID.vcf$ALT,
                                ID.vcf$seq.context.width + 1,
                                strand = NA)
  
  if (any(is.na(canon.ID))) warning("NA ID categories ignored")
  idx <- which(is.na(canon.ID))
  canon.ID <- canon.ID[!is.na(canon.ID)]
  if (length(idx) > 0) {
    ID.vcf <- ID.vcf[-idx, ]
  }
  
  out.ID.vcf <- cbind(ID.vcf, ID.class = canon.ID)
  
  # Create the ID catalog matrix
  tab.ID <- table(canon.ID)
  
  row.order <- data.table(rn = ICAMS::catalog.row.order$ID115)
  
  ID.dt <- as.data.table(tab.ID)
  # ID.dt has two columns, names cannon.dt (from the table() function
  # and N (the count)
  
  ID.dt2 <-
    merge(row.order, ID.dt, by.x = "rn", by.y = "canon.ID", all = TRUE)
  ID.dt2[ is.na(N) , N := 0]
  stopifnot(setequal(unlist(ID.dt2$rn), ICAMS::catalog.row.order$ID115))
  
  ID.mat <- as.matrix(ID.dt2[ , 2])
  rownames(ID.mat) <- ID.dt2$rn
  return(list(catalog = ID.mat[ICAMS::catalog.row.order$ID115, , drop = FALSE],
              annotated.VCF = out.ID.vcf))
}


#' Create a catalog from a \code{matrix}, \code{data.frame}, or \code{vector}
#'
#' @param object A numeric \code{matrix}, numeric \code{data.frame},
#' or \code{vector}.
#' If a \code{vector}, converted to a 1-column \code{matrix}
#' with rownames taken from the element names of the \code{vector}
#' and with column name \code{"Unknown"}.
#' If argument \code{infer.rownames}
#'  is \code{FALSE} than this argument must have
#'   rownames to denote the mutation types. See \code{\link{CatalogRowOrder}}
#'   for more details.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param region A character string designating a region, one of
#' \code{genome}, \code{transcript}, \code{exome}, \code{unknown};
#' see \code{\link{ICAMS}}.
#'
#' @param catalog.type One of "counts", "density", "counts.signature",
#'   "density.signature".
#'
#' @param abundance If \code{NULL}, then
#'  inferred if \code{ref.genome}
#'  is one of
#'  the reference genomes known to ICAMS and \code{region}
#'  is not \code{unknown}. See \code{\link{ICAMS}}.
#'  The argument \code{abundance} should
#'  contain the counts of different source sequences for mutations
#'  in the same format as the numeric vectors in \code{\link{all.abundance}}.
#'  
#' @param infer.rownames If \code{TRUE}, and \code{object} has no
#' rownames, then assume the rows of \code{object} are in the
#' correct order and add the rownames implied by the number of rows
#' in \code{object} (e.g. rownames for SBS 192 if there are 192 rows).
#' If \code{TRUE}, \strong{be sure the order of rows is correct.}
#'
#' @return A catalog as described in \code{\link{ICAMS}}.
#'
#' @export
as.catalog.for.ID115 <- function(object, 
                                 ref.genome = NULL, 
                                 region = "unknown", 
                                 catalog.type = "counts", 
                                 abundance = NULL,
                                 infer.rownames = FALSE) {
  if (!is.matrix(object)) {
    if (is.data.frame(object)) {
      object <- as.matrix(object)
    } else if (is.vector(object)) {
      obj2 <- matrix(object, ncol = 1)
      rownames(obj2) <- names(object)
      colnames(obj2) <- "Unknown"
      object <- obj2
    } else {
      stop("object must be numeric matrix, vector, or data frame")
    }
  }
  stopifnot(mode(object) == "numeric")
  
  if (is.null(rownames(object))) {
    if (!infer.rownames) {
      stop("Require correct rownames on object unless infer.rownames == TRUE")
    }
    rownames(object) <- InferRownames(object)
  } else {
    #removed call to checkandreorder rownames
    object <- object[ICAMS::catalog.row.order$ID115, ,drop=FALSE] 
  }
  
  # StopIfRegionIllegal(region)
  
  class.string  <- "ID115Catalog"
  
  StopIfCatalogTypeIllegal(catalog.type)
  
  if (!is.null(ref.genome)) {
    ref.genome <- NormalizeGenomeArg(ref.genome)
  }
  
  attr(object, "ref.genome") <- ref.genome
  
  attr(object, "catalog.type") <- catalog.type
  
  attr(object, "abundance") <- NULL
  
  class(object) <- append(class(object), class.string, after = 0)
  class(object) <- unique(attributes(object)$class)
  
  attr(object, "region") <- region
  
  return(object)
}


###############################################################################
# Catalog functions 
###############################################################################

#' Create ID (small insertion and deletion) catalog from ID VCFs
#'
#' @param list.of.vcfs List of in-memory ID VCFs. The list names will be
#' the sample ids in the output catalog.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param region A character string acting as a region identifier, one of
#' "genome", "exome".
#' 
#' @param flag.mismatches Optional. If > 0, then if there are mismatches to
#'   references in the ID (insertion/deletion) VCF, generate messages showing
#'   the mismatched rows and continue. Otherwise \code{stop} if there are
#'   mismatched rows. See \code{\link{AnnotateIDVCF}} for more details.
#'
#' @return A list of two elements. 1st element \code{catalog} is the ID (small
#'   insertion and deletion) catalog with attributes added. See
#'   \code{\link{as.catalog}} for more details. 2nd element
#'   \code{annotated.vcfs} is a list of data frames which contain the original
#'   VCF with three additional columns \code{seq.context.width},
#'   \code{seq.context} and \code{ID.class} added. The category assignment of
#'   each ID mutation in VCF can be obtained from \code{ID.class} column.
#'   
#' @note In ID (small insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation
#'   deletion repeat sizes range from 1 to 6+.
#'   
#' @export
VCFsToID115Catalogs<-function(list.of.vcfs, ref.genome, region = "unknown",
                              flag.mismatches = 0){
  ncol <- length(list.of.vcfs)
  
  # Create a 0-column matrix with the correct row labels.
  catID <- matrix(0, nrow = length(ICAMS::catalog.row.order$ID115), ncol = 0)
  rownames(catID) <- ICAMS::catalog.row.order$ID115
  out.list.of.vcfs <- list()
  
  vcf.names <- names(list.of.vcfs)
  for (i in 1:ncol) {
    ID <- list.of.vcfs[[i]]
    list <- AnnotateIDVCF(ID, ref.genome = ref.genome,
                          flag.mismatches = flag.mismatches,
                          name.of.VCF = vcf.names[i])
    # Unlike the case for SBS and DBS, we do not
    # add transcript information.
    tmp <- CreateOneColID115Matrix(list$annotated.vcf)
    one.ID.column <- tmp[[1]]
    out.list.of.vcfs <- c(out.list.of.vcfs, list(tmp[[2]]))
    rm(ID)
    catID <- cbind(catID, one.ID.column)
  }
  
  colnames(catID) <- names(list.of.vcfs)
  names(out.list.of.vcfs) <- names(list.of.vcfs)
  
  return(list(catalog = 
                as.catalog.for.ID115(catID, ref.genome = ref.genome,
                                     region = region, catalog.type = "counts"),
              annotated.vcfs = out.list.of.vcfs))
}

#' "Collapse" a catalog
#' 
#' @description 
#' 
#' \code{Collapse115CatalogTo83} Collapse a ID 115 catalog
#' to a ID 83 catalog.
#'
#' @param catalog A catalog as defined in \code{\link{ICAMS}}.
#'
#' @return A catalog as defined in \code{\link{ICAMS}}.
#' 
#' @export
Collapse115CatalogTo83 <- function(catalog) {
  dt115 <- data.table(catalog)
  list83 <- list()
  list83[1] <- sum(dt115[1:9])
  list83[7] <- sum(dt115[15:23])
  list83[13] <- sum(dt115[29:37])
  list83[19] <- sum(dt115[43:51])
  for (i in 0:4){
    list83[i+2] <- dt115[i+10]
    list83[i+8] <- dt115[i+24]
    list83[i+14] <- dt115[i+38]
  }
  for (i in 0:63){
    list83[i+20] <- dt115[i+52] 
  }
  mat83 <- matrix(unlist(list83),ncol=1,byrow=TRUE)
  rownames(mat83) <- ICAMS::catalog.row.order$ID
  mat83 <- mat83[ICAMS::catalog.row.order$ID, , drop = FALSE]
  
  cat83 <-
    as.catalog(
      object = mat83,
      ref.genome = attributes(catalog)$ref.genome,
      region = attributes(catalog)$region,
      catalog.type = attributes(catalog)$catalog.type,
      abundance = NULL)
  
  #hdp0 etc... otherwise plot will not have names
  attributes(cat83)$dimnames[[2]] <- attributes(catalog)$dimnames[[2]]
  
  return(cat83)
}

#' "Collapse" a matrix of ID 115 catalogs to ID 83 catalog
#'
#' @param catalogs A catalog as defined in \code{\link{ICAMS}}.
#'
#' @return A catalog as defined in \code{\link{ICAMS}}.
#' @export
CollapseID115CatalogsToID83s <- function(catalogs){
  n <- ncol(catalogs)
  cat115s <-catalogs[, 1, drop=FALSE]
  cat83s <- Collapse115CatalogTo83(cat115s)
  for (i in 2:n){
    cat115_i <- catalogs[, i, drop=FALSE]
    cat83_i <-Collapse115CatalogTo83(cat115_i)
    cat83s <- cbind(cat83s,cat83_i)
  }
  return (cat83s)
}

###############################################################################
# Read and write functions
###############################################################################

#' @title Write a catalog to a file.
#' 
#' @param catalog A catalog as defined in \code{\link{ICAMS}} with attributes added.
#' See \code{\link{as.catalog}} for more details.
#'
#' @param file The path of the file to be written.
#' 
#' @param strict If TRUE, then stop if additional checks on the input fail.
#'
#' @export
WriteID115Catalog <- function(catalog, file, strict = TRUE) {
  WriteCat(catalog, file, 115, ICAMS::catalog.row.order$ID115,
           catalog.row.headers.ID115, strict)
}

#' Read catalog
#'
#' Read a catalog in standardized format from path.
#'
#' See also \code{\link{WriteCatalog}}
#'
#' @param file Path to a catalog on disk in the standardized format.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param region region A character string designating a genomic region;
#'  see \code{\link{as.catalog}} and \code{\link{ICAMS}}.
#'
#' @param catalog.type One of "counts", "density", "counts.signature",
#'   "density.signature".
#'
#' @return A catalog as an S3 object; see \code{\link{as.catalog}}.
#'
#' @note In ID (small insertion and deletion) catalogs, deletion repeat sizes
#'   range from 0 to 5+, but for plotting and end-user documentation
#'   deletion repeat sizes range from 1 to 6+.
#'
#' @section Comments:
#' To add or change attributes of the catalog, you can use function
#' \code{\link[base]{attr}}. \cr For example, \code{attr(catalog, "abundance")
#' <- custom.abundance}.
#' 
#' @export
ReadID115Catalog <- function(file, ref.genome = NULL, region = "unknown", 
                             catalog.type = "counts") {
  cos <- data.table::fread(file)
  stopifnot(nrow(cos) == 115)
  
  rn <- apply(cos[ , 1:4], MARGIN = 1, paste, collapse = ":")
  out <- as.matrix(cos[ , -(1:4), drop = FALSE])
  
  stopifnot(setdiff(rn, ICAMS::catalog.row.order$ID115) == c())
  stopifnot(setdiff(ICAMS::catalog.row.order$ID115, rn) == c())
  stopifnot(rn == ICAMS::catalog.row.order$ID115)
  rownames(out) <- rn
  #   if (ncol(out) == 1) colnames(out) <- colnames(cos)[3] 
  out <- out[ICAMS::catalog.row.order$ID115, , drop = FALSE]
  return(as.catalog.for.ID115(out, ref.genome, region, catalog.type))
}

###############################################################################
# Plotting functions for ID(insertion and deletion) catalog start here
###############################################################################

#' Plot \strong{one} spectrum or signature
#' 
#' Plot the spectrum of \strong{one} sample or plot \strong{one} signature. The
#' type of graph is based on one attribute("catalog.type") of the input catalog.
#' You can first use \code{\link{TransformCatalog}} to get different types of
#' catalog and then do the plotting.
#' 
#' @param catalog A catalog as defined in \code{\link{ICAMS}} with attributes added.
#' See \code{\link{as.catalog}} for more details.
#'
#' @param ylim Has the usual meaning. Only implemented for SBS96Catalog and
#'   IndelCatalog.
#'
#' @export
PlotID115Catalog <- function(catalog, ylim = NULL){
  stopifnot(dim(catalog) == c(115, 1))
  
  indel.class.col <- c("#fdbe6f",
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
  
  num.classes <- length(catalog)
  cols <- rep(indel.class.col,
              c(14, 14, 14, 14,
                6, 6, 6, 6,
                6, 6, 6, 6,
                1, 2, 3, 5))
  
  to.plot <- catalog[, 1]
  
  catalog.type <- attributes(catalog)$catalog.type
  if (catalog.type == "counts") {
    # Set a minimum value for ymax to make the plot more informative
    ymax <- 4 * ceiling(max(max(to.plot) * 1.3, 10) / 4)
  } else if (catalog.type == "counts.signature") {
    ymax <- ifelse(max(to.plot) * 1.3 > 1, 1, max(to.plot) * 1.3)
  } else {
    stop('\nCan only plot IndelCatalog with "counts" or "counts.signature" catalog.type.')
  }
  if (is.null(ylim)) {
    ylim <- c(0, ymax)
  } else {
    ymax <- ylim[2]
  }
  
  # Barplot
  bp <- barplot(catalog[, 1], ylim = c(0, ymax), axes = FALSE, xaxt = "n",
                lwd = 3, space = 1.35, border = NA, col = cols, xpd = NA,
                xaxs = "i", yaxt = "n")
  
  if (catalog.type == "counts") {
    # Calculate and draw the total counts for each major type
    counts <- integer(16)
    idx <- c(14, 28, 42, 56, seq(62, 104, 6), 105, 107, 110, 115)
    # midpoints of major types
    # c(7, 21, 35, 49, 59.5, 65.5, 71.5, 77.5, 83.5, 89.5, 95.5, 101.5, 105, 106.5, 109, 112.5)
    midx <- c(seq(7, 49, 14), seq(59.5, 101.5, 6), 105, 106.5, 109, 112.5)
    # scale midpoints to horizontal width of plot
    idx2 <- midx * 270 / 115
    for (i in 1:16) {
      if (i == 1) {
        counts[i] <- sum(catalog[1:idx[1], 1])
      } else {
        counts[i] <- sum(catalog[(idx[i - 1] + 1):idx[i], 1])
      }
      text(idx2[i], ymax * 0.6, labels = counts[i],
           cex = 0.68, adj = 1, xpd = NA)
    }
  } 
  
  # Draw box and grid lines
  rect(xleft = bp[1] - 1.0, 0, xright = bp[num.classes] + 1.0, ymax, col = NA,
       border = "grey60", lwd = 0.5, xpd = NA)
  segments(bp[1] - 1.0, seq(0, ymax, ymax / 4), bp[num.classes] + 1,
           seq(0, ymax, ymax / 4), col = "grey60", lwd = 0.5, xpd = NA)
  
  # Draw mutation class labels and lines above each class
  maj.class.names <- c("1bp deletion", "1bp insertion",
                       ">1bp deletions at repeats\n(Deletion length)",
                       ">1bp insertions at repeats\n(Insertion length)",
                       "Deletions with microhomology\n(Deletion length)")
  x.left <- bp[c(seq(1, 57, 14), seq(63, 99, 6), 105, 106, 108, 111)] - 0.8
  x.right <- bp[c(seq(14, 56, 14), seq(62, 104, 6), 105, 107, 110, 115)] + 0.8
  class.pos <- c((x.left[seq(1, 4, 2)] + x.right[seq(2, 5, 2)]) / 2,
                 (x.left[c(6, 10)] + x.right[c(8, 12)] - 12) / 2,
                 (x.left[13] + x.right[length(x.left)]) / 2)
  category.lab <- c(rep(c("C", "T"), 2), rep(c("2", "3", "4", "5+"), 3))
  category.col <- c(rep(c("black", "white"), 2),
                    rep(c("black", "black", "black", "white"), 3))
  
  # Draw lines above each class
  # draw at ymax instead of ymax
  rect(xleft = x.left, ymax * 1.02, xright = x.right, ymax * 1.09,
       col = indel.class.col, border = NA, xpd = NA)
  text((x.left + x.right) / 2, ymax * 1.05, labels = category.lab,
       cex = 0.65, col = category.col, xpd = NA)
  
  # Draw mutation class labels at the top of the figure
  text(class.pos, ymax * 1.2, labels = maj.class.names, cex = 0.75, xpd = NA)
  
  # Draw the sample name information of the sample
  text(1.5, ymax * 7 / 8, labels = colnames(catalog),
       adj = 0, cex = 0.8, font = 2)
  
  # Draw y axis
  y.axis.values <- seq(0, ymax, ymax / 4)
  if (attributes(catalog)$catalog.type != "counts") {
    y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
    text(-9, ymax / 2, labels = "counts proportion",
         srt = 90, xpd = NA, cex = 0.8)
  } else {
    y.axis.labels <- y.axis.values 
    text(-8, ymax / 2, labels = "counts",
         srt = 90, xpd = NA, cex = 0.8)
  }
  text(0, y.axis.values, labels = y.axis.labels,
       las = 1, adj = 1, xpd = NA, cex = 0.75)
  
  # Draw x axis labels
  mut.type <- c("AA","AT","AG","TA","TT","TG","GA","GT","GG",
                "2", "3", "4", "5", "6+",
                "AA","AC","AG","CA","CC","CG","GA","GC","GG",
                "2", "3", "4", "5", "6+",
                "AA","AT","AG","TA","TT","TG","GA","GT","GG",
                "1", "2", "3", "4", "5+",
                "AA","AC","AG","CA","CC","CG","GA","GC","GG",
                "1", "2", "3", "4", "5+",
                rep(c("1", "2", "3", "4", "5", "6+"), 4),
                rep(c("0", "1", "2", "3", "4", "5+"), 4),
                "1", "1", "2", "1", "2", "3", "1", "2", "3", "4", "5+")
  bottom.pos <- c((x.left[1] + x.right[2]) / 2, (x.left[3] + x.right[4]) / 2,
                  class.pos[3 : length(class.pos)])
  bottom.lab <- c("Homopolymer length", "Homopolymer length",
                  "Number of repeat units", "Number of repeat units",
                  "Microhomology length")
  rect(xleft = x.left, -ymax * 0.09, xright = x.right, -ymax * 0.02,
       col = indel.class.col, border = NA, xpd = NA)
  text(bp, -ymax * 0.11, labels = mut.type, cex = 0.65, xpd = NA, srt = 270, adj = 0)
  text(bottom.pos, -ymax * 0.27, labels = bottom.lab, cex = 0.75, xpd = NA)
  
  return(list(plot.success = TRUE))
}

#' Plot catalog to a PDF file
#'
#' Plot catalog to a PDF file. The type of graph is based on one
#' attribute("catalog.type") of the input catalog. You can first use
#' \code{\link{TransformCatalog}} to get different types of catalog and then do
#' the plotting.
#' 
#' @param catalog A catalog as defined in \code{\link{ICAMS}} with attributes added.
#' See \code{\link{as.catalog}} for more details.
#'
#' @param file The name of the PDF file to be produced.
#' @param ylim Has the usual meaning. Only implemented for SBS96Catalog and IndelCatalog.
#'
#' @return A list whose first element is a logic value indicating whether the
#'   plot is successful. For \strong{SBS192Catalog} with "counts" catalog.type
#'   and non-null abundance and \code{plot.SBS12 = TRUE}, the list will have a
#'   second element which is a list containing the strand bias statistics.
#'   
#' @note The sizes of repeats involved in deletions range from 0 to 5+ in the
#'   mutational-spectra and signature catalog rownames, but for plotting and
#'   end-user documentation deletion repeat sizes range from 1 to 6+.
#' 
#' @export
PlotID115CatalogToPdf <-
  function(catalog, file, ylim = NULL) {
    # Setting the width and length for LANDSCAPE A4 size plotting
    grDevices::pdf(file, width = 11.6929, height = 8.2677, onefile = TRUE)
    
    n <- ncol(catalog)
    opar <- par(mfrow = c(4, 1), mar = c(3, 4, 2.5, 2), oma = c(3, 3, 2, 2))
    on.exit(par(opar))
    
    for (i in 1 : n) {
      cat <- catalog[, i, drop = FALSE]
      PlotID115Catalog(cat, ylim = ylim)
    }
    grDevices::dev.off()
    return(list(plot.success = TRUE))
  }


#' @title Plot an ID 115 signatures (default) or catalog as standard ID83 and save as pdf file
#'
#' @param file The name of the PDF file to be produced.
#' 
#' @param ylim Has the usual meaning. Only implemented for SBS96Catalog and
#'   IndelCatalog.
#'   
#' @param catalog A catalog as defined in \code{\link{ICAMS}}.
#' 
#' @return A list whose first element is a logic value indicating whether the
#'   plot is successful. For \strong{SBS192Catalog} with "counts" catalog.type
#'   and non-null abundance and \code{plot.SBS12 = TRUE}, the list will have a
#'   second element which is a list containing the strand bias statistics.
#' 
#' @export 
PlotID115AsID83ToPdf <-
  function(catalog, file, ylim = NULL) {
    # Setting the width and length for LANDSCAPE A4 size plotting
    grDevices::pdf(file, width = 11.6929, height = 8.2677, onefile = TRUE)
    
    n <- ncol(catalog)
    opar <- par(mfrow = c(4, 1), mar = c(3, 4, 2.5, 2), oma = c(3, 3, 2, 2))
    on.exit(par(opar))
    
    for (i in 1 : n) {
      cat115 <- catalog[, i, drop = FALSE]
      cat83 <- Collapse115CatalogTo83(cat115)
      PlotCatalog.IndelCatalog(cat83, ylim = ylim)
    }
    grDevices::dev.off()
    return(list(plot.success = TRUE))
  }

###############################################################################
# Combined functions
###############################################################################

#' @title Read a list of vcfs and plot ID115 catalogs as pdf
#' 
#' @param list.of.vcfs List of in-memory ID VCFs. The list names will be
#' the sample ids in the output catalog.
#'
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'   
#' @param region A character string acting as a region identifier, one of
#' "genome", "exome".
#' 
#' @param flag.mismatches Optional. If > 0, then if there are mismatches to
#'   references in the ID (insertion/deletion) VCF, generate messages showing
#'   the mismatched rows and continue. Otherwise \code{stop} if there are
#'   mismatched rows. See \code{\link{AnnotateIDVCF}} for more details.
#'   
#' @param file The name of the PDF file to be produced.
#' 
#' @param ylim Has the usual meaning. Only implemented for SBS96Catalog and
#'   IndelCatalog.
#'
#' @export
VCFsToID115CatalogsAndPlotToPdf <-
  function(list.of.vcfs, ref.genome, region = "unknown",
           flag.mismatches = 0, file, ylim = NULL){
    catalogs <- VCFsToID115Catalogs(list.of.vcfs, ref.genome, region)$catalog
    PlotID115CatalogToPdf(catalogs, file, ylim)
  }

###############################################################################
# Strand bias functions
###############################################################################

#' @keywords internal
revcID115 <- function(string) {
  if (string %in% c(target_pooled, reverse_pooled)){
    idx_py <- which(grepl(string, target_pooled, fixed=TRUE))
    idx_pu <- which(grepl(string, reverse_pooled, fixed=TRUE))
    if (length(idx_py)==1){
      return(reverse_pooled[idx_py])
    }
    if (length(idx_pu)==1){
      return(target_pooled[idx_pu])
    }
    stop(paste0("Tried to reverse ", string," but it was not recognized\n"))
  } else {
    idx_py <- which(grepl(string, target, fixed=TRUE))
    idx_pu <- which(grepl(string, reverse, fixed=TRUE))
    if (length(idx_py)==1){
      return(reverse[idx_py])
    }
    if (length(idx_pu)==1){
      return(target[idx_pu])
    }
    stop(paste0("Tried to reverse ", string," but it was not recognized\n"))
  }
}

#' @keywords internal
ID115_StrandBiasGetPlottables <- 
  function(annotated.ID.vcf, damaged.base, vcf.name, pool) {
    df <- annotated.ID.vcf
    df <- df[!is.na(trans.gene.symbol), ]
    
    df1 <- df[, .(REF = REF[1], ALT= ALT[1], trans.strand = trans.strand[1],
                  seq.context = seq.context[1], seq.context.width = seq.context.width[1],
                  trans.gene.symbol = trans.gene.symbol[1],
                  trans.start.pos = trans.start.pos[1], 
                  trans.end.pos = trans.end.pos[1]), 
              by = .(CHROM, POS)]
    
    if (pool) {
      target_i <- target_pooled
      reverse_i <- reverse_pooled
    } else {
      target_i <- target
      reverse_i <- reverse
    }
    
    # Construct a matrix that later can be used to plot transcriptional strand bias
    result <- matrix(data = 0, nrow = 1, ncol = 2*length(target_i))
    rownames(result) <- c(1)
    
    mutation.type <- c(target_i, reverse_i)
    colnames(result) <- mutation.type
    
    if (nrow(df1) > 0){
      df1$mutation <- CanonicalizeID115(df1$seq.context, df1$REF, df1$ALT, df1$seq.context.width + 1, strand = df1$trans.strand)
      
      # filter for non-homopolymer single base indels
      df1 <- df1[which(grepl(pattern = "trans", x = df1$mutation)),,]
      
      # if aren't any non-homopolymer single base indels
      # then return value is same as if original df1 == 0
      if (nrow(df1) > 0){
        if (pool){
          df1$mutation1 <- substr(df1$mutation, 1, 5)
          df1$mutation2 <- substr(df1$mutation, 14, nchar(df1$mutation))
          df1$mutation <- paste0(df1$mutation1, ":", df1$mutation2)
          df1 <- df1[, c(-12,-13)]
        }
        
        for (i in 1:length(target_i)) {
          type <- mutation.type[i]
          idx <- which(unname(mapply(grepl, type, df1$mutation)))
          df2 <- df1[idx,,drop=FALSE]
          if (nrow(df2) == 0) {
            result[, type] <- 0
            result[, revcID115(type)] <- 0
          } else {
            nrow.nontrans <- length(which(grepl("nontrans", df2$mutation,fixed=TRUE)))
            result[1, revcID115(type)] <- nrow.nontrans
            result[1, type] <- nrow(df2) - nrow.nontrans
          }
        }  
      }
    }
    # Carry out logistic regression and get the p-values
    p.values <- ID115_CalculatePValues(df1, pool = pool)
    
    return(list(plotmatrix = result, vcf.df = df1, 
                p.values = p.values, vcf.name = vcf.name))
  }

#' @importFrom stats binom.test p.adjust
#' @keywords internal
ID115_CalculatePValues <- function(dt, pool) {
  if (pool){
    target_i <- target_pooled
  } else {
    target_i <- target
  }
  p.values <- numeric(length(target_i))
  names(p.values) <- target_i
  
  #construct a matrix with 4 (or 36) columns and 2 rows
  mat <- matrix(0L, nrow = 2, ncol = length(target_i))
  colnames(mat) <- target_i
  rownames(mat) <- c("trans", "nontrans")
  if (nrow(dt)> 0){
    for (i in 1:nrow(dt)){
      mutation <- dt[i, mutation]
      if (pool){
        col.name <- substr(mutation, 1, 5)
        row.name <- substr(mutation, 7, nchar(mutation))
      } else {
        col.name <- substr(mutation, 1, 12)
        row.name <- substr(dt[i, mutation], 14, nchar(mutation))
      }
      mat[row.name, col.name] <- mat[row.name, col.name] + 1
    }
  }
  
  # proportion of indels on transcribed strand of genes
  prop.transcribed <- 0.5
  
  for (type in target_i){
    htest <- binom.test(x = mat[, type], p = prop.transcribed, 
                        alternative = "two.sided")
    p.values[type] <- htest$p.value
  }
  
  # Adjust p-values for multiple comparisons to Benjamini-Hochberg false discovery rate
  q.values <- p.adjust(p.values, method = "BH")
  return(q.values)
}

#' @keywords internal
ID115_PlotTransBiasInternal <- function(list, ymax = NULL) {
  result <- list$plotmatrix
  # determined whether to plot pooled version based on ncol of plotmatrix
  if (ncol(result) == 8){
    pool = TRUE
  } else {
    if (ncol(result) == 72){
      pool = FALSE
    }  else {
      stop("Error in `list` argument passed to ID115_PlotTransBiasInternal. ncol of plotmatrix is neither 8 nor 72")
    }
  }
  if (is.null(ymax)) {
    ymax <- max(result)
  }
  if (pool){
    tmp <- matrix(0, nrow = 2, ncol = 4)
    tmp[1,] <- result[1:4]
    tmp[2,] <- result[5:8]
    bpcol = rep(c('#394398', '#e83020'), 4)
    bplabels = colnames(list$plotmatrix)[1:4]
  } else {
    tmp <- matrix(0, nrow = 2, ncol = 36)
    tmp[1,] <- result[1:36]
    tmp[2,] <- result[37:72]
    bpcol = rep(c('#394398', '#e83020'), 36)
    bplabels = colnames(list$plotmatrix)[1:36]
  }
  
  bp <- barplot(tmp, beside = TRUE, main = list$vcf.name,
                ylim = c(0, ifelse(ymax != 0, 1.5 * ymax, 5)), 
                col = bpcol, axisnames = FALSE,
                ylab = "counts")
  colnames(bp) <- bplabels
  rownames(bp) <- c("Transcribed", "Untranscribed")
  
  legend.list <- legend("topright", legend = c("Transcribed", "Untranscribed"), 
                        fill = c('#394398', '#e83020'), bty = "n", cex = 0.8)
  
  if (pool){
    text(bp, x = (bp[1,]+bp[2,])/2, y = -0.1 * ymax, labels = bplabels, 
         pos = 1, xpd = NA, cex = 1)  
  } else {
    text(bp, x = (bp[1,]+bp[2,])/2, y = - 0.3 * ymax, labels = bplabels, 
         pos = 1, xpd = NA, cex = 0.5, srt = 270)
  }
  
  
  for (type in bplabels){
    p.value <- list$p.values[type]
    if (!is.na(p.value)) {
      # draw asterisk if p.value < 0.05
      if (p.value < 0.05){
        # Get the x coordinates of the line segment to be drawn
        x1 <- bp[1, type]
        x2 <- bp[2, type]
        
        # Get the y coordinates of the line segment to be drawn
        y1 <- y2 <- max(tmp) + 0.03 * max(result)
        
        # Draw the line segment
        segments(x1, y1, x2, y2)
        
        # Draw the asterisk on top of line segment
        x3 <- mean(c(x1, x2))
        y3 <- y1 + max(result) * 0.05
        label <- AssignNumberOfAsterisks(p.value) #p.values actually are q.value
        text(x3, y3, label)
      }
    }
  }
}

#' @keywords internal
CheckSeqContextInIDVCF <- function(vcf, column.to.use) {
  if (0 == nrow(vcf)) return()
  stopifnot(!any(vcf$REF == '-'))
  stopifnot(!any(vcf$ALT == '-'))
  cut.pos <- 1 + (nchar(vcf$column.to.use) - 1) / 2
  stopifnot(cut.pos == round(cut.pos))
  cut.from.ref <- substr(vcf$column.to.use, cut.pos,
                         (cut.pos + nchar(vcf$REF)) - 1)
  error.rows <- which(vcf$REF != cut.from.ref)
  if (any(error.rows > 0)) {
    temp <- tempfile(fileext = ".csv")
    write.csv(vcf[error.rows, ], file = temp)
    stop("Seqence context of reference allele is inconsistent,",
         "see file ", temp)
  }
}

#' Plot transcription strand bias 
#'
#' @param annotated.ID.vcf An ID VCF annotated by
#'   \code{AnnotateIDVCFsWithTransRanges}. It \strong{must} have transcript range
#'   information added.
#'   
#' @param pool TODO Jia Geng
#'   
#' @param damaged.base One of \code{NULL}, \code{"purine"} or
#'   \code{"pyrimidine"}. This function allocates approximately
#'   equal numbers of mutations from \code{damaged.base} into
#'   each of \code{num.of.bins} bin by expression level. E.g.
#'   if \code{damaged.base} is \code{"purine"}, then mutations from
#'   A and G will be allocated in approximately equal numbers to
#'   each expression-level bin. The rationale for the name \code{damaged.base}
#'   is that the direction of strand bias is a result of whether the damage
#'   occurs on a purine or pyrimidine.
#'   If \code{NULL}, the function attempts to infer the \code{damaged.base}
#'   based on mutation counts.
#'   
#' @param ymax Limit for the y axis. If not specified, it defaults to NULL and
#'   the y axis limit equals 1.5 times of the maximum mutation counts in a
#'   specific mutation type.
#'   
#' @return A list whose first element is a logic value indicating whether the
#'   plot is successful. The second element is a named numeric vector containing
#'   the p-values printed on the plot.
#'   
#' @section Note: 
#' The strand bias statistics are Benjamini-Hochberg q-values based on two-sided
#' binomial tests of the mutation counts on the transcribed and untranscribed strands
#' relaitve to the actual abundances of C and T on the transcribed strand. On the
#' plot, asterisks indicate q-values as follows *, \eqn{Q<0.05}; **, \eqn{Q<0.01}; ***,
#'  \eqn{Q<0.001}.
#'  
#' @import data.table
#' 
#' @export
#'
#' @examples 
#' file <- c(system.file("extdata/Strelka-ID-vcf/",
#'                       "Strelka.ID.GRCh37.s1.vcf",
#'                       package = "ICAMSxtra"))
#' ID.vcf <- ReadStrelkaIDVCFs(file)
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   annotated.ID.vcf <- AnnotateIDVCFsWithTransRanges(ID.vcf, ref.genome = "hg19",
#'                                                     trans.ranges = trans.ranges.GRCh37,
#'                                                     vcf.names = "Strelka.ID.GRCh37.s1.vcf")
#'   #' alternatively run below code to skip call to AnnotateIDVCFsWithTransRanges 
#'   #' load(c(system.file("extdata/annotated.ID.vcf.rda",package = "ICAMSxtra")))
#'   ID115_PlotTransBias(annotated.ID.vcf = annotated.ID.vcf[[1]],
#'                       pool = TRUE)
#' }
ID115_PlotTransBias <-
  function(annotated.ID.vcf, pool, damaged.base = NULL, ymax = NULL) { 
    if (is.null(names(annotated.ID.vcf))){
      stop("annotated.ID.vcf does not have 'name' property")
    }
    list <- ID115_StrandBiasGetPlottables(annotated.ID.vcf, damaged.base, pool=pool, vcf.name=name(annotated.ID.vcf))
    ID115_PlotTransBiasInternal(list = list, ymax = ymax)
    return(list(plot.success = TRUE, p.values = list$p.values))
  }

#' Plot transcription strand bias to a PDF file
#' 
#' @inheritParams ID115_PlotTransBias
#' 
#' @param annotated.ID.vcfs TODO Jia Geng
#' 
#' @param file The name of output file.
#'   
#' @importFrom stats glm
#' 
#' @inherit ID115_PlotTransBias return
#' 
#' @inheritSection ID115_PlotTransBias Note
#' 
#' @export
#'
#' @examples
#' library(ICAMS)
#' file <- c(system.file("extdata/Strelka-ID-vcf/",
#'                       "Strelka.ID.GRCh37.s1.vcf",
#'                       package = "ICAMSxtra"))
#' list.of.vcfs <- ReadStrelkaIDVCFs(file)
#' ID.vcf <- list.of.vcfs[[1]]             
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   annotated.ID.vcf <- AnnotateIDVCFsWithTransRanges(ID.vcf, ref.genome = "hg19",
#'                                       trans.ranges = trans.ranges.GRCh37, vcf.names = "Strelka.ID.GRCh37.s1")
#'   ID115_PlotTransBiasToPdf(annotated.ID.vcfs = annotated.ID.vcf,
#'                            file = file.path(tempdir(), "test.pdf"),
#'                            pool = TRUE)
#' }
ID115_PlotTransBiasToPdf <- function(annotated.ID.vcfs, file, pool, damaged.base = NULL) {
  # Setting the width and length for A4 size plotting
  grDevices::pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)
  
  opar <- par(mfrow = c(4, 1), mar = c(8, 5.5, 2, 1), oma = c(1, 1, 2, 1))
  on.exit(par(opar))
  
  plot.type <- ifelse(pool, target_pooled, target)
  n <- length(annotated.ID.vcfs)
  vcf.names <- names(annotated.ID.vcfs)
  results <- list()
  for (i in 1:n){
    list <- ID115_StrandBiasGetPlottables(annotated.ID.vcf = annotated.ID.vcfs[[i]], 
                                          vcf.name = vcf.names[[i]], 
                                          pool = pool, 
                                          damaged.base = damaged.base)
    ymax <- max(list$plotmatrix[, c(plot.type, unname(mapply(revcID115, plot.type)))])
    ID115_PlotTransBiasInternal(list = list, ymax = NULL)
    results[[vcf.names[i]]] <- list(plot.success = TRUE, p.values = list$p.values)
  }
  grDevices::dev.off()
  return(results)
}

#' Add sequence context and transcript information to an in-memory ID VCF
#' 
#' @param ID.vcfs A list of in-memory ID VCF as a \code{data.frame}.
#' 
#' @param ref.genome A \code{ref.genome} argument as described in
#'   \code{\link{ICAMS}}.
#'
#' @param trans.ranges Optional. If \code{ref.genome} specifies one of the
#'   \code{\link{BSgenome}} object 
#'   \enumerate{
#'     \item \code{\link[BSgenome.Hsapiens.1000genomes.hs37d5]{BSgenome.Hsapiens.1000genomes.hs37d5}}
#'     \item \code{\link[BSgenome.Hsapiens.UCSC.hg38]{BSgenome.Hsapiens.UCSC.hg38}}
#'     \item \code{\link[BSgenome.Mmusculus.UCSC.mm10]{BSgenome.Mmusculus.UCSC.mm10}}
#'   }
#'   then the function will infer \code{trans.ranges} automatically. Otherwise,
#'   user will need to provide the necessary \code{trans.ranges}. Please refer to
#'   \code{\link{TranscriptRanges}} for more details.
#'   If \code{is.null(trans.ranges)} do not add transcript range
#'   information.
#'   
#' @param vcf.names list of names of the vcfs
#'
#' @return A list of in-memory ID VCFs each a \code{data.table}. These have been annotated
#'   with the sequence context (column name \code{seq.context}) and with
#'   transcript information in the form of a gene symbol (e.g. \code{"TP53"})
#'   and transcript strand. This information is in the columns
#'   \code{trans.start.pos}, \code{trans.end.pos} , \code{trans.strand},
#'   \code{trans.Ensembl.gene.ID} and \code{trans.gene.symbol} in the output.
#'   These columns are not added if \code{is.null(trans.ranges)}.
#'   
#' @export
#' 
#' @examples 
#' file <- c(system.file("extdata/Strelka-ID-vcf",
#'                       "Strelka.ID.GRCh37.s1.vcf",
#'                       package = "ICAMSxtra"))
#' list.of.vcfs <- ReadStrelkaIDVCFs(file)
#' if (requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5", quietly = TRUE)) {
#'   annotated.ID.vcfs <- AnnotateIDVCFsWithTransRanges(list.of.vcfs, ref.genome = "hg19",
#'                                                     trans.ranges = trans.ranges.GRCh37)
#'   }
AnnotateIDVCFsWithTransRanges <- function(ID.vcfs, ref.genome, trans.ranges = NULL, vcf.names) {
  n <- length(ID.vcfs)
  retval <- list()
  for (i in 1:n){
    ID.vcf <- AnnotateIDVCF(ID.vcfs[[i]], ref.genome=ref.genome)[[1]]
    CheckSeqContextInIDVCF(ID.vcf, "seq.context")
    trans.ranges <- ICAMS:::InferTransRanges(ref.genome, trans.ranges)
    if (!is.null(trans.ranges)) {
      ID.vcf <- ICAMS:::AddTranscript(ID.vcf, trans.ranges)
    }
    retval[[vcf.names[[i]]]] <- (as.data.table(ID.vcf))
  }
  return(retval)
}