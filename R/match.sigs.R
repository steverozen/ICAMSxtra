# Functions to find best matches (by cosine similarity) between
# sets of mutational signatures.

#' Find signatures in \code{other.sigs} with the highest cosine similarity to \code{query.sig}.
#'
#' @param query.sig A single signature.
#'
#' @param other.sigs Matrix with each column being one signature.
#'
#' @return The maximum similarity between \code{query.sig}
#' and any signature in
#' \code{other.sigs}; the name of the single element
#' in the vector is the name
#' of a signature with the maximum similarity.
#'
#' @export
#'
Match1Sig <- function(query.sig, other.sigs) {
  sims <-
    apply(other.sigs,
          MARGIN = 2,
          FUN = function(other.sig) {
            return(lsa::cosine(query.sig, other.sig))
            })
  max.sim <- max(sims)
  max.name <- names(sims)[sims == max.sim]
  if(length(max.sim) > 1)
    message("There is more than one signature in other.sigs with the highest cosine similarity to query.sig.\n")
  names(max.sim) <- max.name[1] # There might be mulitple signatures
                                # with the maximum similarity.
  return(max.sim)
}

#' Find the closest match in \code{other.sigs} for each signature in \code{query.sigs}
#'
#' @param query.sigs A signature matrix; signatures for which to find the
#'             closest match in \code{other.sigs}. The colnames are used
#'               as the identifiers of the signatures.
#'
#'
#' @param other.sigs A signature matrix; find the closest matches to
#'  a signature in this matrix.
#'  The colnames are used as the identifiers of the signatures.
#'
#'
#' @return A list with one element for each signature in
#'   \code{query.sigs}. The names
#'   of the list elements are the colnames of \code{query.sigs}.
#'   Each list element
#'   is a vector of length 1, and the name of the vector element is
#'   the name of the closest matching signature in \code{other.sigs},
#'   and the value
#'   is the cosine similarity between the given signature in
#'   \code{query.sigs} and
#'   the matching signature in \code{other.sigs}.
#'
#' @export
#'
MatchSigs1Direction <- function(query.sigs, other.sigs) {

  stopifnot(!is.null(colnames(query.sigs)))

  Match1SigInternal <- function(query.sig.name) {
    query.sig <- query.sigs[ , query.sig.name]
    return(Match1Sig(query.sig, other.sigs))
  }

  ret <-
    lapply(X = colnames(query.sigs),
           FUN = Match1SigInternal)
  names(ret) <- colnames(query.sigs)
  return(ret)
}

#' Bidirectional closest similarities between two sets of signatures.
#'
#' @param sigs1 Matrix of signatures; colnames are used as signature
#'   identifiers, and the colnames in \code{sigs1} should be distinguishable from
#'   those in \code{sigs2}.
#'
#' @param sigs2 Matrix of signatures; colnames are used as signature identifiers.
#'
#' @return A list with the elements:
#'
#' \code{averCosSim}: the average of the cosine similarities between each signature in
#' \code{sigs1} and its closest match in \code{sigs2} and the closest match
#' between each signature in \code{sigs2} and its closest match in \code{sigs1}.
#'
#' \code{match1}: a data frame with rownames being signature identifiers from
#' \code{sigs1}, the signature identifier of the closest match in \code{sigs1}
#' in the 1st column, and the cosine similarity between them in the 2nd column.
#'
#' \code{match2}: a data frame with the rownames being signature identifiers
#' from \code{sigs2}, the signature identifier of the closest match in
#' \code{sigs1} in the 1st column, and the cosine similarity between them in the
#' 2nd column.
#'
#' \code{match1} and \code{match2} might not have the same number of rows.
#'
#' @export
#'
#' @examples
#' seta <- matrix(c(1, 3,   4, 1, 2, 4), ncol = 2)
#' setb <- matrix(c(1, 3.1, 4, 5, 1, 1, 1, 2.8, 4), ncol = 3)
#' colnames(seta) <- c("A.1", "A.2")
#' colnames(setb) <- c("B.1", "B.2", "B.3")
#' MatchSigs2Directions(seta, setb)
#'
MatchSigs2Directions <- function(sigs1, sigs2) {
  # TODO(Steve): match1 and match2 can be simplified

  if (is.null(colnames(sigs1))) {
    colnames(sigs1) <- paste0("Ex.", 1:ncol(sigs1))
  }

  match1 <- MatchSigs1Direction(sigs1, sigs2)
  match2 <- MatchSigs1Direction(sigs2, sigs1)
  averCosSim <-
    (sum(unlist(match1)) + sum(unlist(match2))) /
    (length(match1) + length(match2))

  table1 <-
    data.frame(from=names(match1),
               to=unlist(lapply(match1, names)),
               sim=unlist(match1),
               stringsAsFactors = FALSE)
  rownames(table1) <- table1$from
  table1 <-table1[ , -1, drop = FALSE]

  table2 <-
    data.frame(from=names(match2),
               to=unlist(lapply(match2, names)),
               sim=unlist(match2),
               stringsAsFactors = FALSE)
  rownames(table2) <- table2$from
  table2 <-table2[ , -1, drop = FALSE]

  return(list(averCosSim=averCosSim, match1=table1, match2=table2))
}

#' Get the numerical parts of identifiers.
#'
#' @param s A character vector.
#'
#' @return A vector, each element of which is the integer
#' corresponding to the first string of digits of an element of s.
#'
#' @details Not very sophisticated.
#'
#' @export
#'
#' @examples
#' x<- c("SBS22", "SBS2", "SBS7b", "SBS7a")
#' NumFromId(x)
#' x[order(NumFromId(x))]
#'
NumFromId<- function(s) {
  return(
    as.numeric(
      sub("[^0123456789]*(\\d+).*", "\\1", s, perl = TRUE)))
}

#' A somewhat asymmetrical analysis of a set of "ground truth" and "extracted" signatures.
#'
#' @param ex.sigs Newly extracted signatures to be compared to gt.sigs
#            (actually, this is more general).
#
#' @param gt.sigs "Ground truth" signatures.
#'
#' @param exposure If \code{NULL}, then match
#'   \code{ex.sigs} against all signatures in \code{gt.sigs}.
#'   Otherwise this should be ground-truth exposures used generate the
#'   synthetic spectra from which \code{ex.sigs} were extracted.
#'   In this case we do not
#'   match to ground-truth signatures to that were not in the ground
#'   truth exposure.
#'
#' @return A list with the elements \code{averCosSim}, \code{match1},
#' \code{match2} as for \code{\link{MatchSigs2Directions}},
#' with \code{match1} being matches for the the extracted signatures
#' (\code{ex.sigs}) and \code{match2} being the matches for
#' the ground truth signatures (\code{gt.sigs}). The return list
#' also echos the input arguments \code{ex.sigs} and \code{gt.sigs}.
#'
#' @export
#'
#' @examples
#' gt.sigs <- matrix(c(1, 3,   4, 1, 2, 4), ncol = 2)
#' ex.sigs <- matrix(c(1, 3.1, 4, 5, 1, 1, 1, 2.8, 4), ncol = 3)
#' colnames(gt.sigs) <- c("gt.1", "gt.2")
#' colnames(ex.sigs) <- c("ex.1", "ex.2", "ex.3")
#' tout <- MatchSigsAndRelabel(gt.sigs = gt.sigs, ex.sigs = ex.sigs)
#' tout
#'

MatchSigsAndRelabel <- function(ex.sigs, gt.sigs, exposure = NULL) {

    if (is.null(colnames(ex.sigs))) {
      colnames(ex.sigs) <- paste0("Ex.", 1:ncol(ex.sigs))
    }

    if (!is.null(exposure)) {
      # Remove signatures that are not present in
      # the exposure from which the synthetic data were
      # generated
      exposed.sig.names <- rownames(exposure)[rowSums(exposure) > 0]
      # Make sure we do not have an signatures in exposures that
      # are not in gt.sigs.
      stopifnot(
        setequal(setdiff(exposed.sig.names, colnames(gt.sigs)), c()))
      gt.sigs <- gt.sigs[  , exposed.sig.names]
    }

    sim <- MatchSigs2Directions(ex.sigs, gt.sigs)

    ## Software-reported signatures with a best cosine similarity lower than 0.90
    ## is not considered an "extracted signature", it will rather be regarded as
    ## an artefact.
    true.match1 <- sim$match1
    true.match1 <- true.match1[true.match1$sim >= 0.9, ]
    true.match2 <- sim$match2
    true.match2 <- true.match2[true.match2$sim >= 0.9, ]

    sim$extracted.with.no.best.match <-
      setdiff(colnames(ex.sigs), true.match2$to)

    sim$ground.truth.with.no.best.match <-
      setdiff(colnames(gt.sigs), true.match1$to)
    # TODO(Steve) Review documentation / explanation. Note that
    # e.g. SBS29 might have a best match (BI_COMPOSITE_SBS18_P)
    # but no BI signatures has SBS29 as its best match

    # TODO(Steve) Document the complexity below; mostly it deals
    # with setting up plotting that is easy(?) to interpret.
    labels <- character(ncol(ex.sigs))
    names(labels) <- colnames(ex.sigs)
    nums <- NumFromId(sim$match1$to)
    reordered.ex <- colnames(ex.sigs)[order(nums)]
    ex.sigs.x <- ex.sigs[ , order(nums),drop = FALSE]
    bestmatch.id <- sim$match1[reordered.ex, "to"]
    bestmatch.sim <- sim$match1[reordered.ex, "sim"]
    bestmatch.sim <- round(bestmatch.sim, digits=4)
    init.labels <-
      paste0(reordered.ex, " (", bestmatch.id, " ", bestmatch.sim, ")")
    names(init.labels) <- reordered.ex
    laggards <- setdiff(rownames(sim$match2), bestmatch.id)
    # Falling back to a loop here:
    for (lag in laggards) {
      my.ex.id  <- sim$match2[lag, "to"]
      my.ex.sim <- round(sim$match2[lag, "sim"], digits = 4)
      init.labels[my.ex.id] <-
        paste0(init.labels[my.ex.id],
               " (", lag, " ", my.ex.sim, ")")
    }
    colnames(ex.sigs.x) <- init.labels

    sim$ex.sigs <- ex.sigs.x
    sim$gt.sigs <- gt.sigs

    # Calculate cosine similarity between all extracted signatures, and each of
    # the ground-truth signatures. For example, first calculate the cosine
    # similarity between ground-truth SBS5 and all extracted signatures most
    # similar to SBS5; then calculate the cosine similarity between SBS1 and all
    # extracted signatures most similar to SBS1.
    x2gt <- sim$match1
    # sim$match1 maps each extracted signature to the ground truth
    # signature it is most similar to and to the cosine similarity to this signature.
    sim$cosSim <- list()
    for (gtSigName in rownames(sim$match2)) {

      # In x2gt, Find all extracted signatures similar to
      # ground-truth signature gtSigName.
      values <- x2gt[which(x2gt[,1] == gtSigName), 2]

      if(is.nan(mean(values))) {
        # None of the extracted signatures were most similar to gtSigName.
        # Therefore, we go to sim$match2 and look for the extracted
        # signature that gtSigName is most similar to.
        #
        # This can show how an extracted signature is
        # similar to multiple ground-truth signatures.
        gt2x <- sim$match2
        value <- gt2x[gtSigName, 2]
        sim$cosSim[[gtSigName]] <- value
      } else{
        # There are some extracted signatures most similar to gtSigName.
        sim$cosSim[[gtSigName]] <- mean(values)
      }
    }

    sim$gt.mean.cos.sim <- sim$cosSim
    sim$cosSim <- NULL

    invisible(sim)
  }
