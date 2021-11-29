# Functions to find best matches (by cosine similarity) between
# sets of mutational signatures.

#' Compute a matrix of distances / similarities between two sets of signatures.
#' 
#' @param x1 The first set of signatures (a positive matrix in which each column is a signature).
#'     The elements of \code{x1} will be the rows of the output matrix
#' 
#' @param x2 The second set of signatures, similar data type to \code{x1}.
#'     The elements of \code{x2} will be the columns of the output matrix
#'
#' @param method (as for the \code{philentropy::distance}) function.
#' 
#' @return A matrix with dimensions \code{ncol(x1)} X \code{ncol(x2)} with
#'   each element representing the distance or similarity (depending on \code{method})
#'   between the corresponding elments of \code{x1} and \code{x2}
#'
#' @export
#' 
sig_dist_matrix <- function(x1, x2, method = "cosine") {
  mm <- cbind(x1, x2)
  dd <- suppressMessages(
    philentropy::distance(t(mm), method = method, use.row.names = TRUE))
  dd2 <- dd[1:ncol(x1), , drop = FALSE] # Use the rows that represent the elements of x1
  dd3 <- dd2[, -(1:ncol(x1)), drop = FALSE] # Use that columns that represent elements of x2
  return(dd3)
  
}

if (FALSE) {
  sbs <-PCAWG7::signature$genome$SBS96
  match_two_sig_sets(sbs[ , 1:3], sbs[ , 4:5], cutoff = .7)
  match_two_sig_sets(sbs[ , 4:5], sbs[ , 1:3], cutoff = .7)
  ex <- sbs[ , 4:5]
  colnames(ex) <- c("ex1", "ex2")
  TP_FP_FN_avg_sim(extracted.sigs     = ex,
                   ground.truth.sigs = sbs[ , 1:3],
                   similarity.cutoff = .7)
  ex.sigs <- matrix(c(0.2, 0.8, 0.3, 0.7), nrow = 2)
  colnames(ex.sigs) <- c("ex1", "ex2")
  gt.sigs <- matrix(c(0.21, 0.79, 0.19, 0.81), nrow = 2)
  colnames(gt.sigs) <- c("gt1", "gt2")
  sig_dist_matrix(ex.sigs, gt.sigs)
  TP_FP_FN_avg_sim(extracted.sigs     = ex.sigs,
                   ground.truth.sigs = gt.sigs,
                   similarity.cutoff = .985)
  ex.sigs2 <- cbind(ex.sigs, ex3 = c(0.18, 0.82))
  TP_FP_FN_avg_sim(extracted.sigs     = ex.sigs2,
                   ground.truth.sigs = gt.sigs,
                   similarity.cutoff = .985)
  }

#' Return #TP, #FP, #FN, average cosine similarity between extracted
#'  and ground truth signatures.
#' 
#' @details Match signatures in \code{extracted.sigs} to
#'    signatures in \code{ground.truth.sigs} using the function
#'    \code{\link[clue]{solve_LSAP}}, which used the 
#'    "Hungarian" (a.k.a "Kuhn–Munkres") algorithm.
#'    \url{https://en.wikipedia.org/wiki/Hungarian_algorithm},
#'    which optimizes the total cost associated with the links
#'    between nodes.
#'    The function first computes the
#'    all-pairs cosine similarity matrix between the two
#'    sets of signatures, then converts cosine similarities
#'    to cosine distances (including \code{similarity.cutoff})
#'    by subtracting from 1, then
#'    sets distances > the converted cutoff to very large values.
#'    The applies \code{\link[clue]{solve_LSAP}} to the resulting
#'    matrix to compute an optimal matching between 
#'    \code{extracted.sigs} and \code{ground.truth.sigs}.
#' 
#' @param extracted.sigs Mutational signatures discovered by some analysis.
#'    A numerical-matrix-like object with columns as signatures.
#' 
#' @param ground.truth.sigs Ground-truth mutational signatures from
#'    a synthetic data set. A numerical-matrix-like object with columns
#'    as signatures.
#'    
#' @param similarity.cutoff A signature in \code{ground.truth.sigs}
#'    must be matched
#'    by \code{>= similarity.cutoff} by a signature in \code{extracted.sigs}
#'    to be considered detected.

TP_FP_FN_avg_sim <- 
  function(extracted.sigs, ground.truth.sigs, similarity.cutoff = 0.9) {
  tt.and.matrix <- 
    match_two_sig_sets(extracted.sigs, 
                       ground.truth.sigs, 
                       cutoff = similarity.cutoff)
  # browser()
  tt <- tt.and.matrix$table
  TP.sigs <- tt[ , 1]
  FP.sigs <- setdiff(colnames(extracted.sigs), TP.sigs)
  FN.sigs <- setdiff(colnames(ground.truth.sigs), tt[ , 2])
  tt[ , 3] <- 1 - tt[ , 3]
  colnames(tt) <- c("ex.sig", "gt,sig", "sim")
  return(list(TP          = length(TP.sigs),
              FP          = length(FP.sigs),
              FN          = length(FN.sigs),
              avg.cos.sim =  mean(tt[ , 3]),
              table       = tt,
              sim.matrix  = 1 - tt.and.matrix$orig.matrix))
}


#' @export
#' 
match_two_sig_sets <- 
  function(x1, 
           x2, 
           method              = "cosine", 
           convert.sim.to.dist = function(x) {return(1 -x)},
           cutoff              = 0.9) {
    # browser()
    dd  <- sig_dist_matrix(x1, x2, method = method)
    if (!is.null(convert.sim.to.dist)) {
    dd <- convert.sim.to.dist(dd)
    cutoff <- convert.sim.to.dist(cutoff)
  }
  
  return(internal_matches(dd, cutoff))
}

if (FALSE) {
  tmp.dd <- matrix(
    c(0.99, 0.95, 0.89, 0.91),
    nrow = 2)
  rownames(tmp.dd) <- c("x1", "x2")
  colnames(tmp.dd) <- c("g1", "g2")
  tmp.dd <- 1 - tmp.dd
  internal_matches(tmp.dd, cutoff = 0.1)
  
}

# We break this out as a separate function to simplify testing.
internal_matches <- function(dd, cutoff) {

  # Cannot use Inf in the foreign function call clue::solve_LSAP(dd) (below)
  my.inf <- 9e99 
  
  # browser()
  
  was.transformed <- FALSE
  if (nrow(dd) > ncol(dd)) {
    # Add more columns
    dd <- t(dd)
    was.transformed <- TRUE
  }
  
  dd[dd > cutoff] <- my.inf
  
  solution <- clue::solve_LSAP(dd)
  
  cost.one.pair <- function(rowi) {
    # browser()
    colj <- solution[rowi]
    dist <- dd[rowi, colj]
    name1 <- rownames(dd)[rowi]
    name2 <- colnames(dd)[colj]
    return(list(name1 = name1, name2 = name2, dist = dist))
  }
  
  proto.table <- lapply(1:length(solution), FUN = cost.one.pair)
  
  table1 <- matrix(unlist(proto.table), ncol = 3, byrow = TRUE)
  if (was.transformed) {
    table1 <- table1[ , c(2,1,3)]
    dd <- t(dd)
  }
  
  colnames(table1) <- c("x1", "x2", "dist")
  # browser()
  table1 <- data.frame(table1)
  table1$dist <- as.numeric(table1$dist)
    
  table <- table1[table1$dist < my.inf, , drop = FALSE ]

  return(list(table = table, orig.matrix = dd))
      
}

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
#' \code{match2}: a data frame with the row names being signature identifiers
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

MatchSigs2Directions <- function(sigs1, sigs2) {

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

#' An asymmetrical analysis of a set of "ground truth" and "extracted" signatures.
#' 
#' This function is deprecated. You probably want to use
#' \code{\link{TP_FP_FN_avg_sim}}, \code{\link{sig_dist_matrix}},
#' or \code{\link{match_two_sig_sets}}.
#'
#' @param ex.sigs Newly extracted signatures to be compared to \code{gt.sigs}.
#'
#' @param gt.sigs "Ground truth" signatures.
#'
#' @param exposure If \code{NULL}, then match
#'   \code{ex.sigs} against all signatures in \code{gt.sigs}.
#'   Otherwise this should be ground-truth exposures used generate the
#'   synthetic spectra from which \code{ex.sigs} were extracted.
#'   In this case we do not
#'   match to ground-truth signatures that were not in the ground
#'   truth exposure.
#'   
#' @param similarity.cutoff A ground-truth signature must have been
#'   the best match of an extracted signature with a cosine 
#'   similarity \eqn{$ge$} this value to be considered a true positive.
#'   Otherwise we consider the ground-truth signature to be a false
#'   negative.
#'
#' @return A list with the elements 
#' 
#' * \code{averCosSim} The average of the cosine similarities
#'    between each signature in
#'    \code{ex.sigs} and its closest match in \code{gt.sigs}
#'    and the closest match
#'    between each signature in \code{gt.sigs}
#'    and its closest match in \code{ex.sigs}. 
#'    This may not be what you want. Often one wants
#'    the average of the cosine similarities of the true
#'    positives to their matching ground-truth signatures.
#'                 
#' * \code{match1} The match from extracted signatures to ground truth
#'         signatures. Associated with each extracted signature is
#'         a ground truth signature with best cosine similarity.
#'         Ties are resolved arbitrarily.
#'         
#' * \code{match2} The match from ground truth signatures to extracted
#'         signatures. Associated with each extracted signature is
#'         a ground truth signature with best cosine similarity.
#'         Ties are resolved arbitrarily.
#' 
#' * \code{extracted.with.no.best.match}
#' * \code{ground.truth.with.no.best.match}
#' * \code{ex.sigs} Echo input argument
#' * \code{gt.sigs} Echo input argument
#' * \code{gt.mean.cos.sim}
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

MatchSigsAndRelabel <- function(ex.sigs, 
                                gt.sigs, 
                                exposure          = NULL, 
                                similarity.cutoff = 0.9) {

    if (is.null(colnames(ex.sigs))) {
      colnames(ex.sigs) <- paste0("Ex.", 1:ncol(ex.sigs))
    }

    if (!is.null(exposure)) {
      # Remove signatures that are not present in
      # the exposure from which the synthetic data were
      # generated
      exposed.sig.names <- rownames(exposure)[rowSums(exposure) > 0]
      # Make sure we do not have any signatures in exposures that
      # are not in gt.sigs.
      stopifnot(
        setequal(setdiff(exposed.sig.names, colnames(gt.sigs)), c()))
      gt.sigs <- gt.sigs[  , exposed.sig.names]
    }

    gt.sigs.all.sim <- 
      suppressMessages(philentropy::distance(t(gt.sigs), method = "cosine"))
    if (is.null(dim(gt.sigs.all.sim))) {
      if (gt.sigs.all.sim >= similarity.cutoff) {
        warning("The two ground truth signatures have cosine similarity >= ",
                similarity.cutoff)
      }
    } else {
        # browser()
        gt.sigs.all.sim[lower.tri(gt.sigs.all.sim, diag = TRUE)] <- 0
        if (any(gt.sigs.all.sim >= similarity.cutoff)) {
          rownames(gt.sigs.all.sim) <- colnames(gt.sigs)
          colnames(gt.sigs.all.sim) <- colnames(gt.sigs)
          warning(
            "Some ground truth signatures have cosine similarity >= ",
            similarity.cutoff)
          print(gt.sigs.all.sim)
        }
    }

    sim <- MatchSigs2Directions(ex.sigs, gt.sigs)

    true.match1 <- sim$match1
    # This is the match from extracted signatures to ground-truth signature
    true.match1 <- true.match1[true.match1$sim >= similarity.cutoff, ]
  
    true.match2 <- sim$match2 
    # This is the match from ground-truth signature to extracted signature
    true.match2 <- true.match2[true.match2$sim >= similarity.cutoff, ]

    sim$extracted.with.no.best.match <-
      setdiff(colnames(ex.sigs), true.match2$to)
    # These are extracted signatures that are not the best match
    # of any ground truth signature.

    sim$ground.truth.with.no.best.match <-
      setdiff(colnames(gt.sigs), true.match1$to)
    # These are the ground truth signatures that are not the 
    # best match of any extracted signatures.
    
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
