% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ID115_fns.R
\name{as.catalog.for.ID115}
\alias{as.catalog.for.ID115}
\title{Create a catalog from a \code{matrix}, \code{data.frame}, or \code{vector}}
\usage{
as.catalog.for.ID115(
  object,
  ref.genome = NULL,
  region = "unknown",
  catalog.type = "counts",
  abundance = NULL,
  infer.rownames = FALSE
)
}
\arguments{
\item{object}{A numeric \code{matrix}, numeric \code{data.frame},
or \code{vector}.
If a \code{vector}, converted to a 1-column \code{matrix}
with rownames taken from the element names of the \code{vector}
and with column name \code{"Unknown"}.
If argument \code{infer.rownames}
is \code{FALSE} than this argument must have
rownames to denote the mutation types. See \code{\link{CatalogRowOrder}}
for more details.}

\item{ref.genome}{A \code{ref.genome} argument as described in
\code{\link{ICAMS}}.}

\item{region}{A character string designating a region, one of
\code{genome}, \code{transcript}, \code{exome}, \code{unknown};
see \code{\link{ICAMS}}.}

\item{catalog.type}{One of "counts", "density", "counts.signature",
"density.signature".}

\item{abundance}{If \code{NULL}, then
inferred if \code{ref.genome}
is one of
the reference genomes known to ICAMS and \code{region}
is not \code{unknown}. See \code{\link{ICAMS}}.
The argument \code{abundance} should
contain the counts of different source sequences for mutations
in the same format as the numeric vectors in \code{\link{all.abundance}}.}

\item{infer.rownames}{If \code{TRUE}, and \code{object} has no
rownames, then assume the rows of \code{object} are in the
correct order and add the rownames implied by the number of rows
in \code{object} (e.g. rownames for SBS 192 if there are 192 rows).
If \code{TRUE}, \strong{be sure the order of rows is correct.}}
}
\value{
A catalog as described in \code{\link{ICAMS}}.
}
\description{
Create a catalog from a \code{matrix}, \code{data.frame}, or \code{vector}
}
