% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match.sigs.R
\name{Match1Sig}
\alias{Match1Sig}
\title{Find signatures in \code{other.sigs} with the highest cosine similarity to \code{query.sig}.}
\usage{
Match1Sig(query.sig, other.sigs)
}
\arguments{
\item{query.sig}{A single signature.}

\item{other.sigs}{Matrix with each column being one signature.}
}
\value{
The maximum similarity between \code{query.sig}
and any signature in
\code{other.sigs}; the name of the single element
in the vector is the name
of a signature with the maximum similarity.
}
\description{
Find signatures in \code{other.sigs} with the highest cosine similarity to \code{query.sig}.
}
