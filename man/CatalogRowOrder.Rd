% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CatalogRowOrder}
\alias{CatalogRowOrder}
\alias{catalog.row.order}
\title{Standard order of row names in a catalog}
\format{
A list of character vectors indicating the standard
orders of row names in catalogs.
}
\usage{
catalog.row.order
}
\description{
This data is designed for those
who need to create their own catalogs from formats not
supported by this package. The rownames denote the mutation
types.  For example, for SBS96 catalogs, the rowname
AGAT represents a mutation from AGA > ATA.
}
\note{
In ID (small insertion and deletion) catalogs, deletion repeat sizes
range from 0 to 5+, but for plotting and end-user documentation deletion
repeat sizes range from 1 to 6+. In ID83 catalogs, deletion repeat sizes
range from 0 to 5.
}
\examples{
catalog.row.order$ID115
# There are altogether 115 row names to denote the mutation types
# in ID115 catalog.
# The difference from the .$ID class in \code{\link{ICAMS::catalog.row.order}} is that 
# single base nonhomopolymer indels have trinucleotide context added to them in the format
# INS/DEL:C/T:1:0_PF where P is the preceding base and F is the following base.
}
\keyword{datasets}
