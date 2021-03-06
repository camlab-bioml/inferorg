% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inferorg.R
\name{summarize_organism}
\alias{summarize_organism}
\title{Summarize for a given organism}
\usage{
summarize_organism(df, gene_vec)
}
\value{
A vector of length 3 giving the proportion of elements
in \code{gene_vec} that match the symbols, ensembl IDs and entrez
IDs for the organism that \code{gene_vec} corresponds to.
}
\description{
Summarize for a given organism
}
\keyword{internal}
