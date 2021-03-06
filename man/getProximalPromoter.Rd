% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TrenaProject.R
\name{getProximalPromoter,TrenaProject-method}
\alias{getProximalPromoter,TrenaProject-method}
\alias{getProximalPromoter}
\title{Get the chromosomal region surrounding the current targetGene, with a flanking percentage added up and downstream}
\usage{
\S4method{getProximalPromoter}{TrenaProject}(obj, upstream, downstream)
}
\arguments{
\item{obj}{An object of class TrenaProject}

\item{upstream}{numeric, number of bases before (strand-aware) the TSS}

\item{downstream}{numeric, number of bases after (strand-aware) the TSS}
}
\value{
a chrom.loc (chrom:start-end) string
}
\description{
Get the chromosomal region surrounding the current targetGene, with a flanking percentage added up and downstream
}
