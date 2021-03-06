% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TrenaProject.R
\docType{class}
\name{TrenaProject-class}
\alias{TrenaProject-class}
\alias{.TrenaProject}
\alias{TrenaProject}
\title{TrenaProject}
\usage{
TrenaProject(
  projectName,
  supportedGenes,
  genomeName,
  geneInfoTable.path,
  footprintDatabaseHost,
  footprintDatabaseNames,
  footprintDatabasePort = 5432,
  packageDataDirectory,
  quiet
)
}
\arguments{
\item{supportedGenes}{a vector of character strings}

\item{footprintDatabaseHost}{Character string (e.g., "khaleesi.systemsbiology.net")}

\item{footprintDatabaseNames}{Character string (e.g., "hint_brain_20")}

\item{quiet}{A logical indicating whether or not the Trena object should print output}

\item{expressionDirectory}{A string pointing to a collection of RData expression matrices}

\item{genomicRegionsDirectory}{A string pointing to a collection of RData files describing genomic regions of interest}

\item{variantsDirectory}{A string pointing to a collection of RData variant files}

\item{covariatesFile}{the (optional) name of a covariates files}
}
\value{
An object of the TrenaProject class
}
\description{
TrenaProject and its (projected) subclasses provide convenient containers in which to collect
 trena-related aggregation of a gene's (a hierarchy of classes) including expression data,
transcript and variant info, genomic and epigenomic context, trena models and/or the means to create them
}
