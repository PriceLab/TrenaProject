#' @import trena
#' @importFrom DBI dbConnect dbListTables dbGetQuery dbListConnections dbDisconnect
#' @importFrom RPostgreSQL dbConnect dbListTables dbGetQuery dbListConnections dbDisconnect
#' @import RPostgreSQL
#' @importFrom methods new
#'
#' @title TrenaProject
#------------------------------------------------------------------------------------------------------------------------
#' @name TrenaProject-class
#' @rdname TrenaProject-class
#' @aliases TrenaProject
#'
## @import methods

.TrenaProject <- setClass ("TrenaProject",
                        representation = representation(
                           projectName="character",
                           supportedGenes="character",
                           geneInfoTable="data.frame",
                           genomeName="character",
                           footprintDatabaseHost="character",
                           footprintDatabaseNames="character",
                           footprintDatabasePort="numeric",
                           packageDataDirectory="character",
                           state="environment",
                           quiet="logical"
                           )
                         )


#------------------------------------------------------------------------------------------------------------------------
setGeneric('getProjectName',                signature='obj', function(obj) standardGeneric('getProjectName'))
setGeneric('getSupportedGenes',             signature='obj', function(obj) standardGeneric('getSupportedGenes'))
setGeneric('setTargetGene',                 signature='obj', function(obj, targetGene, curatedGenesOnly=FALSE) standardGeneric('setTargetGene'))
setGeneric('getGenome',                     signature='obj', function(obj) standardGeneric('getGenome'))
setGeneric('getTargetGene',                 signature='obj', function(obj) standardGeneric('getTargetGene'))
setGeneric('getGeneInfoTable',              signature='obj', function(obj) standardGeneric('getGeneInfoTable'))
setGeneric('getFootprintDatabaseHost',      signature='obj', function(obj) standardGeneric ('getFootprintDatabaseHost'))
setGeneric('getFootprintDatabasePort',      signature='obj', function(obj) standardGeneric ('getFootprintDatabasePort'))
setGeneric('getFootprintDatabaseNames',     signature='obj', function(obj) standardGeneric ('getFootprintDatabaseNames'))
setGeneric('getTranscriptsTable',           signature='obj', function(obj, targetGene=NA, all=FALSE) standardGeneric ('getTranscriptsTable'))
setGeneric('getDataDirectory',              signature='obj', function(obj) standardGeneric ('getDataDirectory'))
#setGeneric('getExpressionDirectory',        signature='obj', function(obj) standardGeneric ('getExpressionDirectory'))
setGeneric('getExpressionMatrixNames',      signature='obj', function(obj) standardGeneric ('getExpressionMatrixNames'))
setGeneric('getExpressionMatrix',           signature='obj', function(obj, matrixName) standardGeneric ('getExpressionMatrix'))
setGeneric('getVariantDatasetNames',        signature='obj', function(obj) standardGeneric ('getVariantDatasetNames'))
setGeneric('getVariantDataset',             signature='obj', function(obj, datasetName) standardGeneric ('getVariantDataset'))
setGeneric('getCovariateDatasetNames',      signature='obj', function(obj) standardGeneric ('getCovariateDatasetNames'))
setGeneric('getCovariateDataset',           signature='obj', function(obj, datasetName) standardGeneric ('getCovariateDataset'))
setGeneric('getGenomicRegionsDatasetNames', signature='obj', function(obj) standardGeneric ('getGenomicRegionsDatasetNames'))
setGeneric('getGenomicRegionsDataset',      signature='obj', function(obj, datasetName) standardGeneric ('getGenomicRegionsDataset'))
#' @export
setGeneric('getEnhancers',                  signature='obj', function(obj, targetGene=NA, tissues="all", maxSize=10000) standardGeneric ('getEnhancers'))
#' @export
setGeneric('getEnhancerTissues',            signature='obj', function(obj) standardGeneric ('getEnhancerTissues'))
#' @export
setGeneric('getGeneRegulatoryRegions',      signature='obj', function(obj, targetGene=NA, tissues="all",
                                                                      geneHancerMissing.promoter.upstream=0,
                                                                      geneHancerMissing.promoter.downstream=0,
                                                                      geneHancerSupplemental.promoter.upstream=0,
                                                                      geneHancerSupplemental.promoter.downstream=0)
                                                               standardGeneric ('getGeneRegulatoryRegions'))
#' @export
setGeneric('getClassicalGenePromoter',      signature='obj', function(obj, targetGene=NA, upstream=5000, downstream=5000)
               standardGeneric ('getClassicalGenePromoter'))
#' @export
setGeneric('getEncodeDHS',                  signature='obj', function(obj, targetGene=NA) standardGeneric ('getEncodeDHS'))
#' @export
setGeneric('getChipSeq',                    signature='obj', function(obj, chrom, start, end, tfs=NA) standardGeneric ('getChipSeq'))
#' @export
setGeneric('getGeneRegion',                 signature='obj', function(obj, flankingPercent=0) standardGeneric ('getGeneRegion'))
#' @export
setGeneric('getGeneEnhancersRegion',        signature='obj', function(obj, flankingPercent=0) standardGeneric ('getGeneEnhancersRegion'))
#' @export
setGeneric('getProximalPromoter',           signature='obj', function(obj, upstream, downstream) standardGeneric ('getProximalPromoter'))
setGeneric('recognizedGene',                signature='obj', function(obj, geneName) standardGeneric ('recognizedGene'))
setGeneric('getAllTranscriptionFactors',    signature='obj', function(obj, source) standardGeneric ('getAllTranscriptionFactors'))
#------------------------------------------------------------------------------------------------------------------------
#' Define an object of class Trena
#'
#' @description
#' TrenaProject and its (projected) subclasses provide convenient containers in which to collect
#'  trena-related aggregation of a gene's (a hierarchy of classes) including expression data,
#' transcript and variant info, genomic and epigenomic context, trena models and/or the means to create them
#'
#' @rdname TrenaProject-class
#'
#' @param supportedGenes a vector of character strings
#' @param footprintDatabaseHost Character string (e.g., "khaleesi.systemsbiology.net")
#' @param footprintDatabaseNames Character string (e.g., "hint_brain_20")
#' @param expressionDirectory A string pointing to a collection of RData expression matrices
#' @param genomicRegionsDirectory A string pointing to a collection of RData files describing genomic regions of interest
#' @param variantsDirectory A string pointing to a collection of RData variant files
#' @param covariatesFile  the (optional) name of a covariates files
#' @param quiet A logical indicating whether or not the Trena object should print output
#'
#' @return An object of the TrenaProject class
#'
#' @export
#'
#'
TrenaProject <- function(projectName,
                         supportedGenes,
                         genomeName,
                         geneInfoTable.path,
                         footprintDatabaseHost,
                         footprintDatabaseNames,
                         footprintDatabasePort=5432,
                         packageDataDirectory,
                         quiet)
{

   state <- new.env(parent=emptyenv())
   state$targetGene <- NULL
   state$tbl.transcripts <- NULL

   if(is.na(geneInfoTable.path)){
      tbl.geneInfo <- data.frame()
   } else {
      stopifnot(file.exists(geneInfoTable.path))
      tbl.geneInfo <- get(load(geneInfoTable.path))
     }

     # TODO: (pshannon 14jun2019)
     # TODO: the location of many of the standard data directories are simply calculated.
     # TODO: should be no need to specify them in explicit constructor arguments.
     # TODO: can the identidy and directory of the ultimate derived class be learned
     # TODO: here at runtime?

   .TrenaProject(projectName=projectName,
                 supportedGenes=supportedGenes,
                 genomeName=genomeName,
                 geneInfoTable=tbl.geneInfo,
                 footprintDatabaseHost=footprintDatabaseHost,
                 footprintDatabaseNames=footprintDatabaseNames,
                 footprintDatabasePort=footprintDatabasePort,
                 packageDataDirectory=packageDataDirectory,
                 state=state,
                 quiet=quiet)


} # ctor
#------------------------------------------------------------------------------------------------------------------------
#' get a summary of the objects
#'
#' @rdname show
#' @aliases show
#'
#' @param obj An object of class TrenaProject
#'
#' @export

setMethod('show', 'TrenaProject',

    function(object) {
       cat(sprintf ('--- TrenaProject'), '\n', sep='')
       cat(sprintf("projectName: '%s'", getProjectName(object)), "\n", sep="")
       cat(sprintf("genomeName: '%s'", getGenome(object)), "\n", sep="")
       cat(sprintf("reg regions db host: %s", getFootprintDatabaseHost(object)), "\n", sep='')
       cat(sprintf("reg regions db port: %d", getFootprintDatabasePort(object)), "\n", sep='')
       cat(sprintf("reg regions db names: %s", paste(getFootprintDatabaseNames(object), collapse=", ")), "\n", sep='')
       cat(sprintf("data directory: %s", getDataDirectory(object)), "\n", sep='')
       })

#------------------------------------------------------------------------------------------------------------------------
#' get the project name - to be set by derived classes
#'
#' @rdname getProjectName
#' @aliases getProjectName
#'
#' @param obj An object of class TrenaProject
#'
#' @export

setMethod('getProjectName', 'TrenaProject',

   function(obj) {
      obj@projectName
      })

#------------------------------------------------------------------------------------------------------------------------
#' get the genome name - to be set by derived classes
#'
#' @rdname getGenome
#' @aliases getGenome
#'
#' @param obj An object of class TrenaProject
#'
#' @export

setMethod('getGenome', 'TrenaProject',

   function(obj) {
      obj@genomeName
      })

#------------------------------------------------------------------------------------------------------------------------
#' get the list of genes supported in this project
#'
#' @rdname getSupportedGenes
#' @aliases getSupportedGenes
#'
#' @param obj An object of class TrenaProject
#'
#' @export

setMethod('getSupportedGenes', 'TrenaProject',

   function(obj) {
      obj@supportedGenes
      })

#------------------------------------------------------------------------------------------------------------------------
#' Set a single gene for analysis
#'
#' @rdname setTargetGene
#' @aliases setTargetGene
#'
#' @param obj An object of class TrenaProject
#' @param targetGene a characteor string, the name of the gene
#'
#' @export

setMethod('setTargetGene', 'TrenaProject',

   function(obj, targetGene, curatedGenesOnly=FALSE) {
      if(curatedGenesOnly){
         if(!all(is.na(getSupportedGenes(obj))))
            stopifnot(targetGene %in% getSupportedGenes(obj))
         }
      obj@state$targetGene <- targetGene
      xyz <- "about to set tbl.transcripts"
      targetGene.regex <- sprintf("^%s$", targetGene)
      index <- grep(toupper(targetGene.regex), toupper(obj@geneInfoTable$geneSymbol))
      if(length(index) == 0)
         return(FALSE)
      tbl.tmp <- obj@geneInfoTable[index,]
      obj@state$tbl.transcripts  <- tbl.tmp
      return(TRUE)
      })

#------------------------------------------------------------------------------------------------------------------------
#' get the single gene currently designated for analysis
#'
#' @rdname getTargetGene
#' @aliases getTargetGene
#'
#' @param obj An object of class TrenaProject
#'
#' @export

setMethod('getTargetGene', 'TrenaProject',

   function(obj) {
      obj@state$targetGene
      })

#------------------------------------------------------------------------------------------------------------------------
#' get the host on which the footprints Postgres database is running
#'
#' @rdname getFootprintDatabaseHost
#' @aliases getFootprintDatabaseHost
#'
#' @param obj An object of class TrenaProject
#'
#' @export

setMethod('getFootprintDatabaseHost', 'TrenaProject',

   function(obj) {
      obj@footprintDatabaseHost
      })

#------------------------------------------------------------------------------------------------------------------------
#' get the port  on which the footprints Postgres database is running
#'
#' @rdname getFootprintDatabasePort
#' @aliases getFootprintDatabasePort
#'
#' @param obj An object of class TrenaProject
#'
#' @export

setMethod('getFootprintDatabasePort', 'TrenaProject',

   function(obj) {
      obj@footprintDatabasePort
      })

#------------------------------------------------------------------------------------------------------------------------
#' get the names of the database tables holding footprints
#'
#' @rdname getFootprintDatabaseNames
#' @aliases getFootprintDatabaseNames
#'
#' @param obj An object of class TrenaProject
#'
#' @export

setMethod('getFootprintDatabaseNames', 'TrenaProject',

   function(obj) {
      obj@footprintDatabaseNames
      })

#------------------------------------------------------------------------------------------------------------------------
#' Get the names of the expression directory, wherein serialized expression matrices will be found
#'
#' @rdname getDataDirectory
#' @aliases getDataDirectory
#'
#' @param obj An object of class TrenaProject
#'
#' @export

setMethod('getDataDirectory',  'TrenaProject',

   function(obj) {
      obj@packageDataDirectory
      })

#------------------------------------------------------------------------------------------------------------------------
#' Get the names of the expression matrices - their names with directory and .RData suffix stripped out
#'
#' @rdname getExpressionMatrixNames
#' @aliases getExpressionMatrixNames
#'
#' @param obj An object of class TrenaProject
#'
#' @export

setMethod('getExpressionMatrixNames',  'TrenaProject',

   function(obj) {
     directory <- file.path(getDataDirectory(obj), "expression")
     if(dir.exists(directory)){
        filenames <- grep(".RData", list.files(directory), fixed=TRUE, value=TRUE)
        datasetNames <- sub(".RData", "", filenames)
        return(datasetNames)
        }
      return(list())
      })

#------------------------------------------------------------------------------------------------------------------------
#' Get the a specifc expression matrix
#'
#' @rdname getExpressionMatrix
#' @aliases getExpressionMatrix
#'
#' @param obj An object of class TrenaProject
#' @param matrixName A numeric matrix
#'
#' @export

setMethod('getExpressionMatrix',  'TrenaProject',

    function(obj, matrixName){
       directory <- file.path(getDataDirectory(obj), "expression")
       filename <- sprintf("%s.RData", matrixName)
       full.path <- file.path(directory, filename)
       stopifnot(file.exists(full.path))
       mtx <- NULL
       eval(parse(text=paste("mtx <- ", load(full.path))))
       invisible(mtx)
       })

#------------------------------------------------------------------------------------------------------------------------
#' List the RData files in the variants directory
#'
#' @rdname getVariantDatasetNames
#' @aliases getVariantDatasetNames
#'
#' @param obj An object of class TrenaProject
#'
#' @export

setMethod('getVariantDatasetNames', 'TrenaProject',

      function(obj){
         variantsDirectory <- file.path(obj@packageDataDirectory, "variants")
         if(!file.exists(variantsDirectory))
            return(list())
         filenames <- sub(".RData", "", list.files(variantsDirectory), fixed=TRUE)
         return(filenames)
         })

#------------------------------------------------------------------------------------------------------------------------
#' Get the specified variants table
#'
#' @rdname getVariantDataset
#' @aliases getVariantDataset
#'
#' @param obj An object of class TrenaProject
#' @param datasetName character string, the variant dataset of interest
#'
#' @export

setMethod('getVariantDataset', 'TrenaProject',

    function(obj, datasetName){
        variantsDirectory <- file.path(obj@packageDataDirectory, "variants")
        if(!file.exists(variantsDirectory))
           return(NULL)
        file.name <- sprintf("%s.RData", file.path(variantsDirectory, datasetName))
        tbl <- NULL
        eval(parse(text=sprintf("tbl <- %s", load(file.name))))
        tbl
        })

#------------------------------------------------------------------------------------------------------------------------
#' Get a list of datasets each describing regions of interest, typically regulatory regions, open chromatin, etc
#'
#' @rdname getRegionDatasetNames
#' @aliases getRegionDatasetNames
#'
#' @param obj An object of class TrenaProject
#'
#' @export

setMethod('getGenomicRegionsDatasetNames', 'TrenaProject',

    function(obj){
       directory <- file.path(getDataDirectory(obj), "genomicRegions")
       if(dir.exists(directory)){
          filenames <- grep(".RData", list.files(directory), fixed=TRUE, value=TRUE)
          filenames <- sub(".RData", "", filenames)
          return(filenames)
          }
        return(list())
        })

#------------------------------------------------------------------------------------------------------------------------
#' Get a list of datasets each describing regions of interest, typically regulatory regions, open chromatin, etc
#'
#' @rdname getRegionDatasetNames
#' @aliases getRegionDatasetNames
#'
#' @param obj An object of class TrenaProject
#'
#' @export

setMethod('getGenomicRegionsDataset', 'TrenaProject',

    function(obj, datasetName){
       directory <- file.path(getDataDirectory(obj), "genomicRegions")
       filename <- sprintf("%s.RData", datasetName)
       full.path <- file.path(directory, filename)
       stopifnot(file.exists(full.path))
       tbl <- NULL
       eval(parse(text=paste("tbl <- ", load(full.path))))
       invisible(tbl)
       })

#------------------------------------------------------------------------------------------------------------------------
#' List the RData files in the covariate directory
#'
#' @rdname getCovariateDatasetNames
#' @aliases getCovariateDatasetNames
#'
#' @param obj An object of class TrenaProject
#'
#' @export

setMethod('getCovariateDatasetNames', 'TrenaProject',

      function(obj){
         covariatesDirectory <- file.path(obj@packageDataDirectory, "covariate")
         if(!file.exists(covariatesDirectory))
            return(list())
         filenames <- sub(".RData", "", list.files(covariatesDirectory), fixed=TRUE)
         return(filenames)
         })

#------------------------------------------------------------------------------------------------------------------------
#' Get the specified covariates table
#'
#' @rdname getCovariateDataset
#' @aliases getCovariateDataset
#'
#' @param obj An object of class TrenaProject
#' @param datasetName character string, the covariate dataset of interest
#'
#' @export

setMethod('getCovariateDataset', 'TrenaProject',

    function(obj, datasetName){
        covariatesDirectory <- file.path(obj@packageDataDirectory, "covariates")
        if(!file.exists(covariatesDirectory))
           return(NULL)
        file.name <- sprintf("%s.RData", file.path(covariatesDirectory, datasetName))
        tbl <- NULL
        eval(parse(text=sprintf("tbl <- %s", load(file.name))))
        tbl
        })

#------------------------------------------------------------------------------------------------------------------------
#' return the data.frame with gene ids, chromosome, tss of the primary transcript, and strand for all genes in the project
#'
#' @rdname getGeneInfoTable
#' @aliases getGeneInfoTable
#'
#' @param obj An object of class TrenaProject
#'
#' @return a data.frame with gene names as rownames
#'
#' @export

setMethod('getGeneInfoTable',  'TrenaProject',

   function(obj){
      return(obj@geneInfoTable)
      })

#------------------------------------------------------------------------------------------------------------------------
#' Do we have expression data for the suggested gene? genomic and epigenetic information?
#'
#' @rdname recognizedGene
#' @aliases recognizedGene
#'
#' @param obj An object of class TrenaProject
#' @param geneName A character string in the same protocol as the project's expression matrices
#'
#' @return a chrom.loc (chrom:start-end) string
#'
#' @export

setMethod('recognizedGene',  'TrenaProject',

   function(obj, geneName){
      geneName.regex <- sprintf("^%s$", geneName)
      index <- grep(toupper(geneName.regex), toupper(obj@geneInfoTable$geneSymbol))
      return(length(index) > 0)
      })

#------------------------------------------------------------------------------------------------------------------------
#' Get the transcript info for the gene, just the first row
#'
#' @rdname getTranscriptsTable
#' @aliases getTranscriptsTable
#'
#' @param obj An object of class TrenaProject
#'
#' @export

setMethod('getTranscriptsTable',  'TrenaProject',

   function(obj, targetGene=NA_character_, all=FALSE) {

      tbl <- getGeneInfoTable(obj)
      if(all)
         return(tbl)

      if(nrow(tbl) == 0){
         message("no geneInfoTable for this project, returning empty data.frame")
         return(data.frame())
         }

      if(!is.na(targetGene)){
         return(subset(tbl, geneSymbol==targetGene))
         }

      targetGene <- getTargetGene(obj)
      if(is.null(targetGene)){
         message("no targetGene set for this project, none supplied as argument to this function")
         return(data.frame())
         }

      return(subset(tbl, geneSymbol==targetGene))
      })

#------------------------------------------------------------------------------------------------------------------------
#' Get the chromosomal region surrounding the current targetGene, with a flanking percentage added up and downstream
#'
#' @rdname getGeneRegion
#' @aliases getGeneRegion
#'
#' @param obj An object of class TrenaProject
#' @param flankingPercent a numeric percentage of the gene's total span
#'
#' @return a chrom.loc (chrom:start-end) string
#' @export

setMethod('getGeneRegion',  'TrenaProject',

     function(obj, flankingPercent=0){
        tbl.transcripts <- getTranscriptsTable(obj)[1,]  # currently always nrow of 1
        chrom <- tbl.transcripts$chr
        start <- tbl.transcripts$start
        end   <- tbl.transcripts$end
        span <- 1 + end - start
        flank <- round(span * (flankingPercent/100))
        chromLocString <- sprintf("%s:%d-%d", chrom, start - flank, end + flank)
        list(chrom=chrom, start=start-flank, end=end+flank, chromLocString=chromLocString)
        })

#------------------------------------------------------------------------------------------------------------------------
#' Get the chromosomal region surrounding the current targetGene, with a flanking percentage added up and downstream
#'
#' @rdname getProximalPromoter
#' @aliases getProximalPromoter
#'
#' @param obj An object of class TrenaProject
#' @param upstream numeric, number of bases before (strand-aware) the TSS
#' @param downstream numeric, number of bases after (strand-aware) the TSS
#'
#' @return a chrom.loc (chrom:start-end) string
#' @export

setMethod('getProximalPromoter',  'TrenaProject',

     function(obj, upstream, downstream){
        tbl.geneInfo <- getTranscriptsTable(obj)[1,]  # currently always nrow of 1
        chrom <- tbl.geneInfo$chrom
        tss <- tbl.geneInfo$tss
        start.loc <- tss - upstream
        end.loc   <- tss + downstream
        if(tbl.geneInfo$strand == -1){
          start.loc <- tss - downstream
          end.loc   <- tss + upstream
          }
       loc.string <- sprintf("%s:%d-%d", chrom, start.loc, end.loc)
       list(chrom=chrom, start=start.loc, end=end.loc, chromLocString=loc.string)
       })

#------------------------------------------------------------------------------------------------------------------------
#' get all genes annotated by GO to
#'
#' @rdname getAllTranscriptionFactors
#' @aliases getAllTranscriptionFactors
#'
#' @param obj An object of class TrenaProject
#' @param source A character string e.g., "MotifDb", "Gene Ontology"
#'
#' @export

setMethod('getAllTranscriptionFactors', 'TrenaProject',

   function(obj, source) {
      c()
      })

#------------------------------------------------------------------------------------------------------------------------
