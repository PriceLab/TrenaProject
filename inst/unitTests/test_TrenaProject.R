# test_TrenaProject
#------------------------------------------------------------------------------------------------------------------------
library(TrenaProject)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tp")){

   projectName <- "base class test"
   genomeName <- "unicorn99"

   genes <- c("fu", "bar")
   geneInfoTable.path <- NA_character_

   footprintDatabaseHost <- "fake.host.net"
   footprintDatabaseNames <- c("db1", "db2")

   expressionDirectory <- NA_character_
   genomicRegionsDirectory <- system.file(package="TrenaProject", "extdata", "genomicRegions")
   variantsDirectory <- NA_character_
   covariatesFile <- NA_character_

   tp <- TrenaProject(projectName=projectName,
                      genomeName=genomeName,
                      supportedGenes=genes,
                      geneInfoTable.path=geneInfoTable.path,
                      footprintDatabaseHost=footprintDatabaseHost,
                      footprintDatabaseNames=footprintDatabaseNames,
                      expressionDirectory=expressionDirectory,
                      genomicRegionsDirectory=genomicRegionsDirectory,
                      variantsDirectory=variantsDirectory,
                      covariatesFile=covariatesFile,
                      quiet=TRUE)
   } # creating trenaProj for use in multiple functions below

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   printf("--- test_ctor")

   checkTrue("TrenaProject" %in% is(tp))
   checkEquals(getSupportedGenes(tp), genes)
   checkEquals(getFootprintDatabaseHost(tp), footprintDatabaseHost)
   checkEquals(getFootprintDatabaseNames(tp), footprintDatabaseNames)

   printf("--- testing get/setTargetGene")
   setTargetGene(tp, genes[1])
   checkEquals(getTargetGene(tp), genes[1])

   printf("--- getting transcript info for %s", genes[1])

   suppressMessages(tbl.transcripts <- getTranscriptsTable(tp))
   checkTrue(nrow(tbl.transcripts) == 0)

   names <- getGenomicRegionsDatasetNames(tp)
   checkTrue(all(c("ATAC-seq-erythropoiesis-d04_rep1", "ATAC-seq-erythropoiesis-d12_rep2") %in% names))
   tbl.atac <- getGenomicRegionsDataset(tp, "ATAC-seq-erythropoiesis-d12_rep2")
   checkEquals(dim(tbl.atac), c(6, 4))

} # test_ctor
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
