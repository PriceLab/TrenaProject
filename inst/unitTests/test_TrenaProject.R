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
   variantsDirectory <- NA_character_
   covariatesFile <- NA_character_

   tp <- TrenaProject(projectName=projectName,
                      genomeName=genomeName,
                      supportedGenes=genes,
                      geneInfoTable.path=geneInfoTable.path,
                      footprintDatabaseHost=footprintDatabaseHost,
                      footprintDatabaseNames=footprintDatabaseNames,
                      expressionDirectory=expressionDirectory,
                      variantsDirectory=variantsDirectory,
                      covariatesFile=covariatesFile,
                      quiet=TRUE)
   } # creating trenaProj for use in multiple functions below

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()
   test_ctor_withFootprintDatabasePortSpecified()
   test_getEnhancers()
   test_getPrimaryTranscriptInfo()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   printf("--- test_ctor")

   checkEquals(getSupportedGenes(trenaProj), genes)
   checkEquals(getFootprintDatabaseHost(trenaProj), footprintDatabaseHost)
   checkEquals(getFootprintDatabaseNames(trenaProj), footprintDatabaseNames)

   printf("--- testing get/setTargetGene")
   #checkTrue(is.null(getTargetGene(trenaProj)))
   setTargetGene(trenaProj, genes[1])
   checkEquals(getTargetGene(trenaProj), genes[1])

   printf("--- getting transcript info for %s", genes[1])

   tbl.transcripts <- getTranscriptsTable(trenaProj)
   checkTrue(nrow(tbl.transcripts) == 1)

   printf("--- testing get/setTargetGene")
   #checkTrue(is.null(getTargetGene(trenaProj)))
   setTargetGene(trenaProj, "PIGF")           # a placental gene, not in the IGAP project
   checkEquals(getTargetGene(trenaProj), "PIGF")

   printf("--- getting transcript info for %s", "PIGF")
   tbl.transcripts <- getTranscriptsTable(trenaProj)
   checkTrue(nrow(tbl.transcripts) == 1)

     # return to TREM2, whose coordinates we check below
   setTargetGene(trenaProj, genes[1])

   checkEquals(getExpressionMatrixNames(trenaProj), c("dummyExpressionSet_1", "dummyExpressionSet_2"))

   checkTrue(is.matrix(getExpressionMatrix(trenaProj, "dummyExpressionSet_1")))
   checkTrue(is.matrix(getExpressionMatrix(trenaProj, "dummyExpressionSet_2")))

   expected <- c("someGene.region.vcf", "tbl.snp.gwas.minimal")
   file.list <- getVariantDatasetNames(trenaProj)
   checkTrue(all(expected %in% file.list))
   #checkTrue(file.exists("someGene.region.vcf"))

     # most variant files - other than vcfs - are serialized into .RData files,
     # with that suffix stripped off for human readers (in a presumed Shiny UI)

   #checkTrue(file.exists(sprintf("%s.RData", file.list[["tbl.snp.gwas.minimal"]])))

   tbl.covariates <- getCovariatesTable(trenaProj)

   tbl.enhancers <- getEnhancers(trenaProj)
   checkEquals(colnames(tbl.enhancers), c("chrom", "start", "end", "type", "combinedScore", "geneSymbol"))

   checkTrue(nrow(tbl.enhancers) >= 5)
   checkEquals(unique(tbl.enhancers$geneSymbol), getTargetGene(trenaProj))

   tbl.dhs <- getEncodeDHS(trenaProj)
   checkEquals(colnames(tbl.dhs), c("chrom", "chromStart", "chromEnd", "count", "score"))
      # these open chromatin regions bear no necessary relation to the targetGene
      # instead, they span the region of the targetGene-associated enhancers.  check that
   loc.min <- min(tbl.enhancers$start)
   loc.max <- max(tbl.enhancers$end)
   chromosome <- unique(tbl.enhancers$chrom)
   checkEquals(length(chromosome), 1)
   checkTrue(all(tbl.dhs$chromStart >= loc.min))
   checkTrue(all(tbl.dhs$chromStart <= loc.max))
   checkTrue(all(tbl.dhs$chromEnd >= loc.min))
   checkTrue(all(tbl.dhs$chromEnd <= loc.max))
   checkTrue(all(tbl.dhs$chrom == chromosome))

   tbl.chipseq <- getChipSeq(trenaProj, chrom=chromosome, start=loc.min, end=loc.max, tfs=NA)
   checkTrue(nrow(tbl.chipseq) > 2000)

   checkEquals(colnames(tbl.chipseq), c("chrom", "start", "endpos", "tf", "name", "strand", "peakStart", "peakEnd"))

   checkEquals(getGeneRegion(trenaProj)$chromLocString,      "chr6:41158506-41163186")
   checkEquals(getGeneRegion(trenaProj, flankingPercent=20)$chromLocString, "chr6:41157570-41164122")

   checkEquals(getGeneEnhancersRegion(trenaProj)$chromLocString,                     "chr6:41154324-41210533")
   checkEquals(getGeneEnhancersRegion(trenaProj, flankingPercent=10)$chromLocString, "chr6:41148703-41216154")

   vf <- getVariantDatasetNames(trenaProj)

   checkTrue(nrow(getGeneInfoTable(trenaProj)) > 15000)
   checkTrue(ncol(getGeneInfoTable(trenaProj)) >= 10)

   checkEquals(getFootprintDatabasePort(trenaProj), 5432)


   checkTrue(!recognizedGene(trenaProj, "bogusGene"))      # only genes in the tbl.geneInfo are recognized

} # test_ctor
#------------------------------------------------------------------------------------------------------------------------
test_ctor_withFootprintDatabasePortSpecified <- function()
{
   printf("--- test_ctor_withFootprintDatabasePortSpecified")

   trenaProj <- TrenaProject(supportedGenes=genes,
                             geneInfoTable.path=geneInfoTable.path,
                             footprintDatabaseHost=footprintDatabaseHost,
                             footprintDatabasePort=5433,
                             footprintDatabaseNames=footprintDatabaseNames,
                             expressionDirectory=expressionDirectory,
                             variantsDirectory=variantsDirectory,
                             covariatesFile=covariatesFile,
                             quiet=TRUE)

   checkEquals(getFootprintDatabasePort(trenaProj), 5433)

} # test_ctor_withFootprintDatabsePortSpecified
#------------------------------------------------------------------------------------------------------------------------
test_getEnhancers <- function()
{
   printf("--- test_getEnhancers")

   setTargetGene(trenaProj, "TREM2")
   tbl.trem2 <- getEnhancers(trenaProj)
   checkTrue(all(tbl.trem2$geneSymbol == "TREM2"))

   tbl.mef2c <- getEnhancers(trenaProj, "MEF2C")
   checkTrue(all(tbl.mef2c$geneSymbol == "MEF2C"))

   tbl.trem2.again <- getEnhancers(trenaProj, "TREM2")
   checkEquals(tbl.trem2, tbl.trem2.again)

   tbl.bogus <- getEnhancers(trenaProj, "bogus99")
   checkEquals(nrow(tbl.bogus), 0)

} # test_getEnhancers
#------------------------------------------------------------------------------------------------------------------------
test_getPrimaryTranscriptInfo <- function()
{
   printf("--- test_getPrimaryTranscriptInfo")

   checkEquals(getPrimaryTranscriptInfo(trenaProj, "CRH")$tss, 66178725)

   setTargetGene(trenaProj, "TREM2")
   checkEquals(getPrimaryTranscriptInfo(trenaProj)$tss, 41163176)

} # test_getPrimaryTranscriptInfo
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()