library(biomaRt)
library(org.Mm.eg.db)

#wash7p.unprocessed_pseudogene <- "ENSG00000227232"
#subset(tbl.geneInfo, ensg=="ENSG00000227232")
#           ensg chrom start   end   tss strand geneSymbol entrez    appris   tsl      transcript                   type
#ENSG00000227232  chr1 14404 29570 29570     -1     WASH7P 653635 C_missing tslNA ENST00000488147 unprocessed_pseudogene


ens.genes <- select(org.Mm.eg.db, keys=keys(org.Mm.egENSEMBL), keytype="ENTREZID", columns="ENSEMBL")$ENSEMBL
length(ens.genes) # 68912

deleters <- which(is.na(ens.genes))
length(deleters)
if(length(deleters) > 0)
   ens.genes <- ens.genes[-deleters]
length(ens.genes)
length(unique(ens.genes))
ens.genes <- unique(ens.genes)

ensembl.mm10 <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

coi <- c("ensembl_gene_id", "entrezgene", "chromosome_name", "transcription_start_site",
         "strand", "mgi_symbol", "transcript_biotype", "gene_biotype", "ensembl_transcript_id",
         "transcript_appris", "transcript_tsl", "transcript_gencode_basic")

key <- "ensembl_gene_id"
printf("--------- calling biomart for ensg/symbol/chrom/tss mapping")

tbl.geneInfo <- getBM(attributes=coi, filters=key,
                      values=ens.genes,
                      mart=ensembl.mm10)
dim(tbl.geneInfo)
save(tbl.geneInfo, file="mm10.biomart.table.raw.RData")

appris <- tbl.geneInfo$transcript_appris
appris <- sub("alternative", "B_alternative", appris)
appris <- sub("principal"  , "A_principal",   appris)
missing <- which(appris == "")
appris[missing] <- "C_missing"

tsl <- tbl.geneInfo$transcript_tsl
missing <- which(tsl == "")
tsl[missing] <- "z_missing"

tbl.geneInfo$appris <- appris
tbl.geneInfo$tsl <- tsl

new.order <- with(tbl.geneInfo, order(ensembl_gene_id, appris, tsl))

tbl <- tbl.geneInfo[new.order,]

dups <- which(duplicated(tbl$ensembl_gene_id))
length(dups)
tbl <- tbl[-dups,]
dim(tbl)
symbols.missing <- which(nchar(tbl$mgi_symbol) == 0)
length(symbols.missing)  # 1507
tbl$mgi_symbol[symbols.missing] <- tbl$ensembl_gene_id[symbols.missing]

#  [1] "ensembl_gene_id"
#  [2] "entrezgene"
#  [3] "chromosome_name"
#  [4] "transcription_start_site"
#  [5] "strand"
#  [6] "hgnc_symbol"
#  [7] "transcript_biotype"
#  [8] "gene_biotype"
#  [9] "ensembl_transcript_id"
# [10] "transcript_appris"
# [11] "transcript_tsl"
# [12] "transcript_gencode_basic"
# [13] "appris"
# [14] "tsl"

colnames.oi <- c(1, 3, 4, 5, 6, 2, 13, 14, 9, 7)
tbl.1 <- tbl[, colnames.oi]
colnames(tbl.1) <- c("ensg", "chrom", "tss", "strand", "geneSymbol", "entrez", "appris", "tsl", "transcript", "type")
rownames(tbl.1) <- NULL
tbl.1$chrom <- paste("chr", tbl.1$chrom, sep="")

tbl.geneInfo <- tbl.1

   # keep the work done so far, for possible restart.  final file is written below

save(tbl.geneInfo, file="mm10.biomart.table.cooked.RData")

annotated.genes <- tbl.geneInfo$ensg
coi.2 <- c("ensembl_gene_id", "transcript_start", "transcript_end")

# ensg.bug <- "ENSG00000067601"   # no start/end

tbl.geneBounds <- getBM(attributes=coi.2, filters=key,
                        values=annotated.genes,
                        mart=ensembl.mm10)
dim(tbl.geneBounds)


endpoints <- function(ensg.id) {
   tbl <- subset(tbl.geneBounds, ensembl_gene_id==ensg.id)
   data.frame(ensg=ensg.id,
              start=min(tbl$transcript_start, na.rm=TRUE),
              end=max(tbl$transcript_end, na.rm=TRUE),
              stringsAsFactors=FALSE)
   }

ensg.with.bounds <- unique(tbl.geneBounds$ensembl_gene_id)
length(ensg.with.bounds)

x <- lapply(ensg.with.bounds, endpoints)
tbl.bounds <- do.call(rbind, x)
dim(tbl.bounds)
dim(tbl.geneInfo)

tbl.wide <- merge(tbl.geneInfo, tbl.bounds, by="ensg")
preferred.column.order <- c("ensg",
                            "chrom",
                            "start",
                            "end",
                            "tss",
                            "strand",
                            "geneSymbol",
                            "entrez",
                            "appris",
                            "tsl",
                            "transcript",
                            "type"
                            )
tbl.geneInfo <- tbl.wide[, preferred.column.order]
dim(tbl.geneInfo)
"Abca1" %in% tbl.geneInfo$geneSymbol

save(tbl.geneInfo, file="../../extdata/geneInfoTable_mm10.RData")
