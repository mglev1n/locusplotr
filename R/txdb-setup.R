# # Create TxDb objects for hg19 and hg38
# # retrieve transcript lengths
# suppressPackageStartupMessages({
#   invisible(lapply(c("ggbio", "biovizBase", "data.table",
#                      "TxDb.Hsapiens.UCSC.hg19.knownGene",
#                      "org.Hs.eg.db"),
#                    require, character.only = TRUE))})
# txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#
# txlen <- transcriptLengths(txdb, with.utr5_len=TRUE, with.utr3_len=TRUE)
# setDT(txlen)
# txlen$len <- rowSums(as.matrix(txlen[, .(tx_len, utr5_len, utr3_len)]))
# setkey(txlen, gene_id, len, tx_id)
#
# # filter longesttranscript by gene_id
# ltx <- txlen[!is.na(gene_id)][, tail(.SD,1), by=gene_id]$tx_id
#
# # filter txdb object
# txb <- as.list(txdb)
# txb$transcripts <- txb$transcripts[txb$transcripts$tx_id %in% ltx, ]
# txb$splicings <- txb$splicings[txb$splicings$tx_id %in% ltx,]
# txb$genes <- txb$genes[txb$genes$tx_id %in% ltx,]
# txb_hg19 <- do.call(makeTxDb, txb)
#
# # Save
# AnnotationDbi::saveDb(txb_hg19, "Data/txb_hg19.sqlite")
#
# ## Setup txdb for hg38
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#
# txlen <- transcriptLengths(txdb, with.utr5_len=TRUE, with.utr3_len=TRUE)
# setDT(txlen)
# txlen$len <- rowSums(as.matrix(txlen[, .(tx_len, utr5_len, utr3_len)]))
# setkey(txlen, gene_id, len, tx_id)
#
# # filter longesttranscript by gene_id
# ltx <- txlen[!is.na(gene_id)][, tail(.SD,1), by=gene_id]$tx_id
#
# # filter txdb object
# txb <- as.list(txdb)
# txb$transcripts <- txb$transcripts[txb$transcripts$tx_id %in% ltx, ]
# txb$splicings <- txb$splicings[txb$splicings$tx_id %in% ltx,]
# txb$genes <- txb$genes[txb$genes$tx_id %in% ltx,]
# txb_hg38 <- do.call(makeTxDb, txb)
#
# # Save
# AnnotationDbi::saveDb(txb_hg38, "Data/txb_hg38.sqlite")
