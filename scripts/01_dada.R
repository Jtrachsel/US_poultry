library(dada2); packageVersion("dada2")
path <- './raw_data'
list.files(path)


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=50) # On Windows set multithread=FALSE
head(out)



errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool=TRUE)


mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


seqtab <- makeSequenceTable(mergers)
dim(seqtab)



# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
hist(nchar(getSequences(seqtab)))

# HEREHERE
# maybe this
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]

table(nchar(getSequences(seqtab2)))
hist(nchar(colnames(seqtab2)))

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

table(nchar(getSequences(seqtab.nochim)))

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)



# Original dada2 classifier
# replaced with DECIPHER methods

# taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
# taxa <- addSpecies(taxa, "~/tax/silva_species_assignment_v132.fa.gz")

#
# taxa.print <- taxa # Removing sequence rownames for display only
# rownames(taxa.print) <- NULL
# head(taxa.print)


#######

library(DECIPHER); packageVersion("DECIPHER")

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

taxa <- taxid

### all these unclassified ASVs are too long for the amplicon we tried to generate
colnames(taxa)
unclass <- taxa[is.na(taxa[,"domain"]),]
nrow(unclass) # 28 ASVs at Kingdom level...
rownames(unclass)
unclass_ASVs <- DNAStringSet(x = rownames(unclass))

writeXStringSet(unclass_ASVs, filepath = 'unclass_ASVs.fasta')

#
hist(width(tst))
hist(track[,'nonchim'])

### write out results ###
library(tidyverse)

# track for QC stuff
track <-
  track %>%
  as.data.frame() %>%
  rownames_to_column(var='sample_ID') %>%
  write_tsv('./processed_data/QC_stats_reads_per_sample.tsv')

# taxa for taxonomy

taxa <-
  taxa %>%
  as.data.frame() %>%
  rownames_to_column(var='ASV_seq') %>%
  write_tsv('./processed_data/ASV_taxonomy.tsv')

# seqtab.nochim = ASV count table

seqtab.nochim <-
  seqtab.nochim %>%
  as.data.frame() %>%
  rownames_to_column(var='sample_ID') %>%
  write_tsv('./processed_data/ASV_counts.tsv')

