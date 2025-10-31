##Download Dada2 package if necessary
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ShortRead", version = "3.17")
##Load DaDa2Package
library(dada2); packageVersion("dada2")
#Version 1.28.0

#dplyer is necessary for the error model later on
library(dplyr); packageVersion("dplyr")


#Removing primersites: This is necessary because PCR can induced mistakes/substitutions in primer site
trimLeft = c(19, 20)

#Define path, unzipped fastq files
path <- "~/data/volume_2/16S"
list.files(path)

# Define paths for forward and reverse files seperately. Always check if the format of new data is the same.
fnFs <- sort(list.files(path, pattern="_R1_001", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001", full.names = TRUE))
# Extract sample names, I prefer the sample type over the eNummer, but they should be combined with which plate they were on so Do1T1L_TM033 would

sample.names <- sapply(strsplit(basename(fnFs), "_"), function(x) paste(x[3], x[4], sep = "_"))

## This is to check the quality of the sequences. Takes about a minute for two files. 
# Q score is logaritmic (Q=-10*log10(P)), where P is probability of a wrong base.
# Q40 is 0.0001 (0.01%) wrong bases, Q30 0.001 (0.1%). Q30 is a common border to use for quality.
plotQualityProfile(fnRs[1:2])
plotQualityProfile(fnFs[1:2])
# Plots are a little bit weird, because Baseclear uses binned quality scores for Novaseq. But Q score seems to be oke.

# IF NEW DATA: Place filtered files in filtered/subdirectory. The files will be empty
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
#iF FILTERED DIRECTORY ALREADY EXISTS: The first time i'm doing this I am not managing to do the entire pipeline in one day

#Time to filter and trim based on Q scores.
# MaxEE stands for max errors. So maximum of 2 errors in sequence (both F and R). This is based on an addition of Q scores pre base.
# Trimleft is added to trim the primersides, since they are not trimmed yet.
# Parada and Appril primers for 16s-V4, about 300 bp, so there is a lot of overlap, no problem in big cut in R primer.

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, trimLeft = c(19,20)) # On Windows set multithread=FALSE
head(out)
# This takes about 20/25min for 32 files on my Dell laptop. First step that takes considerable time.


# Dada2's standard error model does not work properly with the binned quality data
# Four different solutions are offered on https://github.com/benjjneb/dada2/issues/1307 by users
# Dada2 does not yet have an official preference, the fourth seems to work best, but has a hard to fix error that also occurs in my testdata.
# Winny tested the other three models, 1 worked best for her. I'll use that one for now.
# But the github should be looked at before every new run, maybe an even better solution will present itself.
# Pay attention: library(dplyr); packageVersion("dplyr") has to be loaded for this error model
loessErrfun_mod1 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# check what this looks like
errF_1 <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod1,
  verbose = TRUE
)

# For Reverse
errR_1 <- learnErrors(
  filtRs,
  multithread = FALSE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod1,
  verbose = TRUE
)

## The error check is the most time consuming step and I am not able to do it on my laptop.

#Time to check te result:
plotErrors(errF_1, nominalQ=TRUE)
plotErrors(errR_1, nominalQ=TRUE)

#Sample inference
### Extensions: By default, the dada function processes each sample 
#independently. However, pooling information across samples can increase
#sensitivity to sequence variants that may be present at very low frequencies
#in multiple samples. The dada2 package offers two types of pooling. 
#dada(..., pool=TRUE) performs standard pooled processing,
#in which all samples are pooled together for sample inference. 
#dada(..., pool="pseudo") performs pseudo-pooling, in which samples
#are processed independently after sharing information between samples,
#approximating pooled sample inference in linear time.
dadaFs <- dada(filtFs, err=errF_1, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR_1, multithread = TRUE)
#With a check
dadaFs[[1]]

## MERGE DATA
#We now merge the forward and reverse reads together to obtain the 
#full denoised sequences. Merging is performed by aligning the denoised 
#forward reads with the reverse-complement of the corresponding denoised
#reverse reads, and then constructing the merged “contig” sequences.
#By default, merged sequences are only output if the forward and reverse reads
#overlap by at least 12 bases, and are identical to each other in the 
#overlap region (but these conditions can be changed via function arguments).
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

###We can now construct an amplicon sequence variant table (ASV) table, 
# a higher-resolution version of the OTU table produced by traditional methods.
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#Dim returens the dimentions of the matrix, to give an indication of differences within.

#This inspects the distribution of sequence lenghts.
#sequencetable(nchar(getSequences(seqtab)))
# Maybe make a nice curve out of this? Expect a bell curve:
sequence_lengths <- table(nchar(getSequences(seqtab)))
sequence_lengths_df <- data.frame(length = as.numeric(names(sequence_lengths)), count = as.numeric(sequence_lengths))
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)
ggplot(sequence_lengths_df, aes(x = length, y = count)) +
  geom_bar(stat = "identity") +
  labs(x = "Sequence Length", y = "Count", title = "Distribution of Sequence Lengths")

### Based on the output of the sequence length, it might be necessary to remove
# seqs that are too long or short (<240, >265 for instance). 
# For now I decided to keep them, because it is a small fraction, but you never
# know what it might entail.
#CODE IF FOR REMOVAL:
#  seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 240:265]). 


#REMOVING CHMIMERAS. This is easier if the ASVs are there than before.
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#Most of your reads should remain after chimera removal 
#(it is not uncommon for a majority of sequence variants to be removed though).
#If most of your reads were removed as chimeric, upstream processing may need 
#to be revisited. In almost all cases this is caused by primer sequences with 
#ambiguous nucleotides that were not removed prior to beginning the DADA2 
#pipeline. To check percentage of removed reads:
sum(seqtab.nochim)/sum(seqtab)

#Last check before assigning taxonomy!
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
### This is a great place to do a last sanity check. 
#Outside of filtering, there should no step in which a majority of reads 
#are lost. If a majority of reads failed to merge, you may need to revisit the
#truncLen parameter used in the filtering step and make sure that the 
#truncated reads span your amplicon. If a majority of reads were removed as
#chimeric, you may need to revisit the removal of primers, as the ambiguous 
#nucleotides in unremoved primers interfere with chimera identification.

##TIME TO ASSIGN TAXONOMY USING SILVA DATABASE!!!!!!!
# First a file with taxonomic info should be downloaded and put in the right directory:
# File can be found here: https://zenodo.org/record/4587955
# When working with short sequences, first we only assign the ASVs to the genus level.
# After this, we can go to species level with the add.species function.
taxa <- assignTaxonomy(seqtab.nochim, "~/data/volume_2/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
# Takes <1h
#add species. Same step as before, download the file first.
# Only looks at exact matches.
taxa <- addSpecies(taxa, "~/data/volume_2/silva_species_assignment_v138.1.fa.gz")

# Up here uses a bayesian method 
