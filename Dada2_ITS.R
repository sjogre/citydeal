#Required libraries
library(dada2)
library(ShortRead)
library(Biostrings)
library(magrittr)
library(dplyr)
library(MASS)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocParallel")

cutadapt <- "C:/Users/sjoerd.gremmen/Documents/R-scripts/cutadapt"
# Dada2 also has functions that you could inport with source("path_to_dada2/FunctionsDADA2.R")
# I have chosen to put the functions in the pipeline, because they are not timeconsuming
# And it gives the option to tweak if necessary.

#set path to fastq files (input)
path <- "C:/Users/sjoerd.gremmen/Documents/ITS1"  ## CHANGE ME to the directory containing the fastq files.
list.files(path)

#Define forward and reverse files in list.
fnFs <- sort(list.files(path, pattern = "_R1_", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_", full.names = TRUE))

#Identify primers: Here is where de ITS pipeline differs from the 16s pipeline.
#Because the length of ITS is very variable, we can't just cut of the primer length/truncate.
# If ITS <250bp, both primers will be in the forward and reverse seq somewhere, and need to be
# actively searched for.

#Define primers
# Because we use multiple forward primers we use the cutadapt function, instead of what is done in the
# Dada2 tutorial. Alternatively we could have used degenerate primers, but Winny had a code
# ready to go.
FWD1 <- "CTAGACTCGTCATCGATGAAGAACGCAG"
FWD2 <- "CTAGACTCGTCAACGATGAAGAACGCAG"
FWD3 <- "CTAGACTCGTCACCGATGAAGAACGCAG"
FWD4 <- "CTAGACTCGTCATCGATGAAGAACGTAG"
FWD5 <- "CTAGACTCGTCATCGATGAAGAACGTGG"
FWDn <- "CTAGACTCGTCANCGATGAAGAACGYRG"
REV <- "TTCCTSCGCTTATTGATATGC"


## Here we make a function to define all orientations of all primers
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
# All the primers (in all orders) are saved in .orients
# To view: //.orients
FWD1.orients <- allOrients(FWD1)
FWD2.orients <- allOrients(FWD2)
FWD3.orients <- allOrients(FWD3)
FWD4.orients <- allOrients(FWD4)
FWD5.orients <- allOrients(FWD5)
REV.orients <- allOrients(REV)

#Found out that cutadapt also works with degenerate primers, this one is all primers combined.
# Should save some computing time (5x)
FWDn.orients <- allOrients(FWDn)

# FILTER 1
# Here we make a new path where we can save the filtered reads.
fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
# We filter the reads for other characters than A/G/T/C
# When running on local Windows pc, set multithread to FALSE.
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE)
#sometimes an error is displayed in the file.write step. This is because of the 
#long filenames (so positive control is the one that is to long)
#Check

##Function for primerhits:
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
# Here we make a table of the primers found and how many
rbind(FWD1.ForwardReads = sapply(FWD1.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD1.ReverseReads = sapply(FWD1.orients, primerHits, fn = fnRs.filtN[[1]]),
      FWD2.ForwardReads = sapply(FWD2.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD2.ReverseReads = sapply(FWD2.orients, primerHits, fn = fnRs.filtN[[1]]),
      FWD3.ForwardReads = sapply(FWD3.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD3.ReverseReads = sapply(FWD3.orients, primerHits, fn = fnRs.filtN[[1]]),
      FWD4.ForwardReads = sapply(FWD4.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD4.ReverseReads = sapply(FWD4.orients, primerHits, fn = fnRs.filtN[[1]]),
      FWD5.ForwardReads = sapply(FWD5.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD5.ReverseReads = sapply(FWD5.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# PRIMER REMOVAL
# Make a new path where the reads can be saved after primer removal
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWDn.RC <- dada2:::rc(FWDn)
REV.RC <- dada2::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
# Only FWD1 is used, error rate default of 0.1 takes out all forward reads
R1.flags <- paste("-g", FWD1, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-g", REV, "-a", FW1.RC) 
# Only FWD2 is used, the default error rate is 0.1 which means that the 
# difference between all different primers are covered

# Removing the primers from the reads
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i],  # output files
                             fnFs[i], fnRs[i], #input files
                             "-m",1)) # filters empty reads
}

# Make a table again to see if it is all removed.
rbind(FWD1.ForwardReads = sapply(FWD1.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD1.ReverseReads = sapply(FWD1.orients, primerHits, fn = fnRs.cut[[1]]),
      FWD2.ForwardReads = sapply(FWD2.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD2.ReverseReads = sapply(FWD2.orients, primerHits, fn = fnRs.cut[[1]]),
      FWD3.ForwardReads = sapply(FWD3.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD3.ReverseReads = sapply(FWD3.orients, primerHits, fn = fnRs.cut[[1]]),
      FWD4.ForwardReads = sapply(FWD4.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD4.ReverseReads = sapply(FWD4.orients, primerHits, fn = fnRs.cut[[1]]),
      FWD5.ForwardReads = sapply(FWD5.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD5.ReverseReads = sapply(FWD5.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
# After running the command above, I saw that FWD3 was still commonly seen in certain sequences.
# I ran cutadapt for a second time, but now with FWD3 as the sequence.Lets see if it works

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_", full.names = TRUE))


# Now you could look at some Quality plots:
plotQualityProfile(cutFs[4:5])
plotQualityProfile(fnRs.cut[1:2])
## For Revers the quality drops drastically after 230bp.



# FILTER AND TRIM
filtFs <- file.path(path, "filtered", basename(cutFs))
filtRs <- file.path(path ,"filtered", basename(cutRs))


out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxEE=c(2,2), truncLen=c(200,200), truncQ=2, maxN=0, rm.phix=TRUE,
                     minLen = 200, compress=TRUE, verbose=TRUE, multithread=FALSE)  # on windows, set multithread = FALSE
head(out)


#Dada2's standard error model does not work properly with the binned quality data
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
  nbases = 1e8,
  errorEstimationFunction = loessErrfun_mod1,
  verbose = TRUE
)

# For Reverse
errR_1 <- learnErrors(
  filtRs,
  multithread = FALSE,
  nbases = 1e8,
  errorEstimationFunction = loessErrfun_mod1,
  verbose = TRUE
)

## The error check is the most time consuming step and I am not able to do it on my laptop.
# I changed nbases to 1e8 to speed stuf up, this works for the test, but maybe not a good
# idea for the full dataset.



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
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 240:265]
#Should delete at least the larger ones. If >300 bp, primers will overlap and are therefore not cut of?


#REMOVING CHMIMERAS. This is easier if the ASVs are there than before.
#If you have removed some fractions, use seqtab2.
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
taxtab <- assignTaxonomy(seqtab.nochim, "C:/Users/sjoerd.gremmen/Documents/R-scripts/D2287E8D3A112827DDED87E4E1AC3370CCA377088CA8EC24ED2C2C09EC0822DE.tar", multithread=TRUE)
# Takes <1h
write.csv(taxtab, "C:/Users/sjoerd.gremmen/Documents/ITS1/Output/ITS_Taxa_test.csv", row.names=TRUE)

#add species. Same step as before, download the file first.
# Only looks at exact matches.
########   taxa_species <- addSpecies(taxa, "~/data/volume_2/silva_species_assignment_v138.1.fa.gz", n = 1e3)
#This is a very RAM intensive step. uses more than this workplace can handle.
# I have tried playing with the n= part, going as low as 1e3
# Up here uses a bayesian method 
# So this doesn't work, we have to split our data

## Check if things are alright
unqs.mock <- seqtab.nochim["PositiveControl_TM030_F_filt.fastq",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

# Convert tax table and sample data to data frames
df.combined <- cbind(taxtab, t(seqtab.nochim))
# Print the output table
write.csv(df.combined, "C:/Users/sjoerd.gremmen/Documents/ITS1/Output/ASV_ITS_final.csv", row.names=TRUE)


