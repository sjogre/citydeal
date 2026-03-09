library(tidyverse)

## Setwd to where te files are
setwd("G:/My Drive/CityDeal/Data_analysis_after_review/Input_data/FA data/")
## Importing all files from fragment analyser, Stolen from Stackoverflow thanks to https://stackoverflow.com/questions/11433432/how-to-import-multiple-csv-files-at-once
temp = list.files(pattern="\\.csv$")
myfiles = lapply(temp, read.delim, sep = ",")
myfiles

FA_all <- bind_rows(myfiles)
FA_all <- separate(FA_all,
              Sample.ID, 
              sep = "_",
              into = c("enumber","plate","wells"))
subset(FA_all, enumber == "e1100038290")

setwd("G:/My Drive/CityDeal/Data_analysis_after_review")
## Now I want to link the enumbers to the samples, to make a nice metadata table
ASV_table_ITS <- read.csv("./Input_data/ITS_spring23_ASV.csv", row.names = 1, header = T)
Sample_names_ITS <- as.data.frame(colnames(ASV_table_ITS))
Sample_names_ITS$names <- Sample_names_ITS$`colnames(ASV_table_ITS)`
Sample_names_ITS$`colnames(ASV_table_ITS)`<- NULL
rm(ASV_table_ITS)
Sample_names_ITS <- separate(Sample_names_ITS,
         names,
         sep = "_",
         into= c("enumber", "NCB-number", "SampleID", "plate_ITS"))
## We need to merge 16S to enumber
ASV_table_16S <-read.csv("./Input_data/Full_ASV_table_16s.csv")
Sample_names_16S <- as.data.frame(colnames(ASV_table_16S))
colnames(Sample_names_16S) <- "sampleID_plate"
Sample_names_16S <- separate(Sample_names_16S,
                             sampleID_plate,
                             sep = "_",
                             into= c("SampleID", "plate_16S"))
Sample_names_16S$SampleID[Sample_names_16S$SampleID == "X"] <- "negcontrol"
Sample_names_16S <- Sample_names_16S[2:length(rownames(Sample_names_16S)),]
rm(ASV_table_16S)

## now merge
Combined_ITS_16S <- merge(Sample_names_ITS, Sample_names_16S, by = "SampleID")
##
Combined_ITS_16S$plate_combo <- paste(Combined_ITS_16S$plate_ITS, Combined_ITS_16S$plate_16S, sep = "_")
Combined_ITS_16S <- separate_rows(Combined_ITS_16S, plate_combo, sep = "_")
Combined_ITS_16S$marker <- rep(c("ITS", "16S"), times = 722)
Combined_ITS_16S$plate_ITS <- NULL
Combined_ITS_16S$plate_16S <- NULL

## Combining plate data with FA data
Combined_ITS_16S$enumber_plate <- paste(Combined_ITS_16S$enumber, Combined_ITS_16S$plate_combo, sep = "_")
FA_all$enumber_plate <- paste(FA_all$enumber, FA_all$plate, sep = "_")

Metadata_FA <- merge(FA_all, Combined_ITS_16S, by = "enumber_plate")

Metadata_FA <- Metadata_FA[!duplicated(Metadata_FA), ]
Metadata_FA$sample_control <- ifelse(Metadata_FA$SampleID == "negcontrol", 
                                     "neg_control", 
                                     "real_sample")
Metadata_FA_16S <- subset(Metadata_FA, marker == "16S")
Metadata_FA_ITS <- subset(Metadata_FA, marker == "ITS")
write.csv(Metadata_FA_16S, "./Input_data/FA_metadata_16S.csv")
write.csv(Metadata_FA_ITS, "./Input_data/FA_metadata_ITS.csv")
write.csv(Metadata_FA, "./Input_data/FA_metadata_Full.csv")
