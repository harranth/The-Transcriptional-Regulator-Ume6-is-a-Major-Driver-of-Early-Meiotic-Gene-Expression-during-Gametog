##### THIS SCRIPT GENERATES DIFFERENTIAL EXPRESSION RESULTS FOR PAIRWISE #####
##### ANALYSIS OF aGFP TAGGED UME6(T99N) STRAINS WITH EITHER GFP TAGGED OR 
##### UNTAGGED IME1 ALLELES AT 2H #####

##### PATH TO WORKING DIRECTORY #####

setwd("/PATH")

##-------Designate File Location-------##
filename <-"PATH/"

##-------Designate File Output---------##
output <- "PATH/"


# Load Data - Select a file and read the data into a data-frame
my_data <- read.csv(filename, sep=",", header=TRUE)

# Rename data files
names(my_data)[2:97] <- c("WT1_0","WT1_2","WT1_4","WT1_6",
                          "WT2_0","WT2_2","WT2_4","WT2_6",
                          "WT3_0","WT3_2","WT3_4","WT3_6",
                          
                          "Ume6Tag1_0","Ume6Tag1_2","Ume6Tag1_4","Ume6Tag1_6",
                          "Ume6Tag2_0","Ume6Tag2_2","Ume6Tag2_4","Ume6Tag2_6",
                          "Ume6Tag3_0","Ume6Tag3_2","Ume6Tag3_4","Ume6Tag3_6",
                          
                          "Ume6T99NTag1_0","Ume6T99NTag1_2","Ume6T99NTag1_4","Ume6T99NTag1_6",
                          "Ume6T99NTag2_0","Ume6T99NTag2_2","Ume6T99NTag2_4","Ume6T99NTag2_6",
                          "Ume6T99NTag3_0","Ume6T99NTag3_2","Ume6T99NTag3_4","Ume6T99NTag3_6",
                          
                          "Ume6T99NTagaGFP1_0","Ume6T99NTagaGFP1_2","Ume6T99NTagaGFP1_4","Ume6T99NTagaGFP1_6",
                          "Ume6T99NTagaGFP2_0","Ume6T99NTagaGFP2_2","Ume6T99NTagaGFP2_4","Ume6T99NTagaGFP2_6",
                          "Ume6T99NTagaGFP3_0","Ume6T99NTagaGFP3_2","Ume6T99NTagaGFP3_4","Ume6T99NTagaGFP3_6",
                          
                          "sfGFPIME1_WT1_0","sfGFPIME1_WT1_2","sfGFPIME1_WT1_4","sfGFPIME1_WT1_6",
                          "sfGFPIME1_WT2_0","sfGFPIME1_WT2_2","sfGFPIME1_WT2_4","sfGFPIME1_WT2_6",
                          "sfGFPIME1_WT3_0","sfGFPIME1_WT3_2","sfGFPIME1_WT3_4","sfGFPIME1_WT3_6",
                         
                          "sfGFPIME1_Ume6Tag1_0","sfGFPIME1_Ume6Tag1_2","sfGFPIME1_Ume6Tag1_4","sfGFPIME1_Ume6Tag1_6",
                          "sfGFPIME1_Ume6Tag2_0","sfGFPIME1_Ume6Tag2_2","sfGFPIME1_Ume6Tag2_4","sfGFPIME1_Ume6Tag2_6",
                          "sfGFPIME1_Ume6Tag3_0","sfGFPIME1_Ume6Tag3_2","sfGFPIME1_Ume6Tag3_4","sfGFPIME1_Ume6Tag3_6",
                          
                          "sfGFPIME1_Ume6T99NTag1_0","sfGFPIME1_Ume6T99NTag1_2","sfGFPIME1_Ume6T99NTag1_4","sfGFPIME1_Ume6T99NTag1_6",
                          "sfGFPIME1_Ume6T99NTag2_0","sfGFPIME1_Ume6T99NTag2_2","sfGFPIME1_Ume6T99NTag2_4","sfGFPIME1_Ume6T99NTag2_6",
                          "sfGFPIME1_Ume6T99NTag3_0","sfGFPIME1_Ume6T99NTag3_2","sfGFPIME1_Ume6T99NTag3_4","sfGFPIME1_Ume6T99NTag3_6",
                          
                          "sfGFPIME1_Ume6T99NTagaGFP1_0","sfGFPIME1_Ume6T99NTagaGFP1_2","sfGFPIME1_Ume6T99NTagaGFP1_4","sfGFPIME1_Ume6T99NTagaGFP1_6",
                          "sfGFPIME1_Ume6T99NTagaGFP2_0","sfGFPIME1_Ume6T99NTagaGFP2_2","sfGFPIME1_Ume6T99NTagaGFP2_4","sfGFPIME1_Ume6T99NTagaGFP2_6",
                          "sfGFPIME1_Ume6T99NTagaGFP3_0","sfGFPIME1_Ume6T99NTagaGFP3_2","sfGFPIME1_Ume6T99NTagaGFP3_4","sfGFPIME1_Ume6T99NTagaGFP3_6"
                          )


######### Subset the 2h columns into a data frame ###########
data <- data.frame(
  my_data$Ume6T99NTagaGFP1_2,my_data$Ume6T99NTagaGFP2_2,my_data$Ume6T99NTagaGFP3_2,
  my_data$sfGFPIME1_Ume6T99NTagaGFP1_2,my_data$sfGFPIME1_Ume6T99NTagaGFP2_2,my_data$sfGFPIME1_Ume6T99NTagaGFP3_2,
  row.names = my_data$gene_id)   

##---------- NAME FOR FILE OUTPUT -------##
condition <- "Ume6_T99N_Rescue_IME1_TagvsUntag_UME6_T99N_3v5_antiGFP_2h"

conditions <- data.frame(row.names=names(data),
                         time=c(replicate(3, c("2h", "2h", "2h"))),
                         IME1_allele=c(replicate(3, "IME1"), 
                                       replicate(3, "sfGFP.IME1")),
                         
                         UME6_allele=c(replicate(6, "UME6.T99N.3v5.antiGFP")), 
                         replicate=c("rep1","rep2","rep3", 
                                     "rep1","rep2","rep3"))


#########DESEQ Setup#########
write.csv(conditions, paste(output, "/Conditions_Setup/", "Conditions.frame_", condition, ".csv", sep= ""))

######### DESEQ ANALYSIS #########
write.csv(conditions, paste(output,"Ume6.conditions.frame_", condition, ".csv", sep= ""))

##########PAIRWISE COMPARISON##########
##---------Generate Your DataMatrix----------###
# BiocManager::install("DESeq2")
# library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = conditions,
                              design = ~ IME1_allele)
##---------(Filtering) Removing Genes--------##
keep <- rowSums(counts(dds)) >=25
dds <- dds[keep,]

dds <- DESeqDataSet(dds, ~ IME1_allele)

##---------Test for Differential Expression---------##
ddsDE <- DESeq(dds)

###Export Normalized Read Counts
##-----------Normalization--------------##
normCounts <- counts(ddsDE, normalized = T)

###---------Export normliazed Read Counts---------###
write.csv(normCounts, paste(output, "normal_", condition, ".csv", sep= ""))

## Sort DEseq2 Results by padj
res <- results(ddsDE, alpha = 0.05)
summary(res)
res
resOrdered <- res[order(res$padj),]

# Save results as .csv
write.csv(resOrdered, paste(output, "deSeq_", condition, ".csv", sep= ""))

# Create variable associated with DESeq2 results
deSeqRes <- read.csv(paste(output, "deSeq_", condition, ".csv", sep= ""), row.names = 1)

# Define cut offs based on padj and log2FC
deSeqRes$diffexpressed <- "Not Significant"
deSeqRes$diffexpressed[deSeqRes$log2FoldChange > 1.5 & deSeqRes$padj < 0.05] <- "Up Regulated"
deSeqRes$diffexpressed[deSeqRes$log2FoldChange < -1.5 & deSeqRes$padj < 0.05] <- "Down Regulated"

# Add cut offs to DESeq2 output
write.csv(deSeqRes, paste(output,"deSeq_", condition, ".csv", sep= ""))


