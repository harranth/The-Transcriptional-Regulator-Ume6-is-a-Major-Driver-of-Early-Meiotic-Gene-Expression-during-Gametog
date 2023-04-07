##### THIS SCRIPT GENERATES DIFFERENTIAL EXPRESSION RESULTS FOR PAIRWISE #####
##### ANALYSIS OF aGFP TAGGED UME6(T99N) STRAINS WITH EITHER GFP-IME1 OR 
##### GFP-B112 ALLELES AT 6H #####

##### PATH TO WORKING DIRECTORY #####
setwd("/PATH")

##-------Designate File Location-------##
filename <-"PATH/"

##-------Designate File Output---------##
output <- "PATH/"


my_data <- read.csv(filename, sep=",", header=TRUE)
names(my_data)[2:73] <- c("sfGFP_B112_0HR_Rep1","sfGFP_B112_2HR_Rep1","sfGFP_B112_4HR_Rep1","sfGFP_B112_6HR_Rep1",
                          "sfGFP_B112_0HR_Rep2","sfGFP_B112_2HR_Rep2","sfGFP_B112_4HR_Rep2","sfGFP_B112_6HR_Rep2",
                          "sfGFP_B112_0HR_Rep3","sfGFP_B112_2HR_Rep3","sfGFP_B112_4HR_Rep3","sfGFP_B112_6HR_Rep3",
                          "sfGFP_GAL4_75to881_0HR_Rep1","sfGFP_GAL4_75to881_2HR_Rep1","sfGFP_GAL4_75to881_4HR_Rep1","sfGFP_GAL4_75to881_6HR_Rep1",
                          "sfGFP_GAL4_75to881_0HR_Rep2","sfGFP_GAL4_75to881_2HR_Rep2","sfGFP_GAL4_75to881_4HR_Rep2","sfGFP_GAL4_75to881_6HR_Rep2",
                          "sfGFP_GAL4_75to881_0HR_Rep3","sfGFP_GAL4_75to881_2HR_Rep3","sfGFP_GAL4_75to881_4HR_Rep3","sfGFP_GAL4_75to881_6HR_Rep3",
                          "sfGFP_IME1_0HR_Rep1","sfGFP_IME1_2HR_Rep1","sfGFP_IME1_4HR_Rep1","sfGFP_IME1_6HR_Rep1",
                          "sfGFP_IME1_0HR_Rep2","sfGFP_IME1_2HR_Rep2","sfGFP_IME1_4HR_Rep2","sfGFP_IME1_6HR_Rep2",
                          "sfGFP_IME1_0HR_Rep3","sfGFP_IME1_2HR_Rep3","sfGFP_IME1_4HR_Rep3","sfGFP_IME1_6HR_Rep3",
                          "GAL4_75to881_0HR_Rep1","GAL4_75to881_2HR_Rep1","GAL4_75to881_4HR_Rep1","GAL4_75to881_6HR_Rep1",
                          "GAL4_75to881_0HR_Rep2","GAL4_75to881_2HR_Rep2","GAL4_75to881_4HR_Rep2","GAL4_75to881_6HR_Rep2",
                          "GAL4_75to881_0HR_Rep3","GAL4_75to881_2HR_Rep3","GAL4_75to881_4HR_Rep3","GAL4_75to881_6HR_Rep3",
                          "B112_0HR_Rep1","B112_2HR_Rep1","B112_4HR_Rep1","B112_6HR_Rep1",
                          "B112_0HR_Rep2","B112_2HR_Rep2","B112_4HR_Rep2","B112_6HR_Rep2",
                          "B112_0HR_Rep3","B112_2HR_Rep3","B112_4HR_Rep3","B112_6HR_Rep3",
                          "IME1_0HR_Rep1","IME1_2HR_Rep1","IME1_4HR_Rep1","IME1_6HR_Rep1",
                          "IME1_0HR_Rep2","IME1_2HR_Rep2","IME1_4HR_Rep2","IME1_6HR_Rep2",
                          "IME1_0HR_Rep3","IME1_2HR_Rep3","IME1_4HR_Rep3","IME1_6HR_Rep3"
)

######### Subset the 2h columns into a data frame ###########
data <- data.frame(
  my_data$sfGFP_IME1_6HR_Rep1,my_data$sfGFP_IME1_6HR_Rep2,my_data$sfGFP_IME1_6HR_Rep3,
  my_data$sfGFP_B112_6HR_Rep1,my_data$sfGFP_B112_6HR_Rep2,my_data$sfGFP_B112_6HR_Rep3,
  row.names = my_data$gene_id)   

##---------- NAME FOR FILE OUTPUT -------##
condition <- "Ume6_T99N_TransAD_Rescue_Tagged_IME1and_B112_at_6h"

conditions <- data.frame(row.names=names(data),
                         time=c(replicate(2, c("6h", "6h", "6h"))
                                ),
                         IME1_allele=c(replicate(3, "sfGFP.IME1"),
                                       replicate(3, "sfGFP.B112")
                                       ),
                         UME6_allele=c(replicate(6, "UME6.T99N.3v5.antiGFP")
                                       ),
                         replicate=c("rep1","rep2","rep3"))


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


