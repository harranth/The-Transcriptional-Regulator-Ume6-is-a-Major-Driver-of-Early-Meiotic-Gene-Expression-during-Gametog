##### THIS SCRIPT GENERATES DIFFERENTIAL EXPRESSION RESULTS FOR PAIRWISE #####
##### ANALYSIS OF CONTROL AND UME6 DEPLETE STRAINS AT-30MIN #####

##### PATH TO WORKING DIRECTORY #####
setwd("/PATH")

##### PATH TO RAW COUNTS DATA #####
filename <-"/PATH"

##### PATH TO OUTPUT FOLDER #####
output <- "/PATH"

##### READ DATA INTO R #####
# Load Data - Select a file and read the data into a data-frame
my_data <- read.csv(filename, sep=",", header=TRUE)

names(my_data)[2:43] <- c("NoTIR1_min30","NoTIR1_0","NoTIR1_15","NoTIR1_30","NoTIR1_60","NoTIR1_120",
                          "TIR1_min30","TIR1_0","TIR1_15","TIR1_30","TIR1_60","TIR1_120",
                          "NoTIR2_min30","NoTIR2_0","NoTIR2_15","NoTIR2_30","NoTIR2_60","NoTIR2_120",
                          "TIR2_min30","TIR2_0","TIR2_15","TIR2_30","TIR2_60","TIR2_120",
                          "NoTIR3_min30","NoTIR3_0","NoTIR3_15","NoTIR3_30","NoTIR3_60","NoTIR3_120",
                          "TIR3_min30","TIR3_0","TIR3_15","TIR3_30","TIR3_60","TIR3_120",
                          "WTrep1","WTrep2","WTrep3",
                          "Delrep1","Delrep2","Delrep3")

data <- data.frame(my_data$TIR1_min30,my_data$TIR2_min30, my_data$TIR3_min30,
                   my_data$NoTIR1_min30,my_data$NoTIR2_min30,my_data$NoTIR3_min30,
                   row.names = my_data$gene_id)




##----------Now set up your DESeq2 comparison-------##
condition <- "CONTROL_VS_UME6_DEPLETE_MIN30"
conditions <- data.frame(row.names=names(data),
                         UME6_allele=c(replicate(3, "TIR"), 
                                       replicate(3, "NoTIR")),
                         time=c(replicate(2,c("min30min","min30min","min30min"))
                                ))

######### DESEQ ANALYSIS #########
write.csv(conditions, paste(output,"Ume6.conditions.frame_", condition, ".csv", sep= ""))

##---------Generate Your DataMatrix----------###
# BiocManager::install("DESeq2")
# library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = data,
                               colData = conditions,
                               design = ~ UME6_allele)
##---------(Filtering) Removing Genes--------##
keep <- rowSums(counts(dds)) >= 25
dds <- dds[keep,]

dds <- DESeqDataSet(dds, ~ UME6_allele)

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
deSeqRes$diffexpressed[deSeqRes$log2FoldChange > 0.6 & deSeqRes$padj < 0.05] <- "Up Regulated"
deSeqRes$diffexpressed[deSeqRes$log2FoldChange < -0.6 & deSeqRes$padj < 0.05] <- "Down Regulated"

# Add cut offs to DESeq2 output
write.csv(deSeqRes, paste(output,"deSeq_", condition, ".csv", sep= ""))