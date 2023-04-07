library(tidyr)
library(dplyr)

### Set Working Directory
setwd("PATH TO WORKING DIRECTORY")

### Set Directory Names
data <- "PATH TO DATA"
output <- paste("PATH TO OUTPUT", sep = "")

### Read in Files
#Average TPM Files
TPM_Avg <- read.csv(paste(data, "TPMAvg_Kallisto_Mitotic_Deplete.csv",sep = ""))#, row.names = 1)
TPM_Avg$gene_name <-  make.unique(TPM_Avg$gene_name, sep = "_")
row.names(TPM_Avg) <- TPM_Avg$gene_name

#DESeq2 Results (Ume6 deplete vs. control; 0 to 120min; padj > 0.05)
DESeq_TC_Targets <- read.csv(paste("deSeq_UME6_AID_TIR_vs_No_TIR_0min_Start_TIMECOURSE.csv", sep = ""))#, row.names = 1)
DESeq_TC_Targets_Sub <- subset(DESeq_TC_Targets, select = c(log2FoldChange, padj, gene_symbol))

#Isolate 177 DEGs from full DESeq2 time series list 
DESeq_TC_Targets_Filt <- filter(DESeq_TC_Targets_Sub, padj < 0.05)
row.names(DESeq_TC_Targets_Filt) <- DESeq_TC_Targets_Filt$gene_symbol

#Read in pairwise DESeq2 analysis (Ume6 deplete vs. control; -30min)
DE_min30 <- read.csv(paste("deSeq_UME6_AID_TIR_vs_No_TIR_min30.csv", sep = ""))#, row.names = 1)

#Read in pairwise DESeq2 analysis (Ume6 deplete vs. control; 0min)
DE_0 <- read.csv(paste("deSeq_UME6_AID_TIR_vs_No_TIR_0min.csv", sep = ""))#, row.names = 1)

#Read in pairwise DESeq2 analysis (Ume6 deplete vs. control; 15min)
DE_15 <- read.csv(paste("deSeq_UME6_AID_TIR_vs_No_TIR_15min.csv", sep = ""))#, row.names = 1)


### Filter Gene Lists (Log2FC and padj)
##################################################################################
# -30min list Filter
DE_min30_Sub <- filter(DE_min30, abs(log2FoldChange) > 1 & padj < 0.05)
DE_min30_Sub$gene_symbol <-  make.unique(DE_min30_Sub$gene_symbol, sep = "_")
row.names(DE_min30_Sub) <- DE_min30_Sub$gene_symbol
# PHO92, HSP82, SAE3 Not included in DESeq min30 results due to low expression

### Merge Files to Find Overlap
DE_min30_Overlap <- merge(DESeq_TC_Targets_Filt, DE_min30_Sub, by = 0)
(n <- nrow(DE_min30_Overlap))

##################################################################################
# 0 min list Filter
DE_0_Sub <- filter(DE_0, abs(log2FoldChange) > 1 & padj < 0.05)
DE_0_Sub$gene_symbol <-  make.unique(DE_0_Sub$gene_symbol, sep = "_")
row.names(DE_0_Sub) <- DE_0_Sub$gene_symbol
# SAE3 Not included in DESeq 0min results due to low expression

### Merge Files to Find Overlap
(DE_0_Overlap <- merge(DESeq_TC_Targets_Filt, DE_0_Sub, by = 0))
(n <- nrow(DE_0_Overlap))


##################################################################################
# 15 min Filter - Identify genes showing acute expression to Ume6 depletion at 15 min
DE_15_Sub <- filter(DE_15, abs(log2FoldChange) > 0.3 & padj < 0.05)
DE_15_Sub$gene_symbol <-  make.unique(DE_15_Sub$gene_symbol, sep = "_")
row.names(DE_15_Sub) <- DE_15_Sub$gene_symbol

### Merge Files to Find Overlap
DE_15_Overlap <- merge(DESeq_TC_Targets_Filt, DE_15_Sub, by = 0)
(n <- nrow(DE_15_Overlap))

#Write to folder (DESeq_Filtered_List)
write.csv(DE_15_Overlap, paste(output, "DESeq2_List_Filtered", ".csv", sep= ""))

################################## TPM CALCULATIONS ##################################################
#Find fold change between control and Ume6 deplete strains
TPM_Avg_Mutate <- mutate(TPM_Avg,
                         "Avg_TIR_OT" = (X.Tag...LexA...LexO.TIR..15min..Avg + 
                                            X.Tag...LexA...LexO.TIR..30min..Avg + 
                                            X.Tag...LexA...LexO.TIR..60min..Avg + 
                                            X.Tag...LexA...LexO.TIR..120min..Avg)/4,
                           
                         "Avg_NoTIR_OT" = (X.Tag...LexA..15min..Avg + 
                                              X.Tag...LexA..30min..Avg + 
                                              X.Tag...LexA..60min..Avg + 
                                              X.Tag...LexA..120min..Avg)/4,
                         
                         "Avg_15min_FC" = X.Tag...LexA...LexO.TIR..15min..Avg / X.Tag...LexA..15min..Avg,
                         "Avg_30min_FC" = X.Tag...LexA...LexO.TIR..30min..Avg / X.Tag...LexA..30min..Avg,
                         "Avg_60min_FC" = X.Tag...LexA...LexO.TIR..60min..Avg / X.Tag...LexA..60min..Avg,
                         "Avg_120min_FC" = X.Tag...LexA...LexO.TIR..120min..Avg / X.Tag...LexA..120min..Avg
                         )

#Find average expression of fold change from 15min to 120min
TPM_Avg_Mutate <- mutate(TPM_Avg_Mutate,
         "Avg_FC_15_to_120min" = (Avg_15min_FC + 
                                     Avg_30min_FC + 
                                     Avg_60min_FC + 
                                     Avg_120min_FC)/4)

#Find genes with high sustained expression
TPM_Avg_Mutate_Filt <- filter(TPM_Avg_Mutate,
                              Avg_FC_15_to_120min > 1.4
                              | Avg_FC_15_to_120min < 0.6)

# DESeq List (177) Overlap
TPM_Avg_Mutate_Filt_Overlap <- merge(DESeq_TC_Targets_Filt, TPM_Avg_Mutate_Filt, by = 0)
(n <- nrow(TPM_Avg_Mutate_Filt_Overlap))

#Write to folder (TPM_Filtered_List)
write.csv(TPM_Avg_Mutate_Filt_Overlap_Mitotic, paste(output, "TPM_List_Filtered", ".csv", sep= ""))


