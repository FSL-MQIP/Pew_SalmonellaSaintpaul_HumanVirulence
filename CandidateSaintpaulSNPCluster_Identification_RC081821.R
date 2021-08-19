##-----------------------------------------------------------
## Project: PEW - Salmonella virulence
## Script purpose: identification of candidate S. Saintpaul SNP clusters for downstream analysis.

## Start Date:  Aug 9, 2021
##-----------------------------------------------------------
## Notes: 
##
##-----------------------------------------------------------
## Load packages.
##-----------------------------------------------------------
# Install essential packages if haven't already.
required_packages <- c('plyr', 'tidyverse', 'readr')
for (p in required_packages) {
  if(!require(p,character.only = TRUE)) {
    install.packages(p, dep = TRUE)
  }
}
# Load the packages.
library(plyr);library(tidyverse);library(readxl)
##-----------------------------------------------------------
# Load the merged metadata.
MergedMetadata2_Accessed062320 <- read_csv("MergedMetadata2_Accessed062320.csv")
##-----------------------------------------------------------
## S. Saintpaul.
Saintpaul_Data <- MergedMetadata2_Accessed062320 %>% filter(grepl("paul", serovar))
dim(Saintpaul_Data) # [1] 1744   15
## Get SNP clusters that include isolates assigned to the serovar "S. Saintpaul".
Saintpaul_SNPCluster <- unique(Saintpaul_Data$PDS_acc)
## Get all isolates in these SNP clusters.
Saintpaul_SNPCluster_Isolates <- MergedMetadata2_Accessed062320 %>% filter(PDS_acc %in% Saintpaul_SNPCluster)

## Get rid of SNP clusters that were likely from other serovars and had isolates misclassified as S. Saintpaul.
Saintpaul_SNPCluster_Isolates_Summary <- ddply(Saintpaul_SNPCluster_Isolates, c("PDS_acc", "serovar"), summarise, N = sum(!is.na(PDS_acc)))
Saintpaul_SNPCluster_Isolates_Summary_Filtered <- Saintpaul_SNPCluster_Isolates_Summary %>% filter(serovar != "")
nrow(Saintpaul_SNPCluster_Isolates_Summary_Filtered) # 538



# Calculate the ratio of Saintpaul isolates in each cluster.
Saintpaul_Ratio_Summary <- data.frame("SNP_Cluster" = rep(NA, length(unique(Saintpaul_SNPCluster_Isolates_Summary_Filtered$PDS_acc))), 
                       "Saintpaul_Ratio" = rep(NA, length(unique(Saintpaul_SNPCluster_Isolates_Summary_Filtered$PDS_acc))),
                       "Saintpaul_Number" = rep(NA, length(unique(Saintpaul_SNPCluster_Isolates_Summary_Filtered$PDS_acc))),
                       "Total_Number" = rep(NA, length(unique(Saintpaul_SNPCluster_Isolates_Summary_Filtered$PDS_acc))))
Saintpaul_Ratio_Summary$SNP_Cluster <- unique(Saintpaul_SNPCluster_Isolates_Summary_Filtered$PDS_acc)
# Manaully check the unique terms used in the serovar column to identify the ones associated with S. Saintpaul
Saintpaul_Serovar_Name <- c("Saint-paul", "Saintpaul", "Saitpaul", "SaintPaul 4, 5 12 : e, h : 1")
for (cluster in Saintpaul_Ratio_Summary$SNP_Cluster) {
  target_content <- Saintpaul_SNPCluster_Isolates_Summary_Filtered %>% filter(PDS_acc == cluster)
  numerator <- sum(target_content$N[which(target_content$serovar %in% Saintpaul_Serovar_Name)])
  denominator <- sum(target_content$N)
  ratio <- numerator/denominator
  Saintpaul_Ratio_Summary$Saintpaul_Ratio[which(Saintpaul_Ratio_Summary$SNP_Cluster==cluster)] <- ratio
  Saintpaul_Ratio_Summary$Saintpaul_Number[which(Saintpaul_Ratio_Summary$SNP_Cluster==cluster)] <- numerator
  Saintpaul_Ratio_Summary$Total_Number[which(Saintpaul_Ratio_Summary$SNP_Cluster==cluster)] <- denominator
}

dim(Saintpaul_Ratio_Summary) # [1] 239   4

Saintpaul_Ratio_Summary_HighConfidence <- subset(Saintpaul_Ratio_Summary, Saintpaul_Ratio > 0.1)

# Candidate SNP clusters to retain.
Candidate_SaintpaulSNPClusters <- Saintpaul_Ratio_Summary_HighConfidence$SNP_Cluster
# These candidate SNP clusters were subject to the final inspection - 
# (1) A representative isolate was selected from each SNP cluster.
# (2) The serotype of the representative isolates was checked using SISTR.
# (3) The SNP clusters whose representative isolate was confirmed to be S. Saintpaul were considered
# S. Saintpaul SNP clusters and kept in the downstream analyses.
