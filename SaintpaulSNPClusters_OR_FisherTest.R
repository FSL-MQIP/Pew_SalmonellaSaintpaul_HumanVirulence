##-----------------------------------------------------------
## Project: PEW - Salmonella virulence
## Script purpose: Odds ratio calculation and Fisher's exact test for
## the S. Saintpaul SNP clusters.

## Start Date:  Aug 9, 2021
##-----------------------------------------------------------
## Notes: 
##
##-----------------------------------------------------------
## Load packages.
##-----------------------------------------------------------
# Install essential packages if haven't already.
required_packages <- c('tidyverse', 'readr')
for (p in required_packages) {
  if(!require(p,character.only = TRUE)) {
    install.packages(p, dep = TRUE)
  }
}
# Load the packages.
library(tidyverse);library(readr)
##-----------------------------------------------------------

## Load packages.
##-----------------------------------------------------------
#library(NCmisc);library(plyr);library(tidyverse);library(ggeffects);library(ggplot2);library(ggpubr);library(stringr);library(stringi)
##-----------------------------------------------------------
# Load the input file.
Saintpaul_SNPClustersStatsV3_RC031121 <- read_csv("Candidate_SaintpaulSNPCluster_Stats.csv")
Saintpaul_SNPClustersStatsV3_RC031121_Filtered <- Saintpaul_SNPClustersStatsV3_RC031121 %>% filter(Serovar == "Saintpaul")

# A loop for calculating OR and performing Fisher's exact tests.
for (i in 1:nrow(Saintpaul_SNPClustersStatsV3_RC031121_Filtered)) {
  
  InSNPCluster_Clinical = Saintpaul_SNPClustersStatsV3_RC031121_Filtered$Clinical_Number[i]
  InSNPCluster_Environmental = Saintpaul_SNPClustersStatsV3_RC031121_Filtered$Environmental_Number[i]
  OutSNPCluster_Clinical = sum(Saintpaul_SNPClustersStatsV3_RC031121_Filtered$Clinical_Number) - InSNPCluster_Clinical
  OutSNPCluster_Environmental = sum(Saintpaul_SNPClustersStatsV3_RC031121_Filtered$Environmental_Number) - InSNPCluster_Environmental
  
  Odds_Ratio = (InSNPCluster_Clinical/InSNPCluster_Environmental)/(OutSNPCluster_Clinical/OutSNPCluster_Environmental)
  Fisher_Table <- matrix(rep(NA,4), nrow = 2)
  colnames(Fisher_Table) <- c("Clinical", "Environmental")
  rownames(Fisher_Table) <- c("In_SNPCluster", "Out_SNPCluster")
  Fisher_Table[1,1] = InSNPCluster_Clinical
  Fisher_Table[1,2] = InSNPCluster_Environmental
  Fisher_Table[2,1] = OutSNPCluster_Clinical
  Fisher_Table[2,2] = OutSNPCluster_Environmental
  
  Results = ifelse(Odds_Ratio > 1, fisher.test(Fisher_Table, alternative = "greater"), fisher.test(Fisher_Table, alternative = "less"))
  
  Saintpaul_SNPClustersStatsV3_RC031121_Filtered$Odds_Ratio[i] = Odds_Ratio
  Saintpaul_SNPClustersStatsV3_RC031121_Filtered$FisherTest_P_Naive[i] = Results$p.value
  
}

Saintpaul_SNPClustersStatsV3_RC031121_Filtered$FisherTest_P_BH <- p.adjust(Saintpaul_SNPClustersStatsV3_RC031121_Filtered$FisherTest_P_Naive, method = "BH")
Saintpaul_SNPClustersStatsV3_RC031121_Filtered$FisherTest_P_BF <- p.adjust(Saintpaul_SNPClustersStatsV3_RC031121_Filtered$FisherTest_P_Naive, method = "bonferroni")



