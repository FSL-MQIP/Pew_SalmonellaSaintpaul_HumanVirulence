##-----------------------------------------------------------
## Project: PEW - Salmonella virulence
## Script purpose: Statistical Analysis for the invasion assay.

## Start Date:  Aug 9, 2021
##-----------------------------------------------------------
## Notes: 
## (1) Candidate factors: epi-type, lineage, and genotype of
## the genomic signatures listed in Tables 3 and 4.
##-----------------------------------------------------------
## Load packages.
##-----------------------------------------------------------
# Install essential packages if haven't already.
required_packages <- c('plyr', 'tidyverse', 'emmeans', 'leaps', 'readr')
for (p in required_packages) {
  if(!require(p,character.only = TRUE)) {
    install.packages(p, dep = TRUE)
  }
}
# Load the packages.
library(plyr);library(tidyverse);library(emmeans);library(leaps);library(readr)
##-----------------------------------------------------------
## Load and process experimental data.
InputFile1 <- read_csv("InvasionExpData_StatisticalAnalysis_RC072521.csv")
InputFile1=na.omit(InputFile1)
InputFile1[,c(1:2,4:21)] <- as.data.frame(lapply(InputFile1[,c(1:2,4:21)], FUN = as.factor))
# Remove genomic signatures that (i) have only one genotype across test strains and (ii) perfectly correlated with lineages.
InputFile1_Selected <- InputFile1 %>% select(-c(STM14_0610, ydiV, orf319, malX, yggR, sspH1))
##-----------------------------------------------------------
## Rank the genomic signatures by their correlation with the response variable.

# Create an empty dataframe.
Corr_Response <- expand.grid(Gene = colnames(InputFile1_Selected)[7:15])
Corr_Response$R_Value <- NA
# Calculate the correlation (R-value) of the genomic signatures.
for (gene in Corr_Response$Gene) {
  mod = lm(InputFile1_Selected$Log10_N_N0 ~ InputFile1_Selected[[gene]])
  mod_smry = summary(mod)
  r_squared = mod_smry$r.squared
  r_value = sqrt(r_squared)
  Corr_Response$R_Value[which(Corr_Response$Gene == gene)] = r_value
}
# Order the genomic signatures by their correlation from highest to lowest.
Corr_Response_Ordered <- Corr_Response %>% arrange(desc(R_Value))

# Manually inspect and re-load the table: combine the genomic signatures sharing identical genotypic patterns across the strains.
write.csv(Corr_Response_Ordered, file = "Corr_Response_Ordered_IA_RC080921.csv", row.names = FALSE)
Corr_Response_Ordered_Filtered <- read.csv("Corr_Response_Ordered_IA_RC080921.csv")

# Also combine the genomic signatures in InputFile1_Selected.
InputFile1_Selected_2 <- InputFile1_Selected %>% select(-c(sseI, bepF, stfC)) %>% rename(gtgE_sseI = gtgE, ygcH_stfC = ygcH, bglA_bepF = bglA)
##-----------------------------------------------------------
## Remove correlated genomic signatures.

# Assess pairwise correlation among the genomic signatures.
Corr_Mat <- matrix(nrow = nrow(Corr_Response_Ordered_Filtered), ncol = nrow(Corr_Response_Ordered_Filtered))
rownames(Corr_Mat) <- colnames(Corr_Mat) <- Corr_Response_Ordered_Filtered$Gene
for (i in 1:nrow(Corr_Mat)) {
  for (j in 1:ncol(Corr_Mat)) {
    if (i >= j) {next} else {
      gene1 = rownames(Corr_Mat)[i]
      gene2 = colnames(Corr_Mat)[j]
      results = fisher.test(InputFile1_Selected_2[[gene1]], InputFile1_Selected_2[[gene2]])
      p = results$p.value
      Corr_Mat[i,j] = p
    }
  }
}
Corr_Mat_DF <- as.data.frame(Corr_Mat)

# Manually transform Corr_Mat_DF for calculating BH-corrected P-values.
write.csv(Corr_Mat_DF, file = "Corr_Mat_DF_IA_RC080921.csv", row.names = TRUE)
# Re-load the transformed data set and calculate BH-corrected P-values.
Corr_Mat_DF_Transformed <- read.csv("Corr_Mat_DF_IA_RC080921.csv")
Corr_Mat_DF_Transformed$BH_P <- p.adjust(Corr_Mat_DF_Transformed$Naive_P, method = "BH")
# Remove correlated genomic signatures.
for (gene in unique(Corr_Mat_DF_Transformed$Gene)) {
  target_content = Corr_Mat_DF_Transformed %>% filter(Gene == gene)
  if (gene == unique(Corr_Mat_DF_Transformed$Gene)[1]) {next} else {
    if (sum(target_content$BH_P < 0.05) != 0) {
      print(gene)
    } else {next}
  }
} # Genes that need to be excluded: gtgA.
InputFile1_Selected_3 = subset(InputFile1_Selected_2, select = -c(gtgA))
levels(InputFile1_Selected_3$Rep) <- c("1", "2", "3", "4")
##-----------------------------------------------------------
## Select variables using the "best subset selection" method.

# Implement the best subset selection algorithm.
Best_Subset <- regsubsets(Log10_N_N0~.-Isolate-Odds_Ratio-Plating_Order,InputFile1_Selected_3, nvmax=11, force.in = c(1:3)) #all subset selection 
Best_Subset_Summary <- summary(Best_Subset)
# Plot the model performances based on different metrics.
par(mfrow=c(2,2))
plot(Best_Subset_Summary$rss,xlab="Number of Variables",ylab="RSS",type="l")
plot(Best_Subset_Summary$adjr2,xlab="Number of Variables",ylab="Adjusted R-Squared",type="l")
which.max(Best_Subset_Summary$adjr2)
points(4,Best_Subset_Summary$adjr2[4], col="red",cex=2,pch=20)
plot(Best_Subset_Summary$cp,xlab="Number of Variables",ylab="Cp",type='l')
which.min(Best_Subset_Summary$cp)
points(4,Best_Subset_Summary$cp[4],col="red",cex=2,pch=20)
which.min(Best_Subset_Summary$bic)
plot(Best_Subset_Summary$bic,xlab="Number of Variables",ylab="BIC",type='l')
points(4,Best_Subset_Summary$bic[4],col="red",cex=2,pch=20)
# BIC, Cp, and adjusted R-squared all suggested p = 4.
# Epi-type, lineage, and bglA_bepF should be included.
##-----------------------------------------------------------
## Build the INV model.
INV_Mod_1 <- lm(Log10_N_N0 ~ 1 + Rep + Epi_Type + Lineage + bglA_bepF, data = InputFile1_Selected_3)
# One-way ANOVA for the INV model.
anova(INV_Mod_1)
##-----------------------------------------------------------
## Model diagnostics for the INV model.
par(mfrow = c(2, 1), mar=c(4,4,1,1))
# Check linearity and homoscedasticity
plot(fitted(IA_Mod_1),residuals(IA_Mod_1), xlab = "Fitted Values", ylab = "Residuals")
# Check the residual normality
qqnorm(residuals(IA_Mod_1), ylab="Residuals", plot.it=TRUE, datax=FALSE, main = NULL)
qqline(resid(IA_Mod_1))
par(mfrow = c(1, 1))
##-----------------------------------------------------------
## Post-hoc analyses.
# Epi-type.
INV_Mod_1_Epi_Type.emm <- emmeans(INV_Mod_1, ~Epi_Type)
INV_Mod_1_Epi_Type.ctr <- contrast(INV_Mod_1_Epi_Type.emm, interaction = "pairwise", simple = "Epi_Type", combine = TRUE, adjust = "none")
# Lineage.
INV_Mod_1_Lineage.emm <- emmeans(INV_Mod_1, ~Lineage)
INV_Mod_1_Lineage.ctr <- contrast(INV_Mod_1_Lineage.emm, interaction = "pairwise", simple = "Lineage", combine = TRUE, adjust = "tukey")
# bglA_bepF
INV_Mod_1_bglA_bepF.emm <- emmeans(INV_Mod_1, ~bglA_bepF)
INV_Mod_1_bglA_bepF.ctr <- contrast(INV_Mod_1_bglA_bepF.emm, interaction = "pairwise", simple = "bglA_bepF", combine = TRUE, adjust = "none")
##-----------------------------------------------------------
