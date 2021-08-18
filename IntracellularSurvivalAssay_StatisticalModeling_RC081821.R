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
InputFile1 <- read_csv("IntracellularSurvivalExpData_StatisticalAnalysis_RC073121.csv")
InputFile1_Filtered_Selected <- InputFile1 %>% filter(Serotype == "Saintpaul", Rep %in% c(1,2,3,4)) %>% select(-c(2, 4, 8:9,11:18, 20:26, 36:37, 31))
InputFile1_Filtered_Selected[,-c(3,7)] <- as.data.frame(lapply(InputFile1_Filtered_Selected[,-c(3,7)], FUN = as.factor))

InputFile1_Filtered_Selected_0hpi <- InputFile1_Filtered_Selected %>% filter(hpi == "0")
InputFile1_Filtered_Selected_2hpi <- InputFile1_Filtered_Selected %>% filter(hpi == "2")
InputFile1_Filtered_Selected_6hpi <- InputFile1_Filtered_Selected %>% filter(hpi == "6")
InputFile1_Filtered_Selected_24hpi <- InputFile1_Filtered_Selected %>% filter(hpi == "24")

# Create separate data sets for 0-2 hpi, 2-6 hpi, and 6-24 hpi.
ICS_Diff_0_2 <- ICS_Diff_2_6 <- ICS_Diff_6_24 <- InputFile1_Filtered_Selected_0hpi[,-c(6:7)]
ICS_Diff_0_2$Diff_0_2 <- InputFile1_Filtered_Selected_2hpi$Log10_CFU_ml - InputFile1_Filtered_Selected_0hpi$Log10_CFU_ml
ICS_Diff_2_6$Diff_2_6 <- InputFile1_Filtered_Selected_6hpi$Log10_CFU_ml - InputFile1_Filtered_Selected_2hpi$Log10_CFU_ml
ICS_Diff_6_24$Diff_6_24 <- InputFile1_Filtered_Selected_24hpi$Log10_CFU_ml - InputFile1_Filtered_Selected_6hpi$Log10_CFU_ml

##-----------------------------------------------------------
# Statistical analysis for 0-2 hpi.
##-----------------------------------------------------------
## Rank the genomic signatures by their correlation with the response variable.
Corr_Response_1 <- expand.grid(Gene = colnames(ICS_Diff_0_2)[6:14])
Corr_Response_1$R_Value <- NA

for (gene in Corr_Response_1$Gene) {
  mod = lm(ICS_Diff_0_2$Diff_0_2 ~ ICS_Diff_0_2[[gene]])
  mod_smry = summary(mod)
  r_squared = mod_smry$r.squared
  r_value = sqrt(r_squared)
  Corr_Response_1$R_Value[which(Corr_Response_1$Gene == gene)] = r_value
}
# Order the genomic signatures by their correlation from highest to lowest.
Corr_Response_1_Ordered <- Corr_Response_1 %>% arrange(desc(R_Value))
# Manually inspect and re-load the table: combine the genomic signatures sharing identical genotypic patterns across the strains.
write.csv(Corr_Response_1_Ordered, file = "Corr_Response_Ordered_ICS1_RC080921.csv", row.names = FALSE)
Corr_Response_1_Ordered_Filtered <- read.csv("Corr_Response_Ordered_ICS1_RC080921.csv")
# Also combine the genomic signatures in ICS_Diff_0_2.
ICS_Diff_0_2_Selected <- ICS_Diff_0_2 %>% select(-c(sseI, bepF, stfC)) %>% rename(gtgE_sseI = gtgE, ygcH_stfC = ygcH, bglA_bepF = bglA)
##-----------------------------------------------------------
## Remove correlated genomic signatures.

# Assess pairwise correlation among the genomic signatures.
Corr_Mat_1 <- matrix(nrow = nrow(Corr_Response_1_Ordered_Filtered), ncol = nrow(Corr_Response_1_Ordered_Filtered))
rownames(Corr_Mat_1) <- colnames(Corr_Mat_1) <- Corr_Response_1_Ordered_Filtered$Gene
for (i in 1:nrow(Corr_Mat_1)) {
  for (j in 1:ncol(Corr_Mat_1)) {
    if (i >= j) {next} else {
      gene1 = rownames(Corr_Mat_1)[i]
      gene2 = colnames(Corr_Mat_1)[j]
      results = fisher.test(ICS_Diff_0_2_Selected[[gene1]], ICS_Diff_0_2_Selected[[gene2]])
      p = results$p.value
      Corr_Mat_1[i,j] = p
    }
  }
}
Corr_Mat_1_DF <- as.data.frame(Corr_Mat_1)

# Manually transform Corr_Mat_1_DF for calculating BH-corrected P-values.
write.csv(Corr_Mat_1_DF, file = "Corr_Mat_DF_ICS1_RC080921.csv", row.names = TRUE)
# Re-load the transformed data set and calculate BH-corrected P-values.
Corr_Mat_DF_1_Transformed <- read.csv("Corr_Mat_DF_ICS1_RC080921.csv")
Corr_Mat_DF_1_Transformed$BH_P <- p.adjust(Corr_Mat_DF_1_Transformed$Naive_P, method = "BH")
# Remove correlated genomic signatures.
for (gene in unique(Corr_Mat_DF_1_Transformed$Gene)) {
  target_content = Corr_Mat_DF_1_Transformed %>% filter(Gene == gene)
  if (gene == unique(Corr_Mat_DF_1_Transformed$Gene)[1]) {next} else {
    if (sum(target_content$BH_P < 0.05) != 0) {
      print(gene)
    } else {next}
  }
} # Genes that need to be excluded: gtgA.
ICS_Diff_0_2_Selected_2 = subset(ICS_Diff_0_2_Selected, select = -c(gtgA))
##-----------------------------------------------------------
## Select variables using the "best subset selection" method.

# Implement the best subset selection algorithm.
Best_Subset_1 <- regsubsets(Diff_0_2~.-Strain-Odds_Ratio,ICS_Diff_0_2_Selected_2, nvmax=10, force.in = c(4:5)) #all subset selection 
Best_Subset_1_Summary <- summary(Best_Subset_1)
# Plot the model performances based on different metrics.
par(mfrow=c(2,2))
plot(Best_Subset_1_Summary$rss,xlab="Number of Variables",ylab="RSS",type="l")
plot(Best_Subset_1_Summary$adjr2,xlab="Number of Variables",ylab="Adjusted R-Squared",type="l")
which.max(Best_Subset_1_Summary$adjr2)
points(4,Best_Subset_1_Summary$adjr2[4], col="red",cex=2,pch=20)
plot(Best_Subset_1_Summary$cp,xlab="Number of Variables",ylab="Cp",type='l')
which.min(Best_Subset_1_Summary$cp)
points(3,Best_Subset_Summary$cp[3],col="red",cex=2,pch=20)
which.min(Best_Subset_1_Summary$bic)
plot(Best_Subset_1_Summary$bic,xlab="Number of Variables",ylab="BIC",type='l')
points(3,Best_Subset_Summary$bic[3],col="red",cex=2,pch=20)
# BIC and Cp suggested p = 3
# Adjusted R-squared suggested p = 4.

# The two models (with p=3 and p=4, respectively) are not nested, so use cv to determine which one is better.
ICS_Diff_0_2_Selected_2$Lineage_IC <- ifelse(ICS_Diff_0_2_Selected_2$Lineage == "IC", 1, 0)
ICS_Diff_0_2_Selected_2$Lineage_IC <- as.factor(ICS_Diff_0_2_Selected_2$Lineage_IC)
# Defining training control as cross-validation (K=5).
set.seed(08092021)
Train_Control <- trainControl(method = "cv", number = 5)
ICS1_Mod_p3_CV <- train(Diff_0_2 ~ 1 + Rep + Lineage_IC + bglA_bepF + ygcH_stfC, data = ICS_Diff_0_2_Selected_2, method = "lm", trControl = Train_Control)
print(ICS1_Mod_p3_CV)
#RMSE       Rsquared   MAE      
#0.2168302  0.5438822  0.1843196
ICS1_Mod_p4_CV <- train(Diff_0_2 ~ 1 + Rep + Epi_Type + sodCI + ygcH_stfC + bglA_bepF, data = ICS_Diff_0_2_Selected_2, method = "lm", trControl = Train_Control)
print(ICS1_Mod_p4_CV)
#RMSE       Rsquared   MAE      
#0.2288953  0.6576763  0.1878462
# The performances were similar, go with p=3 to avoid over-fitting.
##-----------------------------------------------------------
# Build the ICS-1 model.
ICS1_Mod_1 <- lm(Diff_0_2 ~ 1 + Rep + Lineage_IC + bglA_bepF + ygcH_stfC, data = ICS_Diff_0_2_Selected_2)
# One-way ANOVA for the ICS-1 model.
anova(ICS1_Mod_1)
##-----------------------------------------------------------
# Model diagnostics.
par(mfrow = c(2, 1), mar=c(4,4,1,1))
# Check linearity and homoscedasticity
plot(fitted(ICS1_Mod_1),residuals(ICS1_Mod_1), xlab = "Fitted Values", ylab = "Residuals")
# Check the residual normality
qqnorm(residuals(ICS1_Mod_1), ylab="Residuals", plot.it=TRUE, datax=FALSE, main = NULL)
qqline(resid(ICS1_Mod_1))
par(mfrow = c(1, 1))
##-----------------------------------------------------------
## Post-hoc analyses.
# Lineage.
ICS1_Mod_1_Lineage_IC.emm <- emmeans(ICS1_Mod_1, ~Lineage_IC)
ICS1_Mod_1_Lineage_IC.ctr <- contrast(ICS1_Mod_1_Lineage_IC.emm, interaction = "pairwise", simple = "Lineage_IC", combine = TRUE, adjust = "tukey")
##-----------------------------------------------------------
# Statistical analysis for 2-6 hpi.
##-----------------------------------------------------------
## Rank the genomic signatures by their correlation with the response variable.

# Create an empty dataframe.
Corr_Response_2 <- expand.grid(Gene = colnames(ICS_Diff_2_6)[6:14])
Corr_Response_2$R_Value <- NA
# Calculate the correlation (R-value) of the genomic signatures.
for (gene in Corr_Response_2$Gene) {
  mod = lm(ICS_Diff_2_6$Diff_2_6 ~ ICS_Diff_2_6[[gene]])
  mod_smry = summary(mod)
  r_squared = mod_smry$r.squared
  r_value = sqrt(r_squared)
  Corr_Response_2$R_Value[which(Corr_Response_2$Gene == gene)] = r_value
}
# Order the genomic signatures by their correlation from highest to lowest.
Corr_Response_2_Ordered <- Corr_Response_2 %>% arrange(desc(R_Value))
# Manually inspect and re-load the table: combine the genomic signatures sharing identical genotypic patterns across the strains.
write.csv(Corr_Response_2_Ordered, file = "Corr_Response_Ordered_ICS2_RC080921.csv", row.names = FALSE)
Corr_Response_2_Ordered_Filtered <- read.csv("Corr_Response_Ordered_ICS2_RC080921.csv")
# Also combine the genomic signatures in ICS_Diff_2_6.
ICS_Diff_2_6_Selected <- ICS_Diff_2_6 %>% select(-c(sseI, bepF, stfC)) %>% rename(gtgE_sseI = gtgE, ygcH_stfC = ygcH, bglA_bepF = bglA)
##-----------------------------------------------------------
## Remove correlated genomic signatures.

# Assess pairwise correlation among the genomic signatures.
Corr_Mat_2 <- matrix(nrow = nrow(Corr_Response_2_Ordered_Filtered), ncol = nrow(Corr_Response_2_Ordered_Filtered))
rownames(Corr_Mat_2) <- colnames(Corr_Mat_2) <- Corr_Response_2_Ordered_Filtered$Gene
for (i in 1:nrow(Corr_Mat_2)) {
  for (j in 1:ncol(Corr_Mat_2)) {
    if (i >= j) {next} else {
      gene1 = rownames(Corr_Mat_2)[i]
      gene2 = colnames(Corr_Mat_2)[j]
      results = fisher.test(ICS_Diff_2_6_Selected[[gene1]], ICS_Diff_2_6_Selected[[gene2]])
      p = results$p.value
      Corr_Mat_2[i,j] = p
    }
  }
}
Corr_Mat_2_DF <- as.data.frame(Corr_Mat_2)

# Manually transform Corr_Mat_2_DF for calculating BH-corrected P-values.
write.csv(Corr_Mat_2_DF, file = "Corr_Mat_DF_ICS2_RC080921.csv", row.names = TRUE)
# Re-load the transformed data set and calculate BH-corrected P-values.
Corr_Mat_2_DF_Transformed <- read.csv("Corr_Mat_DF_ICS2_RC080921.csv")
Corr_Mat_2_DF_Transformed$BH_P <- p.adjust(Corr_Mat_2_DF_Transformed$Naive_P, method = "BH")
# Remove correlated genomic signatures.
for (gene in unique(Corr_Mat_2_DF_Transformed$Gene)) {
  target_content = Corr_Mat_2_DF_Transformed %>% filter(Gene == gene)
  if (gene == unique(Corr_Mat_2_DF_Transformed$Gene)[1]) {next} else {
    if (sum(target_content$BH_P < 0.05) != 0) {
      print(gene)
    } else {next}
  }
} # Genes that need to be excluded: gtgA.
ICS_Diff_2_6_Selected_2 = subset(ICS_Diff_2_6_Selected, select = -c(gtgA))
##-----------------------------------------------------------
## Select variables using the "best subset selection" method.

# Implement the best subset selection algorithm.
Best_Subset_2 <- regsubsets(Diff_2_6~.-Strain-Odds_Ratio,ICS_Diff_2_6_Selected_2, nvmax=10, force.in = c(4:5)) #all subset selection 
Best_Subset_2_Summary <- summary(Best_Subset_2)
# Plot the model performances based on different metrics.
par(mfrow=c(2,2))
plot(Best_Subset_2_Summary$rss,xlab="Number of Variables",ylab="RSS",type="l")
plot(Best_Subset_2_Summary$adjr2,xlab="Number of Variables",ylab="Adjusted R-Squared",type="l")
which.max(Best_Subset_2_Summary$adjr2)
points(3,Best_Subset_2_Summary$adjr2[3], col="red",cex=2,pch=20)
plot(Best_Subset_2_Summary$cp,xlab="Number of Variables",ylab="Cp",type='l')
which.min(Best_Subset_2_Summary$cp)
points(2,Best_Subset_2_Summary$cp[2],col="red",cex=2,pch=20)
which.min(Best_Subset_2_Summary$bic)
plot(Best_Subset_2_Summary$bic,xlab="Number of Variables",ylab="BIC",type='l')
points(2,Best_Subset_2_Summary$bic[2],col="red",cex=2,pch=20)
# BIC and Cp suggested p = 2
# Adjusted R-squared suggested p = 3.

# The two models are nested, so use LRT to determine which one is better.
ICS2_Mod_p2 <- lm(Diff_2_6 ~ 1 + Rep + Lineage, data = ICS_Diff_2_6_Selected_2)
ICS2_Mod_p3 <- lm(Diff_2_6 ~ 1 + Rep + Lineage + ygcH_stfC, data = ICS_Diff_2_6_Selected_2)
anova(ICS2_Mod_p2, ICS2_Mod_p3)
# P = 0.2566
# We should go with p = 2.
##-----------------------------------------------------------
# Build the ICS-2 model.
ICS2_Mod_1 <- lm(Diff_2_6 ~ 1 + Rep + Lineage, data = ICS_Diff_2_6_Selected_2)
# One-way ANOVA for the ICS-2 model.
anova(ICS2_Mod_1)
##-----------------------------------------------------------
## Model diagnostics for the INV model.
par(mfrow = c(2, 1), mar=c(4,4,1,1))
# Check linearity and homoscedasticity
plot(fitted(ICS2_Mod_1),residuals(ICS2_Mod_1), xlab = "Fitted Values", ylab = "Residuals")
# Check the residual normality
qqnorm(residuals(ICS2_Mod_1), ylab="Residuals", plot.it=TRUE, datax=FALSE, main = NULL)
qqline(resid(ICS2_Mod_1))
par(mfrow = c(1, 1))
##-----------------------------------------------------------
# Post-hoc analyses.
# Lineage.
ICS2_Mod_1_Lineage.emm <- emmeans(ICS2_Mod_1, ~Lineage)
ICS2_Mod_1_Lineage.ctr <- contrast(ICS2_Mod_1_Lineage.emm, interaction = "pairwise", simple = "Lineage", combine = TRUE, adjust = "tukey")
##-----------------------------------------------------------
# Statistical analysis for 6-24 hpi.
##-----------------------------------------------------------
## Rank the genomic signatures by their correlation with the response variable.

# Create an empty dataframe.
Corr_Response_3 <- expand.grid(Gene = colnames(ICS_Diff_6_24)[6:14])
Corr_Response_3$R_Value <- NA
# Calculate the correlation (R-value) of the genomic signatures.
for (gene in Corr_Response_3$Gene) {
  mod = lm(ICS_Diff_6_24$Diff_6_24 ~ ICS_Diff_6_24[[gene]])
  mod_smry = summary(mod)
  r_squared = mod_smry$r.squared
  r_value = sqrt(r_squared)
  Corr_Response_3$R_Value[which(Corr_Response_3$Gene == gene)] = r_value
}
# Order the genomic signatures by their correlation from highest to lowest.
Corr_Response_3_Ordered <- Corr_Response_3 %>% arrange(desc(R_Value))
# Manually inspect and re-load the table: combine the genomic signatures sharing identical genotypic patterns across the strains.
write.csv(Corr_Response_3_Ordered, file = "Corr_Response_Ordered_ICS3_RC080921.csv", row.names = FALSE)
Corr_Response_3_Ordered_Filtered <- read.csv("Corr_Response_Ordered_ICS3_RC080921.csv")
# Also combine the genomic signatures in ICS_Diff_2_6.
ICS_Diff_6_24_Selected <- ICS_Diff_6_24 %>% select(-c(sseI, bepF, stfC)) %>% rename(gtgE_sseI = gtgE, ygcH_stfC = ygcH, bglA_bepF = bglA)
##-----------------------------------------------------------
## Remove correlated genomic signatures.

# Assess pairwise correlation among the genomic signatures.
Corr_Mat_3 <- matrix(nrow = nrow(Corr_Response_3_Ordered_Filtered), ncol = nrow(Corr_Response_3_Ordered_Filtered))
rownames(Corr_Mat_3) <- colnames(Corr_Mat_3) <- Corr_Response_3_Ordered_Filtered$Gene
for (i in 1:nrow(Corr_Mat_3)) {
  for (j in 1:ncol(Corr_Mat_3)) {
    if (i >= j) {next} else {
      gene1 = rownames(Corr_Mat_3)[i]
      gene2 = colnames(Corr_Mat_3)[j]
      results = fisher.test(ICS_Diff_6_24_Selected[[gene1]], ICS_Diff_6_24_Selected[[gene2]])
      p = results$p.value
      Corr_Mat_3[i,j] = p
    }
  }
}
Corr_Mat_3_DF <- as.data.frame(Corr_Mat_3)

# Manually transform Corr_Mat_3_DF for calculating BH-corrected P-values.
write.csv(Corr_Mat_3_DF, file = "Corr_Mat_DF_ICS3_RC080921.csv", row.names = TRUE)
# Re-load the transformed data set and calculate BH-corrected P-values.
Corr_Mat_3_DF_Transformed <- read.csv("Corr_Mat_DF_ICS3_RC080921.csv")
Corr_Mat_3_DF_Transformed$BH_P <- p.adjust(Corr_Mat_3_DF_Transformed$Naive_P, method = "BH")
# Remove correlated genomic signatures.
for (gene in unique(Corr_Mat_3_DF_Transformed$Gene)) {
  target_content = Corr_Mat_3_DF_Transformed %>% filter(Gene == gene)
  if (gene == unique(Corr_Mat_3_DF_Transformed$Gene)[1]) {next} else {
    if (sum(target_content$BH_P < 0.05) != 0) {
      print(gene)
    } else {next}
  }
} # Genes that need to be excluded: gtgA.
ICS_Diff_6_24_Selected_2 = subset(ICS_Diff_6_24_Selected, select = -c(gtgA))
##-----------------------------------------------------------
## Select variables using the "best subset selection" method.

# Implement the best subset selection algorithm.
Best_Subset_3 <- regsubsets(Diff_6_24~.-Strain-Odds_Ratio,ICS_Diff_6_24_Selected_2, nvmax=10, force.in = c(4:5)) #all subset selection 
Best_Subset_3_Summary <- summary(Best_Subset_3)
# Plot the model performances based on different metrics.
par(mfrow=c(2,2))
plot(Best_Subset_3_Summary$rss,xlab="Number of Variables",ylab="RSS",type="l")
plot(Best_Subset_3_Summary$adjr2,xlab="Number of Variables",ylab="Adjusted R-Squared",type="l")
which.max(Best_Subset_3_Summary$adjr2)
points(4,Best_Subset_3_Summary$adjr2[4], col="red",cex=2,pch=20)
plot(Best_Subset_3_Summary$cp,xlab="Number of Variables",ylab="Cp",type='l')
which.min(Best_Subset_3_Summary$cp)
points(1,Best_Subset_3_Summary$cp[1],col="red",cex=2,pch=20)
which.min(Best_Subset_3_Summary$bic)
plot(Best_Subset_3_Summary$bic,xlab="Number of Variables",ylab="BIC",type='l')
points(4,Best_Subset_3_Summary$bic[4],col="red",cex=2,pch=20)
# BIC and Adjusted R-squared suggested p = 4
# Cp suggested p = 1.

# The two models are not nested, so use cv to determine which one is better.
# Defining training control as cross-validation (K=5).
set.seed(08092021)
Train_Control <- trainControl(method = "cv", number = 5)
ICS3_Mod_p1_CV <- train(Diff_6_24 ~ 1 + Rep + sodCI, data = ICS_Diff_6_24_Selected_2, method = "lm", trControl = Train_Control)
print(ICS3_Mod_p1_CV)
ICS3_Mod_p4_CV <- train(Diff_6_24 ~ 1 + Rep + Epi_Type + gtgE_sseI + ygcH_stfC + bglA_bepF + avrA, data = ICS_Diff_6_24_Selected_2, method = "lm", trControl = Train_Control)
print(ICS3_Mod_p4_CV)

# The performances were similar.
# Go with p=1 to avoid overfitting.
##-----------------------------------------------------------
# Build the ICS-3 model.
ICS3_Mod_1 <- lm(Diff_6_24 ~ 1 + Rep + sodCI, data = ICS_Diff_6_24_Selected_2)
# One-way ANOVA for the ICS-3 model.
anova(ICS3_Mod_1)
##-----------------------------------------------------------
# Model diagnostics.
par(mfrow = c(2, 1), mar=c(4,4,1,1))

# Check linearity and homoscedasticity
plot(fitted(ICS3_Mod_1),residuals(ICS3_Mod_1), xlab = "Fitted Values", ylab = "Residuals")
# Check the residual normality
qqnorm(residuals(ICS3_Mod_1), ylab="Residuals", plot.it=TRUE, datax=FALSE, main = NULL)
qqline(resid(ICS3_Mod_1))
par(mfrow = c(1, 1))
##-----------------------------------------------------------
# Post-hoc analyses.
# sodCI.
ICS3_Mod_1_sodCI.emm <- emmeans(ICS3_Mod_1, ~sodCI)
ICS3_Mod_1_sodCI.ctr <- contrast(ICS3_Mod_1_sodCI.emm, interaction = "pairwise", simple = "sodCI", combine = TRUE, adjust = "tukey")
##-----------------------------------------------------------
