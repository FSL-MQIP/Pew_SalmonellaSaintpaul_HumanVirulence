##-----------------------------------------------------------
## Project: PEW trust
## Script purpose: Statistical Analysis for the invasion assay.

## Start Date:  May 10, 2021
##-----------------------------------------------------------
## Notes: 
## (1)
## (2) 
## (3) 

### Load packages.
##-----------------------------------------------------------
#library(ordPens) # Cannot be loaded successfully.
library(lme4);library(lmerTest);library(lattice);library(emmeans)
library(MASS);library(lmtest)
library(plyr);library(tidyverse);library(stringi);library(stringr)
library(ggeffects);library(ggplot2);library(ggpubr);library(multcomp)
##-----------------------------------------------------------
PathToFile1 <- "~/OneDrive - Cornell University/Cornell document/Food safety lab/PEW_project/Phenotypic_Characterization/PhenotypicCharacterization_Results/InvasionAssay/ProcessedData/InvasionAssayData_ForStatisticalAnalysis_RC051021.csv"
InputFile1 <- read.csv(PathToFile1)
InputFile1$Plating_Order <- as.factor(InputFile1$Plating_Order)
InputFile1$Replicate <- as.factor(InputFile1$Replicate)
##-----------------------------------------------------------
# Create a plot showing the invasion assay results.

# Summarize the data (mean + se across bio reps).
InvasionAssay_DataSummary <- ddply(InputFile1, c("Strain"), summarise, 
                                   N = sum(!is.na(LOG10_N_N0)),
                                   mean = mean(LOG10_N_N0, na.rm = TRUE),
                                   sd = sd(LOG10_N_N0, na.rm = TRUE),
                                   se = sd/sqrt(N))
InvasionAssay_DataSummary$Epi_Type <- c("CN", rep("NHA", 3), "CP", rep("HA", 4))
levels(InvasionAssay_DataSummary$Epi_Type) <- c("CN", "CP", "HA", "NHA")

# Create the plot.
InvasionAssay_P <- ggplot(data = InputFile1) + 
  geom_pointrange(data = InvasionAssay_DataSummary, aes(x=Strain, y=mean, ymin=mean-se, ymax=mean+se, fill = Epi_Type),shape=23,size=1,color="black") +
  geom_point(stat = "identity", aes(x=Strain, y=LOG10_N_N0, group=Replicate, shape = Replicate),
             position = position_dodge(0.2), size=2, fill = "grey10", alpha=.9) +
  scale_fill_manual(values=c("springgreen", "orange", "red", "deepskyblue2"), name="Epi-Type",
                    breaks=c("CN", "CP", "HA", "NHA"),
                    labels=c("Negative Control", "Positive Control", "Human-associated", "Non-human-associated")) +
  scale_shape_manual(values=c(3, 7, 8, 11), name = "Biological Replicate", breaks=c("1", "2", "3", "4"), labels=c("Rep 1", "Rep 2", "Rep 3", "Rep 4"))+
  labs(x="Strains", y="Log10(N/N0)") +
  theme(panel.background = element_rect(fill = "grey93", colour = "grey93", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white")) +
  theme(legend.position = "bottom", legend.margin = margin(0.5, 0.5, 0.5, 0.5), legend.text = element_text(size = 10), legend.title = element_text(size = 10))

ggsave("./InvasionAssay_RawResults_RC042921.pdf", InvasionAssay_P, device = "pdf",height = 6, width = 12)

##-----------------------------------------------------------
# Check whether Plating_Order should be included as a random effect.
Mod0.1 <- lmer(LOG10_N_N0 ~ 1 + (1|Plating_Order), data = InputFile1, REML = FALSE)
summary(Mod0.1)
# ICC = 0.005952275
# We don't need to include plating order as a random effect.
Mod0.2 <- lm(LOG10_N_N0 ~ 1, data = InputFile1)
summary(Mod0.2)
##-----------------------------------------------------------
# Epi_Type as the primary variable of interest.
Mod1.1 <- lm(LOG10_N_N0 ~ Epi_Type + Replicate, data = InputFile1)
summary(Mod1.1)
anova(Mod1.1)
anova(Mod0.2, Mod1.1)

Mod1.2 <- lm(LOG10_N_N0 ~ Epi_Type + Replicate + Lineage, data = InputFile1)
summary(Mod1.2)
anova(Mod1.2)
anova(Mod1.1, Mod1.2)

Mod1.3 <- lm(LOG10_N_N0 ~ Epi_Type*Lineage + Replicate, data = InputFile1)
summary(Mod1.3)
anova(Mod1.3)
anova(Mod1.2, Mod1.3)
anova(Mod0.2, Mod1.1, Mod1.2, Mod1.3)
AIC(Mod0.2);AIC(Mod1.1);AIC(Mod1.2);AIC(Mod1.3)

# Mod1.1 should be preferred, because of the one-to-one relationship between Epi_TypeCN and LineageTyphimurium_Negative 
# and between Epi_TypeCP and LineageTyphimurium_Positive.

coef(Mod1.1)
Mod1.1$coefficients
Mod1.1$fitted.values

confint(Mod1.1) #confidence interval for regression coefficient
predict(Mod1.1,data.frame(Epi_Type=(c("CN", "CP", "HA", "NHA")), Replicate=(c("1", "2", "3", "4"))), interval="confidence")  #CI for fitted values
predict(Mod1.1,data.frame(Epi_Type=(c("CN", "CP", "HA", "NHA")), Replicate=(c("1", "2", "3", "4"))), interval="prediction")  # prediction interval for Y


par(mfrow = c(2, 1))
# Scatterplot to check linearity and homoscedasticity
plot(fitted(Mod1.1),residuals(Mod1.1))
# QQ plot to check the residual normality
qqnorm(residuals(Mod1.1), ylab="Residuals", plot.it=TRUE, datax=FALSE)
qqline(resid(Mod1.1))
par(mfrow = c(1, 1))


# Pairwise comparison of the invasion efficiency across different levels of Epi_Type.
Mod1.1.emm <- emmeans(Mod1.1, ~Epi_Type)
Mod1.1.ctr1 <- contrast(Mod1.1.emm, interaction = "pairwise", simple = "Epi_Type", combine = TRUE, adjust = "tukey")
Mod1.1.cld1 <- cld(Mod1.1.emm, single = c("Epi_Type"), Letters = LETTERS)
Mod1.1.cld1$Cld_Letters <- trimws(Mod1.1.cld1$.group)

InvasionAssay_DataSummary_ET <- ddply(InputFile1, c("Epi_Type"), summarise, 
                                   N = sum(!is.na(LOG10_N_N0)),
                                   mean = mean(LOG10_N_N0, na.rm = TRUE),
                                   sd = sd(LOG10_N_N0, na.rm = TRUE),
                                   se = sd/sqrt(N))

InvasionAssay_ModResults_InvET_P <- ggplot(data = Mod1.1.cld1) + 
  geom_point(stat = "identity", aes(x=Epi_Type, y=emmean, color = Epi_Type), alpha = 0.7, size = 6) +
  geom_pointrange(data = InvasionAssay_DataSummary_ET, aes(x=Epi_Type, y=mean, ymin=mean-se, ymax=mean+se, fill = Epi_Type),shape=23,size=0.5,color="black", alpha = 0.9) +
  geom_point(data = InputFile1, aes(x=Epi_Type, y=LOG10_N_N0, group=Replicate, shape = Replicate), position = position_dodge(0.2), size=1.5) +
  geom_text(aes(x=Epi_Type, y=emmean, label=Cld_Letters, color = Epi_Type), vjust=2, hjust=1, size = 5, position = position_dodge(.2), fontface = "bold") +
  scale_fill_manual(values=c("black", "orange", "red", "deepskyblue2"), name="Epi-Type",
                    breaks=c("CN", "CP", "HA", "NHA"),
                    labels=c("Negative Control", "Positive Control", "Human-associated", "Non-human-associated")) +
  scale_color_manual(values=c("black", "orange", "red", "deepskyblue2"), name="Epi-Type",
                     breaks=c("CN", "CP", "HA", "NHA"),
                     labels=c("Negative Control", "Positive Control", "Human-associated", "Non-human-associated")) +
  scale_shape_manual(values=c(3, 7, 8, 11), name = "Biological Replicate", breaks=c("1", "2", "3", "4"), labels=c("Rep 1", "Rep 2", "Rep 3", "Rep 4"))+
  labs(x="Epi-Type", y="Log10(N/N0)") +
  theme(panel.background = element_rect(fill = "grey93", colour = "grey93", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white")) +
  theme(legend.position = "bottom", legend.margin = margin(0.5, 0.5, 0.5, 0.5), legend.text = element_text(size = 10), legend.title = element_text(size = 10))

ggsave("./InvasionAssay_ModResults_InvET_RC051821.pdf", InvasionAssay_ModResults_InvET_P, device = "pdf",height = 6, width = 12)



##-----------------------------------------------------------
# Strain as the primary variable of interest.
Mod2 <- lm(LOG10_N_N0 ~ Strain + Replicate, data = InputFile1)
summary(Mod2)
anova(Mod2)
anova(Mod0.2, Mod2)
AIC(Mod1.1);AIC(Mod2)
# Mod2 should be preferred.

par(mfrow = c(2, 1), mar=c(4,4,1,1))
# The margin indices go from bottom - left - top - right.
# Scatterplot to check linearity and homoscedasticity
plot(fitted(Mod2),residuals(Mod2), xlab = "Fitted Values", ylab = "Residuals")
# QQ plot to check the residual normality
qqnorm(residuals(Mod2), ylab="Residuals", plot.it=TRUE, datax=FALSE, main = NULL)
qqline(resid(Mod2))
par(mfrow = c(1, 1))



# Identify outliers based on the residuals.
InputFile1$Fitted_Value <- as.double(fitted(Mod2))
InputFile1$Residual <- as.double(residuals(Mod2))

# Plot the residual to see if there are outliers.
ResidualBoxPlot <- ggplot(InputFile1) +
  aes(x = "", y = Residual) +
  geom_boxplot(fill = "red", alpha = 0.6) +
  theme_minimal()

# Identify the outliers in the dataset.
Outliers <- boxplot.stats(InputFile1$Residual)$out
OutliersRows <- which(InputFile1$Residual %in% c(Outliers))
OutliersRows
# There are no outliers in terms of the residuals.
# Therefore, the data points for rep 3 are retained.

# Pairwise comparison of the invasion efficiency across different levels of Epi_Type.
Mod2.emm <- emmeans(Mod2, ~Strain)
Mod2.ctr1 <- contrast(Mod2.emm, interaction = "pairwise", simple = "Strain", combine = TRUE, adjust = "tukey")
Mod2.cld1 <- cld(Mod2.emm, single = c("Strain"), Letters = LETTERS)
Mod2.cld1$Cld_Letters <- trimws(Mod2.cld1$.group)
Mod2.cld1$Epi_Type <- c("CN", "HA", "NHA", "HA", "CP", "HA", "HA", "NHA", "NHA")
Mod2.cld1$Epi_Type <- as.factor(Mod2.cld1$Epi_Type)
InvasionAssay_ModResults_P <- ggplot(data = Mod2.cld1) + 
  geom_point(stat = "identity", aes(x=Strain, y=emmean, color = Epi_Type, fill = Epi_Type), alpha = 0.7, size = 6) +
  geom_pointrange(data = InvasionAssay_DataSummary, aes(x=Strain, y=mean, ymin=mean-se, ymax=mean+se, fill = Epi_Type),shape=23,size=0.5,color="black", alpha = 0.9) +
  geom_point(data = InputFile1, aes(x=Strain, y=LOG10_N_N0, group=Replicate, shape = Replicate),
             position = position_dodge(0.2), size=3, fill = "grey10", color = "black") +
  geom_text(aes(x=Strain, y=emmean, label=Cld_Letters, color = Epi_Type), vjust=2, hjust=1, size = 5, position = position_dodge(.2), fontface = "bold") +
  scale_fill_manual(values=c("black", "orange", "red", "deepskyblue2"), name="Epi-Type",
                    breaks=c("CN", "CP", "HA", "NHA"),
                    labels=c("Negative Control", "Positive Control", "Human-associated", "Non-human-associated")) +
  scale_color_manual(values=c("black", "orange", "red", "deepskyblue2"), name="Epi-Type",
                    breaks=c("CN", "CP", "HA", "NHA"),
                    labels=c("Negative Control", "Positive Control", "Human-associated", "Non-human-associated")) +
  scale_shape_manual(values=c(3, 7, 8, 11), name = "Biological Replicate", breaks=c("1", "2", "3", "4"), labels=c("Rep 1", "Rep 2", "Rep 3", "Rep 4"))+
  labs(x="Strains", y="Log10(N/N0)") +
  theme(panel.background = element_rect(fill = "grey93", colour = "grey93", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white")) +
  theme(legend.position = "bottom", legend.margin = margin(0.5, 0.5, 0.5, 0.5), legend.text = element_text(size = 10), legend.title = element_text(size = 10))

ggsave("./InvasionAssay_ModResults_RC051821.pdf", InvasionAssay_ModResults_P, device = "pdf",height = 6, width = 12)

##-----------------------------------------------------------




