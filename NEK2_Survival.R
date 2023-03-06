# CONTENTS OF Code
################################################################
# PROJECT 1: UAMS Survival Data and NEK2 and CD274 (n1136) - (Line 47 - Line 917)
################################################################
   ### SECTION 1: Load Packages
   ### SECTION 2:  Import UAMS data and create Expression Set for Survival Data
      #### Subsection 2A: Finding Event Free Survival (EFS) Cut-point for NEK2 and CD274 
          ##### Sub-subsection 2AI:  Kaplan-Meier(EFS) for NEK2 (204641_at)
          ##### Sub-subsection 2AII: Kaplan-Meier(EFS) for CD274 (223834_at)
      #### Subsection 2B: Finding Overall Survival (OS) Cut-point for NEK2 and CD274 
          ##### Sub-subsection 2BI:  Kaplan-Meier(OS) for NEK2 (204641_at)
          ##### Sub-subsection 2BII: Kaplan-Meier(OS) for CD274 (223834_at)
   ### SECTION 3 - Hazard Ratios for NEK2 and CD274 (EFS)
   ### SECTION 4 - Hazard Ratios for NEK2 and CD274 (OS)
   ### SECTION 5 - Import and prepare NEK2_CD274 SUBGROUP combinations for EFS analysis
      #### Subsection 5A: Hazard Ratio for NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi (EFS)
      #### Subsection 5B: Hazard Ratio for NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi (EFS)
      #### Subsection 5C: Hazard Ratio for NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi (EFS)
   ### SECTION 6 - Import and prepare NEK2_CD274 SUBGROUP combinations for OS analysis
      #### Subsection 6A: Hazard Ratio for NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi (OS)
      #### Subsection 6B: Hazard Ratio for NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi (OS)
      #### Subsection 6C: Hazard Ratio for NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi (OS)

#################################################################
# PROJECT 2: MMRF Survival Data and NEK2 and CD274 (n592) - (Line 921 - 1795)
#################################################################
    ### SECTION 1: Load Packages
    ### SECTION 2:  Import MMRF data and create Expression Set for Survival Data
        #### Subsection 2A: Finding Event Free Survival (EFS) Cut-point for NEK2 and CD274 
            ##### Sub-subsection 2AI:  Kaplan-Meier(EFS) for NEK2 (ENSG00000117650)
            ##### Sub-subsection 2AII: Kaplan-Meier(EFS) for CD274 (ENSG00000120217)
        #### Subsection 2B: Finding Overall Survival (OS) Cut-point for NEK2 and CD274 
            ##### Sub-subsection 2BI:  Kaplan-Meier(OS) for NEK2 (ENSG00000117650)
            ##### Sub-subsection 2BII: Kaplan-Meier(OS) for CD274 (ENSG00000120217)
    ### SECTION 3 - Hazard Ratios for NEK2 and CD274 (EFS)
    ### SECTION 4 - Hazard Ratios for NEK2 and CD274 (OS)
    ### SECTION 5 - Import and prepare NEK2_CD274 SUBGROUP combinations for EFS analysis
       #### Subsection 5A: Hazard Ratio for NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi (EFS)
       #### Subsection 5B: Hazard Ratio for NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi (EFS)
       #### Subsection 5C: Hazard Ratio for NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi (EFS)
    ### SECTION 6 - Import and prepare NEK2_CD274 SUBGROUP combinations for OS analysis
       #### Subsection 6A: Hazard Ratio for NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi (OS)
       #### Subsection 6B: Hazard Ratio for NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi (OS)
       #### Subsection 6C: Hazard Ratio for NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi (OS)  

##################################################################################
##################################################################################
## PROJECT 1: UAMS Survival Data (n1136)
##################################################################################
##################################################################################

##################################################################################
### Section 1: Load Packages
##################################################################################
# if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("Biobase")
# BiocManager::install("GEOquery")
# BiocManager::install("Biobase")
# install.packages(c("survival", "survminer"))
# install.packages("tidyverse")
# install.packages("readxl")
# install.packages("bquote")
library(survival)
library(survminer)
library(ranger)
library(ggplot2)
library(ggfortify)
library(scales)
library(FSA)
library(rcompanion)
library(rstatix)
library(ggpubr)
library(tidyverse)
library(GEOquery)
library(gtsummary)
library(dplyr)
library(Biobase)
library(tibble)
library(readxl)
###################################################################################
### SECTION 2:  Import UAMS data and create Expression Set
###################################################################################
x <- read.table("Assay_NEK2_CD274_UAMS.txt", sep="\t", header = T, row.names = 1)
class(x) # Data frame
x <- as.matrix(x)
class(x) # Matrix array

# Feature Data
f <- read.table("Feature_NEK2_CD274_UAMS.txt", sep="\t", header =T, row.names = 1)
dim(f)
class(f)

# Phenotype data
p <- read.table("Phenotype_NEK2_CD274_UAMS.txt", sep="\t", header =T, row.names = 1)
dim(p)
class(p)

# Create Expression Set
eset <- ExpressionSet(assayData = x,
                      phenoData = AnnotatedDataFrame(p),
                      featureData = AnnotatedDataFrame(f))

#View the number of features/Proteins (rows) and samples (columns)
dim(eset)
# Access data from an ExpressionSet object
x <- exprs(eset) # you can retreive the expression matrix with the function `exprs`
f <- fData(eset) # Retrieve the feature data with `fData'
p <- pData(eset)# Retrieve the phenotype data with`pData`

# extract information of interest from the phenotype data (pdata)
idx <- which(colnames(pData(eset)) %in%
               c('CensOS','CensEFS', 'MonthsOS', 'MonthsEFS'))

metadata <- data.frame(pData(eset)[,idx],
                       row.names = rownames(pData(eset)))

# check that sample names match exactly between pdata and Z-scores 
all((colnames(x) == rownames(metadata)) == TRUE)

# create a merged pdata and Z-scores object
coxdata <- data.frame(metadata, t(x))

coxdata_df <- tibble::rownames_to_column(coxdata, "CHIPID")

# prepare phenotypes
coxdata$CensOS <- as.numeric(coxdata$CensOS) # Censor for the data
coxdata$MonthsOS <- as.numeric(gsub('^KJX|^KJ', '', coxdata$MonthsOS)) # Similar to DaysOS
coxdata$CensEFS <- as.numeric(coxdata$CensEFS) # Censor for the data
coxdata$MonthsEFS <- as.numeric(gsub('^KJX|^KJ', '', coxdata$MonthsEFS)) # Similar to DaysOS

############################################################################################
#### Subsection A: Finding Event Free Survival (EFS) Cut-point for NEK2 and CD274 
############################################################################################
# surv_cutpoint(): Determine the optimal cut-point for each variable using 'maxstat'.
res.cut_EFS <-surv_cutpoint(
  coxdata,
  time = "MonthsEFS",
  event = "CensEFS",
  variables = colnames(coxdata)[5:ncol(coxdata)],
  minprop = 0.10,
  progressbar = TRUE
)
options(max.print = 5000)        # Change global options
summary(res.cut_EFS)

# RES.CAT (categorize by High and Low expression from cutpoints)
res.cat_EFS <- surv_categorize(res.cut_EFS) ### FINE
head(res.cat_EFS)
summary(res.cat_EFS)
str(res.cat_EFS)
#################################################################################
##### Sub-subsection I: Kaplan-Meier(EFS) for NEK2 (204641_at)
#################################################################################
# 4. Fit survival curves and visualize
{fit_NEK2_EFS <- survfit(Surv(MonthsEFS, CensEFS) ~NEK2_AFFY_204641_at, data = res.cat_EFS)
  a <- surv_pvalue(fit_NEK2_EFS)$pval
  b <- unname(summary(fit_NEK2_EFS)$table[,'records'])
  print(a)
  print(b)
  print(fit_NEK2_EFS)
  ggsurvplot(fit_NEK2_EFS, data = res.cat_EFS, risk.table = TRUE, conf.int = F, pval = TRUE)
}

ggsurv_NEK2_EFS <- ggsurvplot(fit_NEK2_EFS, data = res.cat_EFS, 
                         break.time.by = 25,
                               conf.int=TRUE, pval=TRUE, risk.table = TRUE, 
                               surv.median.line = "hv", # Specify median survival
                               tables.theme = theme_survminer(
                                 font.main = c(20, "bold", "darkblue"),
                                 font.submain = c(15, "bold.italic", "purple"),
                                 font.caption = c(14, "plain", "oCBX3ge"),
                                 font.x = c(20, "bold.italic", "red"),
                                 font.y = c(15, "bold.italic", "darkred"),
                                 font.tickslab = c(16, "plain", "darkgreen")),
                               legend.labs=c("High Expression","Low Expression"), legend.title="NEK2",  
                               palette=c("red3","blue"), xlab = "Time (Months)", 
                               title="          NEK2 (204641_at) EFS (UAMS)",
                               caption = "created with survminer",
                               font.title = c(30, "bold", "darkblue"),
                               font.subtitle = c(18, "bold.italic", "purple"),
                               font.legend = c(18),
                               font.x = c(18, "bold.italic", "red"),
                               font.y = c(18, "bold.italic", "darkred"),
                               font.tickslab = c(16, "plain", "darkgreen"),
                               risk.table.fontsize = 6,
                               ########## risk table #########,
                               risk.table.height = 0.25)

ggsurv_NEK2_EFS
#################################################################################
##### Sub-subsection II: Kaplan-Meier(EFS) for CD274 (223834_at)
#################################################################################
# Fit survival curves and visualize
{fit_CD274_EFS <- survfit(Surv(MonthsEFS, CensEFS) ~CD274_AFFY_223834_at, data = res.cat_EFS)
a <- surv_pvalue(fit_CD274_EFS)$pval
b <- unname(summary(fit_CD274_EFS)$table[,'records'])
print(a)
print(b)
print(fit_CD274_EFS)
ggsurvplot(fit_CD274_EFS, data = res.cat_EFS, risk.table = TRUE, conf.int = F, pval = TRUE)
}

ggsurv_CD274_EFS <- ggsurvplot(fit_CD274_EFS, data = res.cat_EFS, 
                         break.time.by = 25,
                         conf.int=TRUE, pval=TRUE, risk.table = TRUE, 
                         surv.median.line = "hv", # Specify median survival
                         tables.theme = theme_survminer(
                           font.main = c(20, "bold", "darkblue"),
                           font.submain = c(15, "bold.italic", "purple"),
                           font.caption = c(14, "plain", "oCBX3ge"),
                           font.x = c(20, "bold.italic", "red"),
                           font.y = c(15, "bold.italic", "darkred"),
                           font.tickslab = c(16, "plain", "darkgreen")),
                         legend.labs=c("High Expression","Low Expression"), legend.title="CD274",  
                         palette=c("red3","blue"), xlab = "Time (Months)", 
                         title="          CD274 (223834_at) EFS (UAMS)",
                         caption = "created with survminer",
                         font.title = c(30, "bold", "darkblue"),
                         font.subtitle = c(18, "bold.italic", "purple"),
                         font.legend = c(18),
                         font.x = c(18, "bold.italic", "red"),
                         font.y = c(18, "bold.italic", "darkred"),
                         font.tickslab = c(16, "plain", "darkgreen"),
                         risk.table.fontsize = 6,
                         ########## risk table #########,
                         risk.table.height = 0.25)
ggsurv_CD274_EFS
################################################################################
#### Subsection B: Finding Overall Survival (OS) Cut-point for NEK2 and CD274
################################################################################
# surv_cutpoint(): Determine the optimal cut-point for each variable using 'maxstat'.
res.cut_OS <-surv_cutpoint(
  coxdata,
  time = "MonthsOS",
  event = "CensOS",
  variables = colnames(coxdata)[5:ncol(coxdata)],
  minprop = 0.1,
  progressbar = TRUE
)
options(max.print = 5000)        # Change global options
summary(res.cut_OS)

# RES.CAT (categorize by High and Low expression from cutpoints)
res.cat_OS <- surv_categorize(res.cut_OS) ### FINE
head(res.cat_OS)
summary(res.cat_OS)
str(res.cat_OS)
#################################################################################
##### Sub-subsection I:  Kaplan-Meier(OS) for NEK2 (204641_at)
#################################################################################
# Fit survival curves and visualize
{fit_NEK2_OS <- survfit(Surv(MonthsOS, CensOS) ~NEK2_AFFY_204641_at, data = res.cat_OS)
a <- surv_pvalue(fit_NEK2_OS)$pval
b <- unname(summary(fit_NEK2_OS)$table[,'records'])
print(a)
print(b)
print(fit_NEK2_OS)
ggsurvplot(fit_NEK2_OS, data = res.cat_OS, risk.table = TRUE, conf.int = F, pval = TRUE)
}

ggsurv_NEK2_OS <- ggsurvplot(fit_NEK2_OS, data = res.cat_OS, 
                         break.time.by = 25,
                         conf.int=TRUE, pval=TRUE, risk.table = TRUE, 
                         surv.median.line = "hv", # Specify median survival
                         tables.theme = theme_survminer(
                           font.main = c(20, "bold", "darkblue"),
                           font.submain = c(15, "bold.italic", "purple"),
                           font.caption = c(14, "plain", "oCBX3ge"),
                           font.x = c(20, "bold.italic", "red"),
                           font.y = c(15, "bold.italic", "darkred"),
                           font.tickslab = c(16, "plain", "darkgreen")),
                         legend.labs=c("High Expression","Low Expression"), legend.title="NEK2",  
                         palette=c("red3","blue"), xlab = "Time (Months)", 
                         title="          NEK2 (204641_at) OS (UAMS)",
                         caption = "created with survminer",
                         font.title = c(30, "bold", "darkblue"),
                         font.subtitle = c(18, "bold.italic", "purple"),
                         font.legend = c(18),
                         font.x = c(18, "bold.italic", "red"),
                         font.y = c(18, "bold.italic", "darkred"),
                         font.tickslab = c(16, "plain", "darkgreen"),
                         risk.table.fontsize = 6,
                         ########## risk table #########,
                         risk.table.height = 0.25)
ggsurv_NEK2_OS
#################################################################################
##### Sub-subsection II: Kaplan-Meier(OS) for CD274 (223834_at)
#################################################################################
# Fit survival curves and visualize
{fit_CD274_OS <- survfit(Surv(MonthsOS, CensOS) ~CD274_AFFY_223834_at, data = res.cat_OS)
a <- surv_pvalue(fit_CD274_OS)$pval
b <- unname(summary(fit_CD274_OS)$table[,'records'])
print(a)
print(b)
print(fit_CD274_OS)
ggsurvplot(fit_CD274_OS, data = res.cat_OS, risk.table = TRUE, conf.int = F, pval = TRUE)
}

ggsurv_CD274_OS <- ggsurvplot(fit_CD274_OS, data = res.cat_OS, 
                         break.time.by = 25,
                         conf.int=TRUE, pval=TRUE, risk.table = TRUE, 
                         surv.median.line = "hv", # Specify median survival
                         tables.theme = theme_survminer(
                           font.main = c(20, "bold", "darkblue"),
                           font.submain = c(15, "bold.italic", "purple"),
                           font.caption = c(14, "plain", "oCBX3ge"),
                           font.x = c(20, "bold.italic", "red"),
                           font.y = c(15, "bold.italic", "darkred"),
                           font.tickslab = c(16, "plain", "darkgreen")),
                         legend.labs=c("High Expression","Low Expression"), legend.title="CD274",  
                         palette=c("red3","blue"), xlab = "Time (Months)", 
                         title="          CD274 (223834_at) OS (UAMS)",
                         caption = "created with survminer",
                         font.title = c(30, "bold", "darkblue"),
                         font.subtitle = c(18, "bold.italic", "purple"),
                         font.legend = c(18),
                         font.x = c(18, "bold.italic", "red"),
                         font.y = c(18, "bold.italic", "darkred"),
                         font.tickslab = c(16, "plain", "darkgreen"),
                         risk.table.fontsize = 6,
                         ########## risk table #########,
                         risk.table.height = 0.25)
ggsurv_CD274_OS
################################################################################################
### SECTION 3 - Hazard Ratios for NEK2 and CD274 (EFS)
################################################################################################
# Import EFS UAMS Data
NEK2_CD274_EFS <- read_excel("NEK2_CD274_SUBGROUPS_UAMS_EFS.xlsx")
NEK2_CD274_EFS
NEK2_CD274_EFS <- as.data.frame(NEK2_CD274_EFS)

## To apply the univariate coxph function to multiple covariates at once for EFS, type this:
covariates_EFS <- c("NEK2", "CD274")
univ_formulas_EFS <- sapply(covariates_EFS,
                            function(x) as.formula(paste('Surv(MonthsEFS, CensEFS)~',x)))

univ_models_EFS <- lapply(univ_formulas_EFS, function(x){coxph(x, data = NEK2_CD274_EFS)})
#Extract data
univ_results_EFS <- lapply(univ_models_EFS,
                           function(x){
                             x <- summary(x)
                             p.value <- signif(x$wald["pvalue"], digits = 3)
                             wald.test <- signif(x$wald["test"],digits = 2)
                             beta <- signif(x$coef[1],digits = 3); #coefficient beta
                             HR <- signif(x$coef[2],digits = 3); #exp(beta)
                             HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                             HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                             HR <- paste0(HR, " (",
                                          HR.confint.lower, "-",HR.confint.upper, ")")
                             res <- c(beta, HR, wald.test, p.value)
                             names(res) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                             return(res)
                             #return(exp(cbind(coef(x), confint(x))))
                           })
res_EFS <- t(as.data.frame(univ_results_EFS, check.names = FALSE))
as.data.frame(res_EFS)

# Multivariate Cox regression analysis
res_cox_EFS <- coxph(Surv(MonthsEFS, CensEFS) ~ 
                       NEK2 + CD274, data = NEK2_CD274_EFS)
summary(res_cox_EFS)
################################################################################################
### SECTION 4 - Hazard Ratios for NEK2 and CD274 (OS)
################################################################################################
# Import OS UAMS Data
NEK2_CD274_OS <- read_excel("NEK2_CD274_SUBGROUPS_UAMS_OS.xlsx")
NEK2_CD274_OS
NEK2_CD274_OS <- as.data.frame(NEK2_CD274_OS)

## To apply the univariate coxph function to multiple covariates at once for OS, type this:
covariates_OS <- c("NEK2", "CD274")
univ_formulas_OS <- sapply(covariates_OS,
                           function(x) as.formula(paste('Surv(MonthsOS, CensOS)~',x)))

univ_models_OS <- lapply(univ_formulas_OS, function(x){coxph(x, data = NEK2_CD274_OS)})
#Extract data
univ_results_OS <- lapply(univ_models_OS,
                          function(x){
                            x <- summary(x)
                            p.value <- signif(x$wald["pvalue"], digits = 3)
                            wald.test <- signif(x$wald["test"],digits = 2)
                            beta <- signif(x$coef[1],digits = 3); #coefficient beta
                            HR <- signif(x$coef[2],digits = 3); #exp(beta)
                            HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                            HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                            HR <- paste0(HR, " (",
                                         HR.confint.lower, "-",HR.confint.upper, ")")
                            res <- c(beta, HR, wald.test, p.value)
                            names(res) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                            return(res)
                            #return(exp(cbind(coef(x), confint(x))))
                          })
res_OS <- t(as.data.frame(univ_results_OS, check.names = FALSE))
as.data.frame(res_OS)

# Multivariate Cox regression analysis
res_cox_OS <- coxph(Surv(MonthsOS, CensOS) ~ 
                      NEK2 + CD274, data = NEK2_CD274_OS)
summary(res_cox_OS)
############################################################################################
### SECTION 5 - Import and prepare NEK2_CD274 SUBGROUP combinations for EFS analysis
#############################################################################################
# Import and prepare data
Subgroups <- read_excel("NEK2_CD274_UAMS_EFS.xlsx")
Subgroups
class(Subgroups)
Subgroups <- as.data.frame(Subgroups)
Subgroups.cox <- coxph(Surv(MonthsEFS, CensEFS) ~ Lo_Hi + Hi_Lo + Hi_Hi, data = Subgroups)
summary(Subgroups.cox)

# Median OS and 95% CI for 4 subgroups
Summarize(MonthsEFS ~ NEK2_CD274, #This needs to be name of column in excel
          data=Subgroups,
          digits=3)
groupwiseMedian(MonthsEFS ~ NEK2_CD274,
                data       = Subgroups,
                conf       = 0.95,
                R          = 5000,
                percentile = TRUE,
                bca        = FALSE,
                digits     = 4)
# Histogram and survival curves for  4 subgroups
class(Subgroups)
Subgroups <- as.data.frame(Subgroups)

# Reorder the Columns first
Risk_order <- c("Lo_Lo", "Lo_Hi", "Hi_Lo", "Hi_Hi") #REORDER the columns!
ggplot(Subgroups, aes(x = factor(NEK2_CD274, Risk_order))) +  # Put Risk_order in to help
  geom_bar(aes(y = (..count..)/sum(..count..))) +
  scale_y_continuous(labels=percent) +
  labs(title = 'Distribution (Using MM Date)', x = 'Score Range', y = 'Percentage %') +
  theme_minimal() +
  theme(axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        text = element_text(size=25))
#Next, we look at survival curves by subgroups.
km_trt_fit_EFS <- survfit(Surv(MonthsEFS, CensEFS) ~ NEK2_CD274, data = Subgroups)
p <- ggsurvplot(km_trt_fit_EFS, conf.int = FALSE,
                break.time.by = 40,
                risk.table = "nrisk_cumevents",
                surv.median.line = "hv", # Specify median survival
                tables.theme = theme_survminer(
                  font.main = c(20, "bold", "black"),
                  font.submain = c(20, "bold.italic", "black"),
                  font.caption = c(10, "plain", "black"),
                  font.x = c(20, "bold.italic", "black"),
                  font.y = c(15, "bold.italic", "black"),
                  font.tickslab = c(12, "plain", "black")),
                legend.labs=c("NEK2-high/CD274-high", "NEK2-high/CD274-low","NEK2-low/CD274-high","NEK2-low/CD274-low"), legend.title="",  
                palette=c("black","red","brown4","blue"), xlab = "Time (Months)", 
                title="                      NEK2/CD274 EFS (UAMS)",
                caption = "created with survminer",
                font.title = c(25, "bold", "black"),
                font.subtitle = c(16, "bold.italic", "black"),
                font.legend = c(12),
                font.x = c(18, "bold.italic", "black"),
                font.y = c(18, "bold.italic", "black"),
                font.tickslab = c(16, "plain", "black"),
                risk.table.fontsize = 4.5,
                ########## risk table #########,
                risk.table.height = 0.25)
p

# Compare the means of the 4 subgroups using ANOVA test
Subgroups.ANOVA <- Subgroups %>% anova_test(MonthsEFS ~ NEK2_CD274)
Subgroups.ANOVA
# Pairwise T-tests for 4 subgroups
pwc <- Subgroups %>%
  # pairwise_t_test(MonthsEFS ~ NEK2_CD274_EFS, p.adjust.method = "bonferroni")
  pairwise_t_test(MonthsEFS ~ NEK2_CD274)
pwc
# Visualization: box plots with p-values of  4 subgroups
Subgroups$NEK2_CD274 <- factor(Subgroups$NEK2_CD274 , levels=c("Hi_Hi","Hi_Lo", "Lo_Hi", "Lo_Lo")) # Reorder Order!!!
pwc <- pwc %>% add_xy_position(x =  "NEK2_CD274")
r <- ggboxplot(Subgroups, x = "NEK2_CD274", y = "MonthsEFS", ylim = c(0,401), fill = "NEK2_CD274", notch = TRUE,
               palette = c("black","red","brown4","blue")) +
  stat_pvalue_manual(pwc, label = "p", size = 5, tip.length = 0.01, step.increase = 0.08) + #step.increase adjusts distance between signif bars
  labs(
    subtitle = get_test_label(Subgroups.ANOVA, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )
r +
  font("title", size = 15, color = "red", face = "bold.italic")+
  font("subtitle", size = 15, color = "orange")+
  font("caption", size = 15, color = "orange")+
  font("xlab", size = 15, color = "blue")+
  font("ylab", size = 17, color = "#993333")+
  font("xy.text", size = 15, color = "gray2", face = "bold")

###############################################################################################
#### Subsection 5A: Hazard Ratio for NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi (EFS)
###############################################################################################

NEK2Hi_CD274Hi_EFS <- read_excel("NEK2_CD274_UAMS_EFS_NEK2Hi_CD274Hi.xlsx")
NEK2Hi_CD274Hi_EFS
NEK2Hi_CD274Hi_EFS <- as.data.frame(NEK2Hi_CD274Hi_EFS)

## To apply the univariate coxph function to multiple covariates at once, type this:
covariates_EFS_Hi_Hi <- c("NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi")
univ_formulas_EFS_Hi_Hi <- sapply(covariates_EFS_Hi_Hi,
                                  function(x) as.formula(paste('Surv(MonthsEFS, CensEFS)~',x)))

univ_models_EFS_Hi_Hi <- lapply(univ_formulas_EFS_Hi_Hi, function(x){coxph(x, data = NEK2Hi_CD274Hi_EFS)})
#Extract data
univ_results_EFS_Hi_Hi <- lapply(univ_models_EFS_Hi_Hi,
                                 function(x){
                                   x <- summary(x)
                                   p.value <- signif(x$wald["pvalue"], digits = 3)
                                   wald.test <- signif(x$wald["test"],digits = 2)
                                   beta <- signif(x$coef[1],digits = 3); #coefficient beta
                                   HR <- signif(x$coef[2],digits = 3); #exp(beta)
                                   HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                                   HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                                   HR <- paste0(HR, " (",
                                                HR.confint.lower, "-",HR.confint.upper, ")")
                                   res_EFS_Hi_Hi <- c(beta, HR, wald.test, p.value)
                                   names(res_EFS_Hi_Hi) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                                   return(res_EFS_Hi_Hi)
                                   #return(exp(cbind(coef(x), confint(x))))
                                 })
res_EFS_Hi_Hi <- t(as.data.frame(univ_results_EFS_Hi_Hi, check.names = FALSE))
as.data.frame(res_EFS_Hi_Hi)


# Multivariate Cox regression analysis
res.cox_EFS_NEK2Hi_CD274Hi <- coxph(Surv(MonthsEFS, CensEFS) ~ 
                                      NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi, data = NEK2Hi_CD274Hi_EFS)
summary(res.cox_EFS_NEK2Hi_CD274Hi)

# Visualizing the estimated distribution of survival times
# Plot the baseline survival function
ggsurvplot(survfit(res.cox_EFS_NEK2Hi_CD274Hi, data = NEK2Hi_CD274Hi_EFS), palette = "#2E9FDF",
           ggtheme = theme_minimal())
# Create the new data  
NEK2Hi_CD274Hi_EFS_df <- with(res.cox_EFS_NEK2Hi_CD274Hi,
                              data.frame(
                                NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi = c(1,0)
                              ))
NEK2Hi_CD274Hi_EFS_df

# Survival curves
NEK2Hi_CD274Hi_EFS_fit <- survfit(res.cox_EFS_NEK2Hi_CD274Hi, newdata = NEK2Hi_CD274Hi_EFS_df)
ggsurvplot(NEK2Hi_CD274Hi_EFS_fit, data = NEK2Hi_CD274Hi_EFS_df, conf.int = TRUE, legend.labs = c("NEK2Hi_CD274Hi_EFS_df = 0", "NEK2Hi_CD274Hi_EFS_df = 1"),
           ggtheme = theme_minimal())
print(NEK2Hi_CD274Hi_EFS_fit)

###############################################################################################
#### Subsection 5B: Hazard Ratio for NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi (EFS)
###############################################################################################
NEK2Hi_CD274Lo_EFS <- read_excel("NEK2_CD274_UAMS_EFS_NEK2Hi_CD274Lo.xlsx")
NEK2Hi_CD274Lo_EFS
NEK2Hi_CD274Lo_EFS <- as.data.frame(NEK2Hi_CD274Lo_EFS)

## To apply the univariate coxph function to multiple covariates at once, type this:
covariates_EFS_Hi_Lo <- c("NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi")
univ_formulas_EFS_Hi_Lo <- sapply(covariates_EFS_Hi_Lo,
                                  function(x) as.formula(paste('Surv(MonthsEFS, CensEFS)~',x)))

univ_models_EFS_Hi_Lo <- lapply(univ_formulas_EFS_Hi_Lo, function(x){coxph(x, data = NEK2Hi_CD274Lo_EFS)})
#Extract data
univ_results_EFS_Hi_Lo <- lapply(univ_models_EFS_Hi_Lo,
                                 function(x){
                                   x <- summary(x)
                                   p.value <- signif(x$wald["pvalue"], digits = 3)
                                   wald.test <- signif(x$wald["test"],digits = 2)
                                   beta <- signif(x$coef[1],digits = 3); #coefficient beta
                                   HR <- signif(x$coef[2],digits = 3); #exp(beta)
                                   HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                                   HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                                   HR <- paste0(HR, " (",
                                                HR.confint.lower, "-",HR.confint.upper, ")")
                                   res_Hi_Lo <- c(beta, HR, wald.test, p.value)
                                   names(res_Hi_Lo) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                                   return(res_Hi_Lo)
                                   #return(exp(cbind(coef(x), confint(x))))
                                 })
res_EFS_Hi_Lo <- t(as.data.frame(univ_results_EFS_Hi_Lo, check.names = FALSE))
as.data.frame(res_EFS_Hi_Lo)

# Multivariate Cox regression analysis
res.cox_EFS_NEK2Hi_CD274Lo <- coxph(Surv(MonthsEFS, CensEFS) ~ 
                                      NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi, data = NEK2Hi_CD274Lo_EFS)
summary(res.cox_EFS_NEK2Hi_CD274Lo)

# Visualizing the estimated distribution of survival times
# Plot the baseline survival function
ggsurvplot(survfit(res.cox_EFS_NEK2Hi_CD274Lo, data = NEK2Hi_CD274Lo_EFS), palette = "#2E9FDF",
           ggtheme = theme_minimal())
# Create the new data  
NEK2Hi_CD274Lo_EFS_df <- with(res.cox_EFS_NEK2Hi_CD274Lo,
                              data.frame(
                                NEK2 = c(0,0),
                                CD274 = c(0,0),
                                NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi = c(1,0)
                              ))
NEK2Hi_CD274Lo_EFS_df

# Survival curves
NEK2Hi_CD274Lo_EFS_fit <- survfit(res.cox_EFS_NEK2Hi_CD274Lo, newdata = NEK2Hi_CD274Lo_EFS_df)
ggsurvplot(NEK2Hi_CD274Lo_EFS_fit, data = NEK2Hi_CD274Lo_EFS_df, conf.int = TRUE, legend.labs = c("NEK2Hi_CD274Lo_EFS_df = 0", "NEK2Hi_CD274Lo_EFS_df = 1"),
           ggtheme = theme_minimal())
print(NEK2Hi_CD274Lo_EFS_fit)

###############################################################################################
#### Subsection 5C: Hazard Ratio for NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi (EFS)
###############################################################################################
NEK2Lo_CD274Lo_EFS <- read_excel("NEK2_CD274_UAMS_EFS_NEK2Lo_CD274Lo.xlsx")
NEK2Lo_CD274Lo_EFS
NEK2Lo_CD274Lo_EFS <- as.data.frame(NEK2Lo_CD274Lo_EFS)

## To apply the univariate coxph function to multiple covariates at once, type this:
covariates_EFS_Lo_Lo <- c("NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi")
univ_formulas_EFS_Lo_Lo <- sapply(covariates_EFS_Lo_Lo,
                                  function(x) as.formula(paste('Surv(MonthsEFS, CensEFS)~',x)))

univ_models_EFS_Lo_Lo <- lapply(univ_formulas_EFS_Lo_Lo, function(x){coxph(x, data = NEK2Lo_CD274Lo_EFS)})
#Extract data
univ_results_EFS_Lo_Lo <- lapply(univ_models_EFS_Lo_Lo,
                                 function(x){
                                   x <- summary(x)
                                   p.value <- signif(x$wald["pvalue"], digits = 3)
                                   wald.test <- signif(x$wald["test"],digits = 2)
                                   beta <- signif(x$coef[1],digits = 3); #coefficient beta
                                   HR <- signif(x$coef[2],digits = 3); #exp(beta)
                                   HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                                   HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                                   HR <- paste0(HR, " (",
                                                HR.confint.lower, "-",HR.confint.upper, ")")
                                   res_EFS_Lo_Lo <- c(beta, HR, wald.test, p.value)
                                   names(res_EFS_Lo_Lo) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                                   return(res_EFS_Lo_Lo)
                                   #return(exp(cbind(coef(x), confint(x))))
                                 })
res_EFS_Lo_Lo <- t(as.data.frame(univ_results_EFS_Lo_Lo, check.names = FALSE))
as.data.frame(res_EFS_Lo_Lo)

# Multivariate Cox regression analysis
res.cox_EFS_NEK2Lo_CD274Lo <- coxph(Surv(MonthsEFS, CensEFS) ~ 
                                      NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi, data = NEK2Lo_CD274Lo_EFS)
summary(res.cox_EFS_NEK2Lo_CD274Lo)

# Visualizing the estimated distribution of survival times
# Plot the baseline survival function
ggsurvplot(survfit(res.cox_EFS_NEK2Lo_CD274Lo, data = NEK2Lo_CD274Lo_EFS), palette = "#2E9FDF",
           ggtheme = theme_minimal())
# Create the new data  
NEK2Lo_CD274Lo_EFS_df <- with(res.cox_EFS_NEK2Lo_CD274Lo,
                              data.frame(
                                NEK2 = c(0,0),
                                CD274 = c(0,0),
                                NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi = c(1,0)
                              ))
NEK2Lo_CD274Lo_EFS_df

# Survival curves
NEK2Lo_CD274Lo_EFS_fit <- survfit(res.cox_EFS_NEK2Lo_CD274Lo, newdata = NEK2Lo_CD274Lo_EFS_df)
ggsurvplot(NEK2Lo_CD274Lo_EFS_fit, data = NEK2Lo_CD274Lo_EFS_df, conf.int = TRUE, legend.labs = c("NEK2Lo_CD274Lo_EFS_df = 0", "NEK2Lo_CD274Lo_EFS_df = 1"),
           ggtheme = theme_minimal())
print(NEK2Lo_CD274Lo_EFS_fit)

###########################################################################################
### SECTION 6 - Import and prepare NEK2_CD274 SUBGROUP combinations for OS analysis
###########################################################################################
# Multivariate Cox regression analysis of  4 subgroups
Subgroups <- read_excel("NEK2_CD274_UAMS_OS.xlsx")
Subgroups
class(Subgroups)
Subgroups <- as.data.frame(Subgroups)
Subgroups.cox <- coxph(Surv(MonthsOS, CensOS) ~ Lo_Hi + Hi_Lo + Hi_Hi, data = Subgroups)
summary(Subgroups.cox)

# Median OS and 95% CI for  4 subgroups
Summarize(MonthsOS ~ NEK2_CD274, #This needs to be name of column in excel
          data=Subgroups,
          digits=3)
groupwiseMedian(MonthsOS ~ NEK2_CD274,
                data       = Subgroups,
                conf       = 0.95,
                R          = 5000,
                percentile = TRUE,
                bca        = FALSE,
                digits     = 4)
# Histogram and survival curves for  4 subgroups
class(Subgroups)
Subgroups <- as.data.frame(Subgroups)

# Reorder the Columns first
Risk_order <- c("Lo_Lo", "Lo_Hi", "Hi_Lo", "Hi_Hi") #REORDER the columns!
ggplot(Subgroups, aes(x = factor(NEK2_CD274, Risk_order))) +  # Put Risk_order in to help
  geom_bar(aes(y = (..count..)/sum(..count..))) +
  scale_y_continuous(labels=percent) +
  labs(title = 'Distribution (Using MM Date)', x = 'Score Range', y = 'Percentage %') +
  theme_minimal() +
  theme(axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        text = element_text(size=25))
#Next, we look at survival curves by subgroups.
km_trt_fit_CD274_OS <- survfit(Surv(MonthsOS, CensOS) ~ NEK2_CD274, data = Subgroups)
q <- ggsurvplot(km_trt_fit_CD274_OS, conf.int = FALSE,
                break.time.by = 40,
                risk.table = "nrisk_cumevents",
                surv.median.line = "hv", # Specify median survival
                tables.theme = theme_survminer(
                  font.main = c(20, "bold", "black"),
                  font.submain = c(20, "bold.italic", "black"),
                  font.caption = c(10, "plain", "black"),
                  font.x = c(20, "bold.italic", "black"),
                  font.y = c(15, "bold.italic", "black"),
                  font.tickslab = c(12, "plain", "black")),
                legend.labs=c("NEK2-high/CD274-high","NEK2-high/CD274-low","NEK2-low/CD274-high","NEK2-low/CD274-low"), legend.title="",  
                palette=c("black","red","brown4","blue"), xlab = "Time (Months)", 
                title="                      NEK2/CD274 OS (UAMS)",
                caption = "created with survminer",
                font.title = c(25, "bold", "black"),
                font.subtitle = c(16, "bold.italic", "black"),
                font.legend = c(12),
                font.x = c(18, "bold.italic", "black"),
                font.y = c(18, "bold.italic", "black"),
                font.tickslab = c(16, "plain", "black"),
                risk.table.fontsize = 4.5,
                ########## risk table #########,
                risk.table.height = 0.25)
q
# Compare the mean of Subgroups using ANOVA test
Subgroups.ANOVA <- Subgroups %>% anova_test(MonthsOS ~ NEK2_CD274)
Subgroups.ANOVA
# Pairwise T-tests for  4 subgroups
pwc <- Subgroups %>%
  # pairwise_t_test(MonthsOS ~ NEK2_CD274_OS, p.adjust.method = "bonferroni")
  pairwise_t_test(MonthsOS ~ NEK2_CD274)
pwc
# Visualization: box plots with p-values of  4 subgroups
Subgroups$NEK2_CD274 <- factor(Subgroups$NEK2_CD274 , levels=c("Hi_Hi","Hi_Lo", "Lo_Hi", "Lo_Lo")) # Reorder Order!!!
pwc <- pwc %>% add_xy_position(x =  "NEK2_CD274")
s <- ggboxplot(Subgroups, x = "NEK2_CD274", y = "MonthsOS", ylim = c(0,400), fill = "NEK2_CD274", notch = TRUE,
               palette = c("black","red","brown4","blue")) +
  stat_pvalue_manual(pwc, label = "p", size = 5, tip.length = 0.01, step.increase = 0.08) + #step.increase adjusts distance between signif bars
  labs(
    subtitle = get_test_label(Subgroups.ANOVA, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )
s +
  font("title", size = 15, color = "red", face = "bold.italic")+
  font("subtitle", size = 15, color = "orange")+
  font("caption", size = 15, color = "orange")+
  font("xlab", size = 15, color = "blue")+
  font("ylab", size = 17, color = "#993333")+
  font("xy.text", size = 15, color = "gray2", face = "bold")

###############################################################################################
#### Subsection 6A: Hazard Ratio for NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi (OS)
###############################################################################################
NEK2Hi_CD274Hi_OS <- read_excel("NEK2_CD274_UAMS_OS_NEK2Hi_CD274Hi.xlsx")
NEK2Hi_CD274Hi_OS
NEK2Hi_CD274Hi_OS <- as.data.frame(NEK2Hi_CD274Hi_OS)

## To apply the univariate coxph function to multiple covariates at once, type this:
covariates_OS_Hi_Hi <- c("NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi")
univ_formulas_OS_Hi_Hi <- sapply(covariates_OS_Hi_Hi,
                                 function(x) as.formula(paste('Surv(MonthsOS, CensOS)~',x)))

univ_models_OS_Hi_Hi <- lapply(univ_formulas_OS_Hi_Hi, function(x){coxph(x, data = NEK2Hi_CD274Hi_OS)})
#Extract data
univ_results_OS_Hi_Hi <- lapply(univ_models_OS_Hi_Hi,
                                function(x){
                                  x <- summary(x)
                                  p.value <- signif(x$wald["pvalue"], digits = 3)
                                  wald.test <- signif(x$wald["test"],digits = 2)
                                  beta <- signif(x$coef[1],digits = 3); #coefficient beta
                                  HR <- signif(x$coef[2],digits = 3); #exp(beta)
                                  HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                                  HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                                  HR <- paste0(HR, " (",
                                               HR.confint.lower, "-",HR.confint.upper, ")")
                                  res_OS_Hi_Hi <- c(beta, HR, wald.test, p.value)
                                  names(res_OS_Hi_Hi) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                                  return(res_OS_Hi_Hi)
                                  #return(exp(cbind(coef(x), confint(x))))
                                })
res_OS_Hi_Hi <- t(as.data.frame(univ_results_OS_Hi_Hi, check.names = FALSE))
as.data.frame(res_OS_Hi_Hi)

# Multivariate Cox regression analysis
res.cox_OS_NEK2Hi_CD274Hi <- coxph(Surv(MonthsOS, CensOS) ~ 
                                     NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi, data = NEK2Hi_CD274Hi_OS)
summary(res.cox_OS_NEK2Hi_CD274Hi)

# Visualizing the estimated distribution of survival times
# Plot the baseline survival function
ggsurvplot(survfit(res.cox_OS_NEK2Hi_CD274Hi, data = NEK2Hi_CD274Hi_OS), palette = "#2E9FDF",
           ggtheme = theme_minimal())
# Create the new data  
NEK2Hi_CD274Hi_OS_df <- with(res.cox_OS_NEK2Hi_CD274Hi,
                             data.frame(
                               NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi = c(1,0)
                             ))
NEK2Hi_CD274Hi_OS_df

# Survival curves
NEK2Hi_CD274Hi_OS_fit <- survfit(res.cox_OS_NEK2Hi_CD274Hi, newdata = NEK2Hi_CD274Hi_OS_df)
ggsurvplot(NEK2Hi_CD274Hi_OS_fit, data = NEK2Hi_CD274Hi_OS_df, conf.int = TRUE, legend.labs = c("NEK2Hi_CD274Hi_OS_df = 0", "NEK2Hi_CD274Hi_OS_df = 1"),
           ggtheme = theme_minimal())
print(NEK2Hi_CD274Hi_OS_fit)

###############################################################################################
#### Subsection 6B: Hazard Ratio for NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi (OS)
###############################################################################################
NEK2Hi_CD274Lo_OS <- read_excel("NEK2_CD274_UAMS_OS_NEK2Hi_CD274Lo.xlsx")
NEK2Hi_CD274Lo_OS
NEK2Hi_CD274Lo_OS <- as.data.frame(NEK2Hi_CD274Lo_OS)

## To apply the univariate coxph function to multiple covariates at once, type this:
covariates_OS_Hi_Lo <- c("NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi")
univ_formulas_OS_Hi_Lo <- sapply(covariates_OS_Hi_Lo,
                                 function(x) as.formula(paste('Surv(MonthsOS, CensOS)~',x)))

univ_models_OS_Hi_Lo <- lapply(univ_formulas_OS_Hi_Lo, function(x){coxph(x, data = NEK2Hi_CD274Lo_OS)})
#Extract data
univ_results_OS_Hi_Lo <- lapply(univ_models_OS_Hi_Lo,
                                function(x){
                                  x <- summary(x)
                                  p.value <- signif(x$wald["pvalue"], digits = 3)
                                  wald.test <- signif(x$wald["test"],digits = 2)
                                  beta <- signif(x$coef[1],digits = 3); #coefficient beta
                                  HR <- signif(x$coef[2],digits = 3); #exp(beta)
                                  HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                                  HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                                  HR <- paste0(HR, " (",
                                               HR.confint.lower, "-",HR.confint.upper, ")")
                                  res_Hi_Lo <- c(beta, HR, wald.test, p.value)
                                  names(res_Hi_Lo) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                                  return(res_Hi_Lo)
                                  #return(exp(cbind(coef(x), confint(x))))
                                })
res_OS_Hi_Lo <- t(as.data.frame(univ_results_OS_Hi_Lo, check.names = FALSE))
as.data.frame(res_OS_Hi_Lo)

# Multivariate Cox regression analysis
res.cox_OS_NEK2Hi_CD274Lo <- coxph(Surv(MonthsOS, CensOS) ~ 
                                     NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi, data = NEK2Hi_CD274Lo_OS)
summary(res.cox_OS_NEK2Hi_CD274Lo)

# Visualizing the estimated distribution of survival times
# Plot the baseline survival function
ggsurvplot(survfit(res.cox_OS_NEK2Hi_CD274Lo, data = NEK2Hi_CD274Lo_OS), palette = "#2E9FDF",
           ggtheme = theme_minimal())
# Create the new data  
NEK2Hi_CD274Lo_OS_df <- with(res.cox_OS_NEK2Hi_CD274Lo,
                             data.frame(
                               NEK2 = c(0,0),
                               CD274 = c(0,0),
                               NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi = c(1,0)
                             ))
NEK2Hi_CD274Lo_OS_df

# Survival curves
NEK2Hi_CD274Lo_OS_fit <- survfit(res.cox_OS_NEK2Hi_CD274Lo, newdata = NEK2Hi_CD274Lo_OS_df)
ggsurvplot(NEK2Hi_CD274Lo_OS_fit, data = NEK2Hi_CD274Lo_OS_df, conf.int = TRUE, legend.labs = c("NEK2Hi_CD274Lo_OS_df = 0", "NEK2Hi_CD274Lo_OS_df = 1"),
           ggtheme = theme_minimal())
print(NEK2Hi_CD274Lo_OS_fit)

###############################################################################################
#### Subsection 6C: Hazard Ratio for NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi (OS)
###############################################################################################
NEK2Lo_CD274Lo_OS <- read_excel("NEK2_CD274_UAMS_OS_NEK2Lo_CD274Lo.xlsx")
NEK2Lo_CD274Lo_OS
NEK2Lo_CD274Lo_OS <- as.data.frame(NEK2Lo_CD274Lo_OS)

## To apply the univariate coxph function to multiple covariates at once, type this:
covariates_OS_Lo_Lo <- c("NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi")
univ_formulas_OS_Lo_Lo <- sapply(covariates_OS_Lo_Lo,
                                 function(x) as.formula(paste('Surv(MonthsOS, CensOS)~',x)))

univ_models_OS_Lo_Lo <- lapply(univ_formulas_OS_Lo_Lo, function(x){coxph(x, data = NEK2Lo_CD274Lo_OS)})
#Extract data
univ_results_OS_Lo_Lo <- lapply(univ_models_OS_Lo_Lo,
                                function(x){
                                  x <- summary(x)
                                  p.value <- signif(x$wald["pvalue"], digits = 3)
                                  wald.test <- signif(x$wald["test"],digits = 2)
                                  beta <- signif(x$coef[1],digits = 3); #coefficient beta
                                  HR <- signif(x$coef[2],digits = 3); #exp(beta)
                                  HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                                  HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                                  HR <- paste0(HR, " (",
                                               HR.confint.lower, "-",HR.confint.upper, ")")
                                  res_OS_Lo_Lo <- c(beta, HR, wald.test, p.value)
                                  names(res_OS_Lo_Lo) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                                  return(res_OS_Lo_Lo)
                                  #return(exp(cbind(coef(x), confint(x))))
                                })
res_OS_Lo_Lo <- t(as.data.frame(univ_results_OS_Lo_Lo, check.names = FALSE))
as.data.frame(res_OS_Lo_Lo)

# Multivariate Cox regression analysis
res.cox_OS_NEK2Lo_CD274Lo <- coxph(Surv(MonthsOS, CensOS) ~ 
                                     NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi, data = NEK2Lo_CD274Lo_OS)
summary(res.cox_OS_NEK2Lo_CD274Lo)

# Visualizing the estimated distribution of survival times
# Plot the baseline survival function
ggsurvplot(survfit(res.cox_OS_NEK2Lo_CD274Lo, data = NEK2Lo_CD274Lo_OS), palette = "#2E9FDF",
           ggtheme = theme_minimal())
# Create the new data  
NEK2Lo_CD274Lo_OS_df <- with(res.cox_OS_NEK2Lo_CD274Lo,
                             data.frame(
                               NEK2 = c(0,0),
                               CD274 = c(0,0),
                               NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi = c(1,0)
                             ))
NEK2Lo_CD274Lo_OS_df

# Survival curves
NEK2Lo_CD274Lo_OS_fit <- survfit(res.cox_OS_NEK2Lo_CD274Lo, newdata = NEK2Lo_CD274Lo_OS_df)
ggsurvplot(NEK2Lo_CD274Lo_OS_fit, data = NEK2Lo_CD274Lo_OS_df, conf.int = TRUE, legend.labs = c("NEK2Lo_CD274Lo_OS_df = 0", "NEK2Lo_CD274Lo_OS_df = 1"),
           ggtheme = theme_minimal())
print(NEK2Lo_CD274Lo_OS_fit)

##################################################################################
##################################################################################
## PROJECT 2: MMRF Survival Data (n592)
##################################################################################
##################################################################################

##################################################################################
### Section 1: Load Packages
##################################################################################
# if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("Biobase")
# BiocManager::install("GEOquery")
# BiocManager::install("Biobase")
# install.packages(c("survival", "survminer"))
# install.packages("tidyverse")
# install.packages("readxl")
# install.packages("bquote")
library(survival)
library(survminer)
library(ranger)
library(ggplot2)
library(ggfortify)
library(scales)
library(FSA)
library(rcompanion)
library(rstatix)
library(ggpubr)
library(tidyverse)
library(GEOquery)
library(gtsummary)
library(dplyr)
library(Biobase)
library(tibble)
library(readxl)
###################################################################################
### SECTION 2:  Import MMRF data and create Expression Set 
###################################################################################
x <- read.table("Assay_NEK2_CD274_MMRF.txt", sep="\t", header = T, row.names = 1)
class(x) # Data frame
x <- as.matrix(x) 
class(x) # Matrix array

# Feature Data
f <- read.table("Feature_NEK2_CD274_MMRF.txt", sep="\t", header =T, row.names = 1)
dim(f)
class(f)

# Phenotype data
p <- read.table("Phenotype_NEK2_CD274_MMRF.txt", sep="\t", header =T, row.names = 1)
dim(p)
class(p)

# Create Expression Set
eset <- ExpressionSet(assayData = x,
                      phenoData = AnnotatedDataFrame(p),
                      featureData = AnnotatedDataFrame(f))

#View the number of features/Proteins (rows) and samples (columns)
dim(eset)
# Access data from an ExpressionSet object
x <- exprs(eset) # you can retreive the expression matrix with the function `exprs`
f <- fData(eset) # Retrieve the feature data with `fData'
p <- pData(eset)# Retrieve the phenotype data with`pData`

# extract information of interest from the phenotype data (pdata)
idx <- which(colnames(pData(eset)) %in%
               c('CensOS','CensEFS', 'MonthsOS', 'MonthsEFS'))

metadata <- data.frame(pData(eset)[,idx],
                       row.names = rownames(pData(eset)))

# check that sample names match exactly between pdata and Z-scores 
all((colnames(x) == rownames(metadata)) == TRUE)

# create a merged pdata and Z-scores object
coxdata <- data.frame(metadata, t(x))

coxdata_df <- tibble::rownames_to_column(coxdata, "CHIPID")

# prepare phenotypes
coxdata$CensOS <- as.numeric(coxdata$CensOS) # Censor for the data
coxdata$MonthsOS <- as.numeric(gsub('^KJX|^KJ', '', coxdata$MonthsOS)) # Similar to DaysOS
coxdata$CensEFS <- as.numeric(coxdata$CensEFS) # Censor for the data
coxdata$MonthsEFS <- as.numeric(gsub('^KJX|^KJ', '', coxdata$MonthsEFS)) # Similar to DaysOS

############################################################################################
#### Subsection A: Finding Event Free Survival (EFS) Cut-point for NEK2 and CD274 
############################################################################################
# surv_cutpoint(): Determine the optimal cut-point for each variable using 'maxstat'.
res.cut_EFS <-surv_cutpoint(
  coxdata,
  time = "MonthsEFS",
  event = "CensEFS",
  variables = colnames(coxdata)[5:ncol(coxdata)],
  minprop = 0.10,
  progressbar = TRUE
)
options(max.print = 5000)        # Change global options
summary(res.cut_EFS)

# RES.CAT (categorize by High and Low expression from cutpoints)
res.cat_EFS <- surv_categorize(res.cut_EFS) ### FINE
head(res.cat_EFS)
summary(res.cat_EFS)
str(res.cat_EFS)
#################################################################################
##### Sub-subsection I: Kaplan-Meier(EFS) for NEK2 (ENSG00000117650)
#################################################################################
# 4. Fit survival curves and visualize
{fit_NEK2_EFS <- survfit(Surv(MonthsEFS, CensEFS) ~NEK2_ENSG00000117650, data = res.cat_EFS)
a <- surv_pvalue(fit_NEK2_EFS)$pval
b <- unname(summary(fit_NEK2_EFS)$table[,'records'])
print(a)
print(b)
print(fit_NEK2_EFS)
ggsurvplot(fit_NEK2_EFS, data = res.cat_EFS, risk.table = TRUE, conf.int = F, pval = TRUE)
}

ggsurv_NEK2_EFS <- ggsurvplot(fit_NEK2_EFS, data = res.cat_EFS, 
                              break.time.by = 25,
                              conf.int=TRUE, pval=TRUE, risk.table = TRUE, 
                              surv.median.line = "hv", # Specify median survival
                              tables.theme = theme_survminer(
                                font.main = c(20, "bold", "darkblue"),
                                font.submain = c(15, "bold.italic", "purple"),
                                font.caption = c(14, "plain", "oCBX3ge"),
                                font.x = c(20, "bold.italic", "red"),
                                font.y = c(15, "bold.italic", "darkred"),
                                font.tickslab = c(16, "plain", "darkgreen")),
                              legend.labs=c("High Expression","Low Expression"), legend.title="NEK2",  
                              palette=c("red3","blue"), xlab = "Time (Months)", 
                              title="          NEK2 (ENSG00000117650) EFS (MMRF)",
                              caption = "created with survminer",
                              font.title = c(30, "bold", "darkblue"),
                              font.subtitle = c(18, "bold.italic", "purple"),
                              font.legend = c(18),
                              font.x = c(18, "bold.italic", "red"),
                              font.y = c(18, "bold.italic", "darkred"),
                              font.tickslab = c(16, "plain", "darkgreen"),
                              risk.table.fontsize = 6,
                              ########## risk table #########,
                              risk.table.height = 0.25)

ggsurv_NEK2_EFS
#################################################################################
##### Sub-subsection II: Kaplan-Meier(EFS) for CD274 (ENSG00000120217)
#################################################################################
# Fit survival curves and visualize
{fit_CD274_EFS <- survfit(Surv(MonthsEFS, CensEFS) ~CD274_ENSG00000120217, data = res.cat_EFS)
a <- surv_pvalue(fit_CD274_EFS)$pval
b <- unname(summary(fit_CD274_EFS)$table[,'records'])
print(a)
print(b)
print(fit_CD274_EFS)
ggsurvplot(fit_CD274_EFS, data = res.cat_EFS, risk.table = TRUE, conf.int = F, pval = TRUE)
}

ggsurv_CD274_EFS <- ggsurvplot(fit_CD274_EFS, data = res.cat_EFS, 
                               break.time.by = 25,
                               conf.int=TRUE, pval=TRUE, risk.table = TRUE, 
                               surv.median.line = "hv", # Specify median survival
                               tables.theme = theme_survminer(
                                 font.main = c(20, "bold", "darkblue"),
                                 font.submain = c(15, "bold.italic", "purple"),
                                 font.caption = c(14, "plain", "oCBX3ge"),
                                 font.x = c(20, "bold.italic", "red"),
                                 font.y = c(15, "bold.italic", "darkred"),
                                 font.tickslab = c(16, "plain", "darkgreen")),
                               legend.labs=c("High Expression","Low Expression"), legend.title="CD274",  
                               palette=c("red3","blue"), xlab = "Time (Months)", 
                               title="          CD274 (ENSG00000120217) EFS (MMRF)",
                               caption = "created with survminer",
                               font.title = c(30, "bold", "darkblue"),
                               font.subtitle = c(18, "bold.italic", "purple"),
                               font.legend = c(18),
                               font.x = c(18, "bold.italic", "red"),
                               font.y = c(18, "bold.italic", "darkred"),
                               font.tickslab = c(16, "plain", "darkgreen"),
                               risk.table.fontsize = 6,
                               ########## risk table #########,
                               risk.table.height = 0.25)
ggsurv_CD274_EFS
################################################################################
#### Subsection B: Finding Overall Survival (OS) Cut-point for NEK2 and CD274
################################################################################
# surv_cutpoint(): Determine the optimal cut-point for each variable using 'maxstat'.
res.cut_OS <-surv_cutpoint(
  coxdata,
  time = "MonthsOS",
  event = "CensOS",
  variables = colnames(coxdata)[5:ncol(coxdata)],
  minprop = 0.1,
  progressbar = TRUE
)
options(max.print = 5000)        # Change global options
summary(res.cut_OS)

# RES.CAT (categorize by High and Low expression from cutpoints)
res.cat_OS <- surv_categorize(res.cut_OS) ### FINE
head(res.cat_OS)
summary(res.cat_OS)
str(res.cat_OS)
#################################################################################
##### Sub-subsection I:  Kaplan-Meier(OS) for NEK2 (ENSG00000117650)
#################################################################################
# Fit survival curves and visualize
{fit_NEK2_OS <- survfit(Surv(MonthsOS, CensOS) ~NEK2_ENSG00000117650, data = res.cat_OS)
a <- surv_pvalue(fit_NEK2_OS)$pval
b <- unname(summary(fit_NEK2_OS)$table[,'records'])
print(a)
print(b)
print(fit_NEK2_OS)
ggsurvplot(fit_NEK2_OS, data = res.cat_OS, risk.table = TRUE, conf.int = F, pval = TRUE)
}

ggsurv_NEK2_OS <- ggsurvplot(fit_NEK2_OS, data = res.cat_OS, 
                             break.time.by = 25,
                             conf.int=TRUE, pval=TRUE, risk.table = TRUE, 
                             surv.median.line = "hv", # Specify median survival
                             tables.theme = theme_survminer(
                               font.main = c(20, "bold", "darkblue"),
                               font.submain = c(15, "bold.italic", "purple"),
                               font.caption = c(14, "plain", "oCBX3ge"),
                               font.x = c(20, "bold.italic", "red"),
                               font.y = c(15, "bold.italic", "darkred"),
                               font.tickslab = c(16, "plain", "darkgreen")),
                             legend.labs=c("High Expression","Low Expression"), legend.title="NEK2",  
                             palette=c("red3","blue"), xlab = "Time (Months)", 
                             title="          NEK2 (ENSG00000117650) EFS (MMRF)",
                             caption = "created with survminer",
                             font.title = c(30, "bold", "darkblue"),
                             font.subtitle = c(18, "bold.italic", "purple"),
                             font.legend = c(18),
                             font.x = c(18, "bold.italic", "red"),
                             font.y = c(18, "bold.italic", "darkred"),
                             font.tickslab = c(16, "plain", "darkgreen"),
                             risk.table.fontsize = 6,
                             ########## risk table #########,
                             risk.table.height = 0.25)
ggsurv_NEK2_OS
#################################################################################
##### Sub-subsection II: Kaplan-Meier(OS) for CD274 (ENSG00000120217)
#################################################################################
# Fit survival curves and visualize
{fit_CD274_OS <- survfit(Surv(MonthsOS, CensOS) ~CD274_ENSG00000120217, data = res.cat_OS)
a <- surv_pvalue(fit_CD274_OS)$pval
b <- unname(summary(fit_CD274_OS)$table[,'records'])
print(a)
print(b)
print(fit_CD274_OS)
ggsurvplot(fit_CD274_OS, data = res.cat_OS, risk.table = TRUE, conf.int = F, pval = TRUE)
}

ggsurv_CD274_OS <- ggsurvplot(fit_CD274_OS, data = res.cat_OS, 
                              break.time.by = 25,
                              conf.int=TRUE, pval=TRUE, risk.table = TRUE, 
                              surv.median.line = "hv", # Specify median survival
                              tables.theme = theme_survminer(
                                font.main = c(20, "bold", "darkblue"),
                                font.submain = c(15, "bold.italic", "purple"),
                                font.caption = c(14, "plain", "oCBX3ge"),
                                font.x = c(20, "bold.italic", "red"),
                                font.y = c(15, "bold.italic", "darkred"),
                                font.tickslab = c(16, "plain", "darkgreen")),
                              legend.labs=c("High Expression","Low Expression"), legend.title="CD274",  
                              palette=c("red3","blue"), xlab = "Time (Months)", 
                              title="          CD274 (ENSG00000120217) EFS (MMRF",
                              caption = "created with survminer",
                              font.title = c(30, "bold", "darkblue"),
                              font.subtitle = c(18, "bold.italic", "purple"),
                              font.legend = c(18),
                              font.x = c(18, "bold.italic", "red"),
                              font.y = c(18, "bold.italic", "darkred"),
                              font.tickslab = c(16, "plain", "darkgreen"),
                              risk.table.fontsize = 6,
                              ########## risk table #########,
                              risk.table.height = 0.25)
ggsurv_CD274_OS
################################################################################################
### SECTION 3 - Hazard Ratios for NEK2 and CD274 (EFS)
################################################################################################
# Import EFS MMRF Data
NEK2_CD274_EFS <- read_excel("NEK2_CD274_MMRF_EFS.xlsx")
NEK2_CD274_EFS
NEK2_CD274_EFS <- as.data.frame(NEK2_CD274_EFS)

## To apply the univariate coxph function to multiple covariates at once for EFS, type this:
covariates_EFS <- c("NEK2", "CD274")
univ_formulas_EFS <- sapply(covariates_EFS,
                            function(x) as.formula(paste('Surv(MonthsEFS, CensEFS)~',x)))

univ_models_EFS <- lapply(univ_formulas_EFS, function(x){coxph(x, data = NEK2_CD274_EFS)})
#Extract data
univ_results_EFS <- lapply(univ_models_EFS,
                           function(x){
                             x <- summary(x)
                             p.value <- signif(x$wald["pvalue"], digits = 3)
                             wald.test <- signif(x$wald["test"],digits = 2)
                             beta <- signif(x$coef[1],digits = 3); #coefficient beta
                             HR <- signif(x$coef[2],digits = 3); #exp(beta)
                             HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                             HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                             HR <- paste0(HR, " (",
                                          HR.confint.lower, "-",HR.confint.upper, ")")
                             res <- c(beta, HR, wald.test, p.value)
                             names(res) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                             return(res)
                             #return(exp(cbind(coef(x), confint(x))))
                           })
res_EFS <- t(as.data.frame(univ_results_EFS, check.names = FALSE))
as.data.frame(res_EFS)

# Multivariate Cox regression analysis
res_cox_EFS <- coxph(Surv(MonthsEFS, CensEFS) ~ 
                       NEK2 + CD274, data = NEK2_CD274_EFS)
summary(res_cox_EFS)
################################################################################################
### SECTION 4 - Hazard Ratios for NEK2 and CD274 (OS)
################################################################################################
# Import OS MMRF Data
NEK2_CD274_OS <- read_excel("NEK2_CD274_MMRF_OS.xlsx")
NEK2_CD274_OS
NEK2_CD274_OS <- as.data.frame(NEK2_CD274_OS)

## To apply the univariate coxph function to multiple covariates at once for OS, type this:
covariates_OS <- c("NEK2", "CD274")
univ_formulas_OS <- sapply(covariates_OS,
                           function(x) as.formula(paste('Surv(MonthsOS, CensOS)~',x)))

univ_models_OS <- lapply(univ_formulas_OS, function(x){coxph(x, data = NEK2_CD274_OS)})
#Extract data
univ_results_OS <- lapply(univ_models_OS,
                          function(x){
                            x <- summary(x)
                            p.value <- signif(x$wald["pvalue"], digits = 3)
                            wald.test <- signif(x$wald["test"],digits = 2)
                            beta <- signif(x$coef[1],digits = 3); #coefficient beta
                            HR <- signif(x$coef[2],digits = 3); #exp(beta)
                            HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                            HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                            HR <- paste0(HR, " (",
                                         HR.confint.lower, "-",HR.confint.upper, ")")
                            res <- c(beta, HR, wald.test, p.value)
                            names(res) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                            return(res)
                            #return(exp(cbind(coef(x), confint(x))))
                          })
res_OS <- t(as.data.frame(univ_results_OS, check.names = FALSE))
as.data.frame(res_OS)

# Multivariate Cox regression analysis
res_cox_OS <- coxph(Surv(MonthsOS, CensOS) ~ 
                      NEK2 + CD274, data = NEK2_CD274_OS)
summary(res_cox_OS)
############################################################################################
### SECTION 5 - Import and prepare NEK2_CD274 SUBGROUP combinations for EFS analysis
#############################################################################################
# Import and prepare data
Subgroups <- read_excel("NEK2_CD274_MMRF_EFS.xlsx")
Subgroups
class(Subgroups)
Subgroups <- as.data.frame(Subgroups)
Subgroups.cox <- coxph(Surv(MonthsEFS, CensEFS) ~ Lo_Hi + Hi_Lo + Hi_Hi, data = Subgroups)
summary(Subgroups.cox)

# Median OS and 95% CI for 4 subgroups
Summarize(MonthsEFS ~ NEK2_CD274, #This needs to be name of column in excel
          data=Subgroups,
          digits=3)
groupwiseMedian(MonthsEFS ~ NEK2_CD274,
                data       = Subgroups,
                conf       = 0.95,
                R          = 5000,
                percentile = TRUE,
                bca        = FALSE,
                digits     = 4)
# Histogram and survival curves for  4 subgroups
class(Subgroups)
Subgroups <- as.data.frame(Subgroups)

# Reorder the Columns first
Risk_order <- c("Lo_Lo", "Lo_Hi", "Hi_Lo", "Hi_Hi") #REORDER the columns!
ggplot(Subgroups, aes(x = factor(NEK2_CD274, Risk_order))) +  # Put Risk_order in to help
  geom_bar(aes(y = (..count..)/sum(..count..))) +
  scale_y_continuous(labels=percent) +
  labs(title = 'Distribution (Using MM Date)', x = 'Score Range', y = 'Percentage %') +
  theme_minimal() +
  theme(axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        text = element_text(size=25))
#Next, we look at survival curves by subgroups.
km_trt_fit_EFS <- survfit(Surv(MonthsEFS, CensEFS) ~ NEK2_CD274, data = Subgroups)
p <- ggsurvplot(km_trt_fit_EFS, conf.int = FALSE,
                break.time.by = 40,
                risk.table = "nrisk_cumevents",
                surv.median.line = "hv", # Specify median survival
                tables.theme = theme_survminer(
                  font.main = c(20, "bold", "black"),
                  font.submain = c(20, "bold.italic", "black"),
                  font.caption = c(10, "plain", "black"),
                  font.x = c(20, "bold.italic", "black"),
                  font.y = c(15, "bold.italic", "black"),
                  font.tickslab = c(12, "plain", "black")),
                legend.labs=c("NEK2-high/CD274-high", "NEK2-high/CD274-low","NEK2-low/CD274-high","NEK2-low/CD274-low"), legend.title="",  
                palette=c("black","red","brown4","blue"), xlab = "Time (Months)", 
                title="                      NEK2/CD274 EFS (MMRF)",
                caption = "created with survminer",
                font.title = c(25, "bold", "black"),
                font.subtitle = c(16, "bold.italic", "black"),
                font.legend = c(12),
                font.x = c(18, "bold.italic", "black"),
                font.y = c(18, "bold.italic", "black"),
                font.tickslab = c(16, "plain", "black"),
                risk.table.fontsize = 4.5,
                ########## risk table #########,
                risk.table.height = 0.25)
p

# Compare the means of the 4 subgroups using ANOVA test
Subgroups.ANOVA <- Subgroups %>% anova_test(MonthsEFS ~ NEK2_CD274)
Subgroups.ANOVA
# Pairwise T-tests for  4 subgroups
pwc <- Subgroups %>%
  # pairwise_t_test(MonthsEFS ~ NEK2_CD274_EFS, p.adjust.method = "bonferroni")
  pairwise_t_test(MonthsEFS ~ NEK2_CD274)
pwc
# Visualization: box plots with p-values of  4 subgroups
Subgroups$NEK2_CD274 <- factor(Subgroups$NEK2_CD274 , levels=c("Hi_Hi","Hi_Lo", "Lo_Hi", "Lo_Lo")) # Reorder Order!!!
pwc <- pwc %>% add_xy_position(x =  "NEK2_CD274")
r <- ggboxplot(Subgroups, x = "NEK2_CD274", y = "MonthsEFS", ylim = c(0,401), fill = "NEK2_CD274", notch = TRUE,
               palette = c("black","red","brown4","blue")) +
  stat_pvalue_manual(pwc, label = "p", size = 5, tip.length = 0.01, step.increase = 0.08) + #step.increase adjusts distance between signif bars
  labs(
    subtitle = get_test_label(Subgroups.ANOVA, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )
r +
  font("title", size = 15, color = "red", face = "bold.italic")+
  font("subtitle", size = 15, color = "orange")+
  font("caption", size = 15, color = "orange")+
  font("xlab", size = 15, color = "blue")+
  font("ylab", size = 17, color = "#993333")+
  font("xy.text", size = 15, color = "gray2", face = "bold")

###############################################################################################
#### Subsection 5A: Hazard Ratio for NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi (EFS)
###############################################################################################

NEK2Hi_CD274Hi_EFS <- read_excel("MMRF_EFS_NEK2Hi_CD274Hi.xlsx")
NEK2Hi_CD274Hi_EFS
NEK2Hi_CD274Hi_EFS <- as.data.frame(NEK2Hi_CD274Hi_EFS)

## To apply the univariate coxph function to multiple covariates at once, type this:
covariates_EFS_Hi_Hi <- c("NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi")
univ_formulas_EFS_Hi_Hi <- sapply(covariates_EFS_Hi_Hi,
                                  function(x) as.formula(paste('Surv(MonthsEFS, CensEFS)~',x)))

univ_models_EFS_Hi_Hi <- lapply(univ_formulas_EFS_Hi_Hi, function(x){coxph(x, data = NEK2Hi_CD274Hi_EFS)})
#Extract data
univ_results_EFS_Hi_Hi <- lapply(univ_models_EFS_Hi_Hi,
                                 function(x){
                                   x <- summary(x)
                                   p.value <- signif(x$wald["pvalue"], digits = 3)
                                   wald.test <- signif(x$wald["test"],digits = 2)
                                   beta <- signif(x$coef[1],digits = 3); #coefficient beta
                                   HR <- signif(x$coef[2],digits = 3); #exp(beta)
                                   HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                                   HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                                   HR <- paste0(HR, " (",
                                                HR.confint.lower, "-",HR.confint.upper, ")")
                                   res_EFS_Hi_Hi <- c(beta, HR, wald.test, p.value)
                                   names(res_EFS_Hi_Hi) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                                   return(res_EFS_Hi_Hi)
                                   #return(exp(cbind(coef(x), confint(x))))
                                 })
res_EFS_Hi_Hi <- t(as.data.frame(univ_results_EFS_Hi_Hi, check.names = FALSE))
as.data.frame(res_EFS_Hi_Hi)


# Multivariate Cox regression analysis
res.cox_EFS_NEK2Hi_CD274Hi <- coxph(Surv(MonthsEFS, CensEFS) ~ 
                                      NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi, data = NEK2Hi_CD274Hi_EFS)
summary(res.cox_EFS_NEK2Hi_CD274Hi)

# Visualizing the estimated distribution of survival times
# Plot the baseline survival function
ggsurvplot(survfit(res.cox_EFS_NEK2Hi_CD274Hi, data = NEK2Hi_CD274Hi_EFS), palette = "#2E9FDF",
           ggtheme = theme_minimal())
# Create the new data  
NEK2Hi_CD274Hi_EFS_df <- with(res.cox_EFS_NEK2Hi_CD274Hi,
                              data.frame(
                                NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi = c(1,0)
                              ))
NEK2Hi_CD274Hi_EFS_df

# Survival curves
NEK2Hi_CD274Hi_EFS_fit <- survfit(res.cox_EFS_NEK2Hi_CD274Hi, newdata = NEK2Hi_CD274Hi_EFS_df)
ggsurvplot(NEK2Hi_CD274Hi_EFS_fit, data = NEK2Hi_CD274Hi_EFS_df, conf.int = TRUE, legend.labs = c("NEK2Hi_CD274Hi_EFS_df = 0", "NEK2Hi_CD274Hi_EFS_df = 1"),
           ggtheme = theme_minimal())
print(NEK2Hi_CD274Hi_EFS_fit)

###############################################################################################
#### Subsection 5B: Hazard Ratio for NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi (EFS)
###############################################################################################

NEK2Hi_CD274Lo_EFS <- read_excel("MMRF_EFS_NEK2Hi_CD274Lo.xlsx")
NEK2Hi_CD274Lo_EFS
NEK2Hi_CD274Lo_EFS <- as.data.frame(NEK2Hi_CD274Lo_EFS)

## To apply the univariate coxph function to multiple covariates at once, type this:
covariates_EFS_Hi_Lo <- c("NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi")
univ_formulas_EFS_Hi_Lo <- sapply(covariates_EFS_Hi_Lo,
                                  function(x) as.formula(paste('Surv(MonthsEFS, CensEFS)~',x)))

univ_models_EFS_Hi_Lo <- lapply(univ_formulas_EFS_Hi_Lo, function(x){coxph(x, data = NEK2Hi_CD274Lo_EFS)})
#Extract data
univ_results_EFS_Hi_Lo <- lapply(univ_models_EFS_Hi_Lo,
                                 function(x){
                                   x <- summary(x)
                                   p.value <- signif(x$wald["pvalue"], digits = 3)
                                   wald.test <- signif(x$wald["test"],digits = 2)
                                   beta <- signif(x$coef[1],digits = 3); #coefficient beta
                                   HR <- signif(x$coef[2],digits = 3); #exp(beta)
                                   HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                                   HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                                   HR <- paste0(HR, " (",
                                                HR.confint.lower, "-",HR.confint.upper, ")")
                                   res_Hi_Lo <- c(beta, HR, wald.test, p.value)
                                   names(res_Hi_Lo) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                                   return(res_Hi_Lo)
                                   #return(exp(cbind(coef(x), confint(x))))
                                 })
res_EFS_Hi_Lo <- t(as.data.frame(univ_results_EFS_Hi_Lo, check.names = FALSE))
as.data.frame(res_EFS_Hi_Lo)

# Multivariate Cox regression analysis
res.cox_EFS_NEK2Hi_CD274Lo <- coxph(Surv(MonthsEFS, CensEFS) ~ 
                                      NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi, data = NEK2Hi_CD274Lo_EFS)
summary(res.cox_EFS_NEK2Hi_CD274Lo)

# Visualizing the estimated distribution of survival times
# Plot the baseline survival function
ggsurvplot(survfit(res.cox_EFS_NEK2Hi_CD274Lo, data = NEK2Hi_CD274Lo_EFS), palette = "#2E9FDF",
           ggtheme = theme_minimal())
# Create the new data  
NEK2Hi_CD274Lo_EFS_df <- with(res.cox_EFS_NEK2Hi_CD274Lo,
                              data.frame(
                                NEK2 = c(0,0),
                                CD274 = c(0,0),
                                NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi = c(1,0)
                              ))
NEK2Hi_CD274Lo_EFS_df

# Survival curves
NEK2Hi_CD274Lo_EFS_fit <- survfit(res.cox_EFS_NEK2Hi_CD274Lo, newdata = NEK2Hi_CD274Lo_EFS_df)
ggsurvplot(NEK2Hi_CD274Lo_EFS_fit, data = NEK2Hi_CD274Lo_EFS_df, conf.int = TRUE, legend.labs = c("NEK2Hi_CD274Lo_EFS_df = 0", "NEK2Hi_CD274Lo_EFS_df = 1"),
           ggtheme = theme_minimal())
print(NEK2Hi_CD274Lo_EFS_fit)

###############################################################################################
#### Subsection 5C: Hazard Ratio for NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi (EFS)
###############################################################################################

NEK2Lo_CD274Lo_EFS <- read_excel("MMRF_EFS_NEK2Lo_CD274Lo.xlsx")
NEK2Lo_CD274Lo_EFS
NEK2Lo_CD274Lo_EFS <- as.data.frame(NEK2Lo_CD274Lo_EFS)

## To apply the univariate coxph function to multiple covariates at once, type this:
covariates_EFS_Lo_Lo <- c("NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi")
univ_formulas_EFS_Lo_Lo <- sapply(covariates_EFS_Lo_Lo,
                                  function(x) as.formula(paste('Surv(MonthsEFS, CensEFS)~',x)))

univ_models_EFS_Lo_Lo <- lapply(univ_formulas_EFS_Lo_Lo, function(x){coxph(x, data = NEK2Lo_CD274Lo_EFS)})
#Extract data
univ_results_EFS_Lo_Lo <- lapply(univ_models_EFS_Lo_Lo,
                                 function(x){
                                   x <- summary(x)
                                   p.value <- signif(x$wald["pvalue"], digits = 3)
                                   wald.test <- signif(x$wald["test"],digits = 2)
                                   beta <- signif(x$coef[1],digits = 3); #coefficient beta
                                   HR <- signif(x$coef[2],digits = 3); #exp(beta)
                                   HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                                   HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                                   HR <- paste0(HR, " (",
                                                HR.confint.lower, "-",HR.confint.upper, ")")
                                   res_EFS_Lo_Lo <- c(beta, HR, wald.test, p.value)
                                   names(res_EFS_Lo_Lo) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                                   return(res_EFS_Lo_Lo)
                                   #return(exp(cbind(coef(x), confint(x))))
                                 })
res_EFS_Lo_Lo <- t(as.data.frame(univ_results_EFS_Lo_Lo, check.names = FALSE))
as.data.frame(res_EFS_Lo_Lo)

# Multivariate Cox regression analysis
res.cox_EFS_NEK2Lo_CD274Lo <- coxph(Surv(MonthsEFS, CensEFS) ~ 
                                      NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi, data = NEK2Lo_CD274Lo_EFS)
summary(res.cox_EFS_NEK2Lo_CD274Lo)

# Visualizing the estimated distribution of survival times
# Plot the baseline survival function
ggsurvplot(survfit(res.cox_EFS_NEK2Lo_CD274Lo, data = NEK2Lo_CD274Lo_EFS), palette = "#2E9FDF",
           ggtheme = theme_minimal())
# Create the new data  
NEK2Lo_CD274Lo_EFS_df <- with(res.cox_EFS_NEK2Lo_CD274Lo,
                              data.frame(
                                NEK2 = c(0,0),
                                CD274 = c(0,0),
                                NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi = c(1,0)
                              ))
NEK2Lo_CD274Lo_EFS_df

# Survival curves
NEK2Lo_CD274Lo_EFS_fit <- survfit(res.cox_EFS_NEK2Lo_CD274Lo, newdata = NEK2Lo_CD274Lo_EFS_df)
ggsurvplot(NEK2Lo_CD274Lo_EFS_fit, data = NEK2Lo_CD274Lo_EFS_df, conf.int = TRUE, legend.labs = c("NEK2Lo_CD274Lo_EFS_df = 0", "NEK2Lo_CD274Lo_EFS_df = 1"),
           ggtheme = theme_minimal())
print(NEK2Lo_CD274Lo_EFS_fit)

###########################################################################################
### SECTION 6 - Import and prepare NEK2_CD274 SUBGROUP combinations for OS analysis
###########################################################################################
# Multivariate Cox regression analysis of  4 subgroups
Subgroups <- read_excel("NEK2_CD274_MMRF_OS.xlsx")
Subgroups
class(Subgroups)
Subgroups <- as.data.frame(Subgroups)
Subgroups.cox <- coxph(Surv(MonthsOS, CensOS) ~ Lo_Hi + Hi_Lo + Hi_Hi, data = Subgroups)
summary(Subgroups.cox)

# Median OS and 95% CI for  4 subgroups
Summarize(MonthsOS ~ NEK2_CD274, #This needs to be name of column in excel
          data=Subgroups,
          digits=3)
groupwiseMedian(MonthsOS ~ NEK2_CD274,
                data       = Subgroups,
                conf       = 0.95,
                R          = 5000,
                percentile = TRUE,
                bca        = FALSE,
                digits     = 4)
# Histogram and survival curves for  4 subgroups
class(Subgroups)
Subgroups <- as.data.frame(Subgroups)

# Reorder the Columns first
Risk_order <- c("Lo_Lo", "Lo_Hi", "Hi_Lo", "Hi_Hi") #REORDER the columns!
ggplot(Subgroups, aes(x = factor(NEK2_CD274, Risk_order))) +  # Put Risk_order in to help
  geom_bar(aes(y = (..count..)/sum(..count..))) +
  scale_y_continuous(labels=percent) +
  labs(title = 'Distribution (Using MM Date)', x = 'Score Range', y = 'Percentage %') +
  theme_minimal() +
  theme(axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        text = element_text(size=25))
#Next, we look at survival curves by subgroups.
km_trt_fit_CD274_OS <- survfit(Surv(MonthsOS, CensOS) ~ NEK2_CD274, data = Subgroups)
q <- ggsurvplot(km_trt_fit_CD274_OS, conf.int = FALSE,
                break.time.by = 40,
                risk.table = "nrisk_cumevents",
                surv.median.line = "hv", # Specify median survival
                tables.theme = theme_survminer(
                  font.main = c(20, "bold", "black"),
                  font.submain = c(20, "bold.italic", "black"),
                  font.caption = c(10, "plain", "black"),
                  font.x = c(20, "bold.italic", "black"),
                  font.y = c(15, "bold.italic", "black"),
                  font.tickslab = c(12, "plain", "black")),
                legend.labs=c("NEK2-high/CD274-high","NEK2-high/CD274-low","NEK2-low/CD274-high","NEK2-low/CD274-low"), legend.title="",  
                palette=c("black","red","brown4","blue"), xlab = "Time (Months)", 
                title="                      NEK2/CD274 OS (MMRF)",
                caption = "created with survminer",
                font.title = c(25, "bold", "black"),
                font.subtitle = c(16, "bold.italic", "black"),
                font.legend = c(12),
                font.x = c(18, "bold.italic", "black"),
                font.y = c(18, "bold.italic", "black"),
                font.tickslab = c(16, "plain", "black"),
                risk.table.fontsize = 4.5,
                ########## risk table #########,
                risk.table.height = 0.25)
q
# Compare the mean of Subgroups using ANOVA test
Subgroups.ANOVA <- Subgroups %>% anova_test(MonthsOS ~ NEK2_CD274)
Subgroups.ANOVA
# Pairwise T-tests for  4 subgroups
pwc <- Subgroups %>%
  # pairwise_t_test(MonthsOS ~ NEK2_CD274_OS, p.adjust.method = "bonferroni")
  pairwise_t_test(MonthsOS ~ NEK2_CD274)
pwc
# Visualization: box plots with p-values of  4 subgroups
Subgroups$NEK2_CD274 <- factor(Subgroups$NEK2_CD274 , levels=c("Hi_Hi","Hi_Lo", "Lo_Hi", "Lo_Lo")) # Reorder Order!!!
pwc <- pwc %>% add_xy_position(x =  "NEK2_CD274")
s <- ggboxplot(Subgroups, x = "NEK2_CD274", y = "MonthsOS", ylim = c(0,400), fill = "NEK2_CD274", notch = TRUE,
               palette = c("black","red","brown4","blue")) +
  stat_pvalue_manual(pwc, label = "p", size = 5, tip.length = 0.01, step.increase = 0.08) + #step.increase adjusts distance between signif bars
  labs(
    subtitle = get_test_label(Subgroups.ANOVA, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )
s +
  font("title", size = 15, color = "red", face = "bold.italic")+
  font("subtitle", size = 15, color = "orange")+
  font("caption", size = 15, color = "orange")+
  font("xlab", size = 15, color = "blue")+
  font("ylab", size = 17, color = "#993333")+
  font("xy.text", size = 15, color = "gray2", face = "bold")

###############################################################################################
#### Subsection 6A: Hazard Ratio for NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi (OS)
###############################################################################################
NEK2Hi_CD274Hi_OS <- read_excel("MMRF_OS_NEK2Hi_CD274Hi.xlsx")
NEK2Hi_CD274Hi_OS
NEK2Hi_CD274Hi_OS <- as.data.frame(NEK2Hi_CD274Hi_OS)

## To apply the univariate coxph function to multiple covariates at once, type this:
covariates_OS_Hi_Hi <- c("NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi")
univ_formulas_OS_Hi_Hi <- sapply(covariates_OS_Hi_Hi,
                                 function(x) as.formula(paste('Surv(MonthsOS, CensOS)~',x)))

univ_models_OS_Hi_Hi <- lapply(univ_formulas_OS_Hi_Hi, function(x){coxph(x, data = NEK2Hi_CD274Hi_OS)})
#Extract data
univ_results_OS_Hi_Hi <- lapply(univ_models_OS_Hi_Hi,
                                function(x){
                                  x <- summary(x)
                                  p.value <- signif(x$wald["pvalue"], digits = 3)
                                  wald.test <- signif(x$wald["test"],digits = 2)
                                  beta <- signif(x$coef[1],digits = 3); #coefficient beta
                                  HR <- signif(x$coef[2],digits = 3); #exp(beta)
                                  HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                                  HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                                  HR <- paste0(HR, " (",
                                               HR.confint.lower, "-",HR.confint.upper, ")")
                                  res_OS_Hi_Hi <- c(beta, HR, wald.test, p.value)
                                  names(res_OS_Hi_Hi) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                                  return(res_OS_Hi_Hi)
                                  #return(exp(cbind(coef(x), confint(x))))
                                })
res_OS_Hi_Hi <- t(as.data.frame(univ_results_OS_Hi_Hi, check.names = FALSE))
as.data.frame(res_OS_Hi_Hi)


# Multivariate Cox regression analysis
res.cox_OS_NEK2Hi_CD274Hi <- coxph(Surv(MonthsOS, CensOS) ~ 
                                     NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi, data = NEK2Hi_CD274Hi_OS)
summary(res.cox_OS_NEK2Hi_CD274Hi)

# Visualizing the estimated distribution of survival times
# Plot the baseline survival function
ggsurvplot(survfit(res.cox_OS_NEK2Hi_CD274Hi, data = NEK2Hi_CD274Hi_OS), palette = "#2E9FDF",
           ggtheme = theme_minimal())
# Create the new data  
NEK2Hi_CD274Hi_OS_df <- with(res.cox_OS_NEK2Hi_CD274Hi,
                             data.frame(
                               NEK2Hi_CD274Hi_vs_NEK2Lo_CD274Hi = c(1,0)
                             ))
NEK2Hi_CD274Hi_OS_df

# Survival curves
NEK2Hi_CD274Hi_OS_fit <- survfit(res.cox_OS_NEK2Hi_CD274Hi, newdata = NEK2Hi_CD274Hi_OS_df)
ggsurvplot(NEK2Hi_CD274Hi_OS_fit, data = NEK2Hi_CD274Hi_OS_df, conf.int = TRUE, legend.labs = c("NEK2Hi_CD274Hi_OS_df = 0", "NEK2Hi_CD274Hi_OS_df = 1"),
           ggtheme = theme_minimal())
print(NEK2Hi_CD274Hi_OS_fit)

###############################################################################################
#### Subsection 6B: Hazard Ratio for NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi (OS)
###############################################################################################

NEK2Hi_CD274Lo_OS <- read_excel("MMRF_OS_NEK2Hi_CD274Lo.xlsx")
NEK2Hi_CD274Lo_OS
NEK2Hi_CD274Lo_OS <- as.data.frame(NEK2Hi_CD274Lo_OS)

## To apply the univariate coxph function to multiple covariates at once, type this:
covariates_OS_Hi_Lo <- c("NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi")
univ_formulas_OS_Hi_Lo <- sapply(covariates_OS_Hi_Lo,
                                 function(x) as.formula(paste('Surv(MonthsOS, CensOS)~',x)))

univ_models_OS_Hi_Lo <- lapply(univ_formulas_OS_Hi_Lo, function(x){coxph(x, data = NEK2Hi_CD274Lo_OS)})
#Extract data
univ_results_OS_Hi_Lo <- lapply(univ_models_OS_Hi_Lo,
                                function(x){
                                  x <- summary(x)
                                  p.value <- signif(x$wald["pvalue"], digits = 3)
                                  wald.test <- signif(x$wald["test"],digits = 2)
                                  beta <- signif(x$coef[1],digits = 3); #coefficient beta
                                  HR <- signif(x$coef[2],digits = 3); #exp(beta)
                                  HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                                  HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                                  HR <- paste0(HR, " (",
                                               HR.confint.lower, "-",HR.confint.upper, ")")
                                  res_Hi_Lo <- c(beta, HR, wald.test, p.value)
                                  names(res_Hi_Lo) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                                  return(res_Hi_Lo)
                                  #return(exp(cbind(coef(x), confint(x))))
                                })
res_OS_Hi_Lo <- t(as.data.frame(univ_results_OS_Hi_Lo, check.names = FALSE))
as.data.frame(res_OS_Hi_Lo)

# Multivariate Cox regression analysis
res.cox_OS_NEK2Hi_CD274Lo <- coxph(Surv(MonthsOS, CensOS) ~ 
                                     NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi, data = NEK2Hi_CD274Lo_OS)
summary(res.cox_OS_NEK2Hi_CD274Lo)

# Visualizing the estimated distribution of survival times
# Plot the baseline survival function
ggsurvplot(survfit(res.cox_OS_NEK2Hi_CD274Lo, data = NEK2Hi_CD274Lo_OS), palette = "#2E9FDF",
           ggtheme = theme_minimal())
# Create the new data  
NEK2Hi_CD274Lo_OS_df <- with(res.cox_OS_NEK2Hi_CD274Lo,
                             data.frame(
                               NEK2 = c(0,0),
                               CD274 = c(0,0),
                               NEK2Hi_CD274Lo_vs_NEK2Lo_CD274Hi = c(1,0)
                             ))
NEK2Hi_CD274Lo_OS_df

# Survival curves
NEK2Hi_CD274Lo_OS_fit <- survfit(res.cox_OS_NEK2Hi_CD274Lo, newdata = NEK2Hi_CD274Lo_OS_df)
ggsurvplot(NEK2Hi_CD274Lo_OS_fit, data = NEK2Hi_CD274Lo_OS_df, conf.int = TRUE, legend.labs = c("NEK2Hi_CD274Lo_OS_df = 0", "NEK2Hi_CD274Lo_OS_df = 1"),
           ggtheme = theme_minimal())
print(NEK2Hi_CD274Lo_OS_fit)

###############################################################################################
#### Subsection 6C: Hazard Ratio for NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi (OS)
###############################################################################################
NEK2Lo_CD274Lo_OS <- read_excel("MMRF_OS_NEK2Lo_CD274Lo.xlsx")
NEK2Lo_CD274Lo_OS
NEK2Lo_CD274Lo_OS <- as.data.frame(NEK2Lo_CD274Lo_OS)

## To apply the univariate coxph function to multiple covariates at once, type this:
covariates_OS_Lo_Lo <- c("NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi")
univ_formulas_OS_Lo_Lo <- sapply(covariates_OS_Lo_Lo,
                                 function(x) as.formula(paste('Surv(MonthsOS, CensOS)~',x)))

univ_models_OS_Lo_Lo <- lapply(univ_formulas_OS_Lo_Lo, function(x){coxph(x, data = NEK2Lo_CD274Lo_OS)})
#Extract data
univ_results_OS_Lo_Lo <- lapply(univ_models_OS_Lo_Lo,
                                function(x){
                                  x <- summary(x)
                                  p.value <- signif(x$wald["pvalue"], digits = 3)
                                  wald.test <- signif(x$wald["test"],digits = 2)
                                  beta <- signif(x$coef[1],digits = 3); #coefficient beta
                                  HR <- signif(x$coef[2],digits = 3); #exp(beta)
                                  HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                                  HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                                  HR <- paste0(HR, " (",
                                               HR.confint.lower, "-",HR.confint.upper, ")")
                                  res_OS_Lo_Lo <- c(beta, HR, wald.test, p.value)
                                  names(res_OS_Lo_Lo) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                                  return(res_OS_Lo_Lo)
                                  #return(exp(cbind(coef(x), confint(x))))
                                })
res_OS_Lo_Lo <- t(as.data.frame(univ_results_OS_Lo_Lo, check.names = FALSE))
as.data.frame(res_OS_Lo_Lo)

# Multivariate Cox regression analysis
res.cox_OS_NEK2Lo_CD274Lo <- coxph(Surv(MonthsOS, CensOS) ~ 
                                     NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi, data = NEK2Lo_CD274Lo_OS)
summary(res.cox_OS_NEK2Lo_CD274Lo)

# Visualizing the estimated distribution of survival times
# Plot the baseline survival function
ggsurvplot(survfit(res.cox_OS_NEK2Lo_CD274Lo, data = NEK2Lo_CD274Lo_OS), palette = "#2E9FDF",
           ggtheme = theme_minimal())
# Create the new data  
NEK2Lo_CD274Lo_OS_df <- with(res.cox_OS_NEK2Lo_CD274Lo,
                             data.frame(
                               NEK2 = c(0,0),
                               CD274 = c(0,0),
                               NEK2Lo_CD274Lo_vs_NEK2Lo_CD274Hi = c(1,0)
                             ))
NEK2Lo_CD274Lo_OS_df

# Survival curves
NEK2Lo_CD274Lo_OS_fit <- survfit(res.cox_OS_NEK2Lo_CD274Lo, newdata = NEK2Lo_CD274Lo_OS_df)
ggsurvplot(NEK2Lo_CD274Lo_OS_fit, data = NEK2Lo_CD274Lo_OS_df, conf.int = TRUE, legend.labs = c("NEK2Lo_CD274Lo_OS_df = 0", "NEK2Lo_CD274Lo_OS_df = 1"),
           ggtheme = theme_minimal())
print(NEK2Lo_CD274Lo_OS_fit)