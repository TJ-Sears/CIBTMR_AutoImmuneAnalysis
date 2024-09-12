
# backwards selection of time-dep covariate controlled rel / OS analysis

library(ggplot2)
library(RColorBrewer)
library(survival)
library(ggpubr)
library(tidyverse)
library(dplyr)
#all samples
library(splitstackshape)
library(Rtsne)
library(ggplot2)
library(umap)
library(readxl)
library(dplyr)
library(gridExtra)
library(survminer)
library(ggpubr)
library(ggfortify)
library(survival)

clinical_HLA_phbr<-read.table("/Users/tjsears/Code/CIBMTR/Tables/complete_dataset_Feb25.txt",sep="\t",header=T)
clinical_HLA_phbr$agvhd24[clinical_HLA_phbr$agvhd34!=0]<-0
clinical_HLA_phbr$agvhd2<-clinical_HLA_phbr$agvhd24

test_df_coxPH <- clinical_HLA_phbr %>%
  select(Individual, dead, intxsurv, rel, intxrel,Age.at.diagnosis, condint, agvhd2, agvhd34, ipssrpr, donorgp,
         AI, dnrsex, graftype, sex, cgvhd, dnrage, AI_plus, AI.I, AI.II, Poor..High..PHBR.I, Poor..High..PHBR.II,PHBR.I.change,    
         PHBR.II.change,intxcgvhd,intxagvhd24,intxagvhd34)

test_df_coxPH$dnrage[is.na(test_df_coxPH$dnrage)]<-mean(test_df_coxPH$dnrage,na.rm = T)
test_df_coxPH[is.na(test_df_coxPH)]<-0
test_df_coxPH <- test_df_coxPH %>% 
  mutate(across(everything(), ~ifelse(. == '.', 0, .)))
test_df_coxPH <- test_df_coxPH %>%
  filter(condint != 99, donorgp != 99)

test_df_coxPH$graftype<-as.factor(test_df_coxPH$graftype)
test_df_coxPH$dnrsex<-as.numeric(test_df_coxPH$dnrsex)

test_df_coxPH$ImprovedPHBR.I<--1*test_df_coxPH$Poor..High..PHBR.I
test_df_coxPH$ImprovedPHBR.II<--1*test_df_coxPH$Poor..High..PHBR.II

test_df_coxPH$WorsePHBR.I<-test_df_coxPH$PHBR.I.change
test_df_coxPH$WorsePHBR.II<-test_df_coxPH$PHBR.II.change

# rename confusing variables
test_df_coxPH$Donor_Age<-test_df_coxPH$dnrage
test_df_coxPH$Donor_Age<-test_df_coxPH$dnrage
test_df_coxPH$cGVHD<-test_df_coxPH$cgvhd
test_df_coxPH$graft_type<-test_df_coxPH$graftype22

########
## OS ##
########

time_dep_df<-test_df_coxPH[test_df_coxPH$intxsurv>=1,]

df_time_dep <-
  tibble(Individual = time_dep_df$Individual,
         time = time_dep_df$intxcgvhd,
         cgvhd = time_dep_df$cgvhd)

df_time_dep_24 <-
  tibble(Individual = time_dep_df$Individual,
         time24 = time_dep_df$intxagvhd24,
         agvhd2 = time_dep_df$agvhd2)

df_time_dep_34 <-
  tibble(Individual = time_dep_df$Individual,
         time34 = time_dep_df$intxagvhd34,
         agvhd34 = time_dep_df$agvhd34)

df_ind <- tmerge(data1=time_dep_df[, c("intxsurv","dead","Individual","condint","ipssrpr","donorgp","AI.I","dnrsex","dnrage","WorsePHBR.I","WorsePHBR.II","ImprovedPHBR.I","ImprovedPHBR.II",'AI.II','sex','graftype',"Age.at.diagnosis")], 
                 data2=time_dep_df[, c("intxsurv","dead","Individual","condint","ipssrpr","donorgp","AI.I","dnrsex","dnrage","WorsePHBR.I","WorsePHBR.II","ImprovedPHBR.I","ImprovedPHBR.II",'AI.II','sex','graftype',"Age.at.diagnosis")], 
                 id=Individual, event=event(intxsurv,dead)) 
df_final <-
  tmerge(data1=df_ind,
         data2=df_time_dep,
         id=Individual,
         cgvhd=tdc(time, cgvhd))
df_final <-
  tmerge(data1=df_final,
         data2=df_time_dep_24,
         id=Individual,
         agvhd2=tdc(time24, agvhd2))
df_final <-
  tmerge(data1=df_final,
         data2=df_time_dep_34,
         id=Individual,
         agvhd34=tdc(time34, agvhd34))

df_final$cgvhd[is.na(df_final$cgvhd)]<-0
df_final$agvhd2[is.na(df_final$agvhd2)]<-0
df_final$agvhd34[is.na(df_final$agvhd34)]<-0

# Define function to fit Cox model and find variable with the highest p-value
fit_cox_model <- function(data, variables) {
  print(variables)
  formula <- as.formula(paste("Surv(tstart, tstop, event) ~",paste(variables, collapse = " + ")))
  print(formula)
  model <- coxph(formula, data = data)
  summary_model <- summary(model)
  #print(summary_model)
  worst_var <- rownames(summary_model$coefficients)[which.max(summary_model$coefficients[, "Pr(>|z|)"])]
  aic <- AIC(model)
  return(list(aic = aic, worst_var = worst_var))
}

all_variables <- c("condint","ipssrpr","donorgp","AI.I","dnrsex","dnrage","ImprovedPHBR.I","ImprovedPHBR.II","agvhd34","agvhd2",'AI.II','sex','graftype',"Age.at.diagnosis")

best_aic <- Inf
best_variable_combination <- NULL

# Iterate through all possible combinations to find the best model based on AIC
rerun <- TRUE
if (rerun) {
  repeat {
    res <- fit_cox_model(df_final, c(all_variables))
    aic <- res$aic
    worst_var <- res$worst_var
    if (aic < best_aic) {
      best_aic <- aic
      all_variables <- setdiff(all_variables, worst_var)
      best_variable_combination <- all_variables
      print(aic)
      print(worst_var)
    } else {
      break
    }
  }
}

cat("Best variable combination:", best_variable_combination, "\n")
cat("Best AIC:", best_aic, "\n")

# Fitting the best model
best_formula <- as.formula(paste("Surv(tstart, tstop, event) ~",paste(best_variable_combination, collapse = " + ")))
final_model <- coxph(best_formula, data = df_final)
os_summary<-summary(final_model)

best_formula2 <- as.formula(paste("Surv(tstart, tstop, event) ~ cgvhd:AI.I +",paste(best_variable_combination, collapse = " + ")))
final_model2 <- coxph(best_formula2, data = df_final)
os_summary2<-summary(final_model2)

CIBMTR_cox <- as.data.frame(cbind(os_summary$coefficients[,c(1,3,4,5)],os_summary$conf.int))
CIBMTR_cox$Variable<-rownames(CIBMTR_cox)
CIBMTR_cox$`lower .95`<-CIBMTR_cox$coef-2*CIBMTR_cox$`se(coef)`
CIBMTR_cox$`upper .95`<-CIBMTR_cox$coef+2*CIBMTR_cox$`se(coef)`
CIBMTR_cox<-CIBMTR_cox[order(CIBMTR_cox$coef),]
CIBMTR_cox$Variable<-factor(CIBMTR_cox$Variable,levels=CIBMTR_cox$Variable)
CIBMTR_cox_dead<-CIBMTR_cox

# Hazard plot of CIBMTR cox

# Assuming CIBMTR_cox is your existing DataFrame
CIBMTR_cox <- CIBMTR_cox %>%
  mutate(
    label_text = ifelse(`Pr(>|z|)` < 0.05, 
                        sprintf("bold('p = %.4f')", `Pr(>|z|)`),
                        sprintf("'p = %.4f'", `Pr(>|z|)`)
    )
  )
CIBMTR_cox$signif <- ifelse(CIBMTR_cox$`Pr(>|z|)` < 0.001, "***",
                            ifelse(CIBMTR_cox$`Pr(>|z|)` < 0.01, "**",
                                   ifelse(CIBMTR_cox$`Pr(>|z|)` < 0.05, "*", "")))

replacement_dict <- c("AI.I" = "Class-I\nautoimmune allele","AI.II" = "Class-II\nautoimmune allele", "graftype22" = "stem cell source", "Age.at.diagnosis" = "age at diagnosis",
                      'dnrage'='donor age','ipssrpr'='disease status at\ntime of transplant','dnrsex'='donor sex',
                      'condint'='conditioning intensity','donorgp'='donor type','ImprovedPHBR.II'='Improved PHBR-II',
                      'ImprovedPHBR.I'='Improved PHBR-I','cgvhd'='cGVHD','agvhd2'='aGVHD grade 2','agvhd34'='aGVHD grade 3-4',
                      'sex'='recipient sex')
replacement_dict <- c("AI.I" = "Class-I\nautoimmune allele","AI.II" = "Class-II\nautoimmune allele", "graftype22" = "stem cell source", "Age.at.diagnosis" = "age at diagnosis",
                      'dnrage'='donor age','ipssrpr'='disease status at\ntime of transplant','dnrsex'='donor sex',
                      'condint'='conditioning intensity','donorgp'='donor type','ImprovedPHBR.II'='Improved Donor\nAntigen Presentation',
                      'ImprovedPHBR.I'='Improved PHBR-I','cgvhd'='cGVHD','agvhd2'='aGVHD grade 2','agvhd34'='aGVHD grade 3-4',
                      'sex'='recipient sex')

CIBMTR_cox <- CIBMTR_cox %>%
  mutate(Variable = recode(Variable, !!!replacement_dict))

CIBMTR_cox <- CIBMTR_cox %>%
  mutate(across(c(`coef`, `lower .95`, `upper .95`), exp))

CIBMTR_cox<-CIBMTR_cox[CIBMTR_cox$Variable=='Improved Donor\nAntigen Presentation',]

# Plotting using ggplot2 with conditional bold labels using plotmath and correct color application
ggplot(data = CIBMTR_cox, aes(x = Variable, y = `coef`, ymin = `lower .95`, ymax = `upper .95`)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black", size=1) +
  geom_errorbar(aes(color = "green"), size=1.2, linetype=1, width=0.2) + # Apply color directly in geom_errorbar
  geom_point(aes(color = "green"), size=5, shape=15) + # Apply color directly in geom_point
  geom_label(aes(label = label_text), parse = TRUE, nudge_x = 0.27, nudge_y = 0) +  
  geom_text(aes(y = `coef`, label = signif), nudge_x = 0.42, nudge_y = 0.0,size=7) + coord_flip() + 
  theme_bw(base_size = 22) + 
  labs(y = "Hazard Ratio (OS)", x = "Clinical Variable") + 
  scale_color_manual(values = c("green" = "#0C7C59")) + # Ensure this matches with assigned colors
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black")) +
  ggtitle("Overall Survival")
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig2/OsHazzy_PHBR_only.pdf",width = 8,height = 3)


#########
## REL ##
#########

time_dep_df<-test_df_coxPH[test_df_coxPH$intxsurv>=2,]

df_time_dep <-
  tibble(Individual = time_dep_df$Individual,
         time = time_dep_df$intxcgvhd,
         cgvhd = time_dep_df$cgvhd)

df_time_dep_24 <-
  tibble(Individual = time_dep_df$Individual,
         time24 = time_dep_df$intxagvhd24,
         agvhd2 = time_dep_df$agvhd2)

df_time_dep_34 <-
  tibble(Individual = time_dep_df$Individual,
         time34 = time_dep_df$intxagvhd34,
         agvhd34 = time_dep_df$agvhd34)

df_ind <- tmerge(data1=time_dep_df[, c("intxsurv","intxrel","rel","Individual","condint","ipssrpr","donorgp","AI.I","dnrsex","dnrage","WorsePHBR.I","WorsePHBR.II","ImprovedPHBR.I","ImprovedPHBR.II",'AI.II','sex','graftype',"Age.at.diagnosis")], 
                 data2=time_dep_df[, c("intxsurv","intxrel","rel","Individual","condint","ipssrpr","donorgp","AI.I","dnrsex","dnrage","WorsePHBR.I","WorsePHBR.II","ImprovedPHBR.I","ImprovedPHBR.II",'AI.II','sex','graftype',"Age.at.diagnosis")], 
                 id=Individual, event=event(intxrel,rel)) 
df_final <-
  tmerge(data1=df_ind,
         data2=df_time_dep,
         id=Individual,
         cgvhd=tdc(time, cgvhd))
df_final <-
  tmerge(data1=df_final,
         data2=df_time_dep_24,
         id=Individual,
         agvhd2=tdc(time24, agvhd2))
df_final <-
  tmerge(data1=df_final,
         data2=df_time_dep_34,
         id=Individual,
         agvhd34=tdc(time34, agvhd34))

df_final$cgvhd[is.na(df_final$cgvhd)]<-0
df_final$agvhd2[is.na(df_final$agvhd2)]<-0
df_final$agvhd34[is.na(df_final$agvhd34)]<-0

# Define function to fit Cox model and find variable with the highest p-value
fit_cox_model <- function(data, variables) {
  print(variables)
  formula <- as.formula(paste("Surv(tstart, tstop, event) ~",paste(variables, collapse = " + ")))
  print(formula)
  model <- coxph(formula, data = data)
  summary_model <- summary(model)
  #print(summary_model)
  worst_var <- rownames(summary_model$coefficients)[which.max(summary_model$coefficients[, "Pr(>|z|)"])]
  aic <- AIC(model)
  return(list(aic = aic, worst_var = worst_var))
}

all_variables <- c('cgvhd',"condint","ipssrpr","donorgp","AI.I","dnrsex","dnrage","ImprovedPHBR.I","ImprovedPHBR.II","agvhd34","agvhd2",'AI.II','sex','graftype',"Age.at.diagnosis")

best_aic <- Inf
best_variable_combination <- NULL

# Iterate through all possible combinations to find the best model based on AIC
rerun <- TRUE
if (rerun) {
  repeat {
    res <- fit_cox_model(df_final, c(all_variables))
    aic <- res$aic
    worst_var <- res$worst_var
    if (aic < best_aic) {
      best_aic <- aic
      all_variables <- setdiff(all_variables, worst_var)
      best_variable_combination_rel <- all_variables
      print(aic)
      print(worst_var)
    } else {
      break
    }
  }
}

cat("Best variable combination:", best_variable_combination_rel, "\n")
cat("Best AIC:", best_aic, "\n")

best_variable_combination_rel<-c(best_variable_combination_rel)

# Fitting the best model
best_formula <- as.formula(paste("Surv(tstart, tstop, event) ~",paste(best_variable_combination_rel, collapse = " + ")))
final_model_REL <- coxph(best_formula, data = df_final)
rel_summary<-summary(final_model_REL)

best_formula2 <- as.formula(paste("Surv(tstart, tstop, event) ~ cgvhd:AI.I +",paste(best_variable_combination_rel, collapse = " + ")))
final_model2 <- coxph(best_formula2, data = df_final)
rel_summary2<-summary(final_model2)

#list of variables for later

CIBMTR_cox <- as.data.frame(cbind(rel_summary$coefficients[,c(1,3,4,5)],rel_summary$conf.int))
CIBMTR_cox$Variable<-rownames(CIBMTR_cox)
CIBMTR_cox$`lower .95`<-CIBMTR_cox$coef-2*CIBMTR_cox$`se(coef)`
CIBMTR_cox$`upper .95`<-CIBMTR_cox$coef+2*CIBMTR_cox$`se(coef)`
CIBMTR_cox<-CIBMTR_cox[order(CIBMTR_cox$coef),]
CIBMTR_cox$Variable<-factor(CIBMTR_cox$Variable,levels=CIBMTR_cox$Variable)
CIBMTR_cox_rel<-CIBMTR_cox
# Hazard plot of CIBMTR cox
# Assuming CIBMTR_cox is your existing DataFrame
CIBMTR_cox <- CIBMTR_cox %>%
  mutate(
    label_text = ifelse(`Pr(>|z|)` < 0.05, 
                        sprintf("bold('p = %.4f')", `Pr(>|z|)`),
                        sprintf("'p = %.4f'", `Pr(>|z|)`)
    )
  )
CIBMTR_cox$signif <- ifelse(CIBMTR_cox$`Pr(>|z|)` < 0.001, "***",
                            ifelse(CIBMTR_cox$`Pr(>|z|)` < 0.01, "**",
                                   ifelse(CIBMTR_cox$`Pr(>|z|)` < 0.05, "*", "")))

replacement_dict <- c("AI.I" = "Class-I\nautoimmune allele","AI.II" = "Class-II\nautoimmune allele", "graftype22" = "stem cell source", "Age.at.diagnosis" = "age at diagnosis",
                      'dnrage'='donor age','ipssrpr'='disease status at\ntime of transplant','dnrsex'='donor sex',
                      'condint'='conditioning intensity','donorgp'='donor type','ImprovedPHBR.II'='Improved PHBR-II',
                      'ImprovedPHBR.I'='Improved PHBR-I','cgvhd'='cGVHD','agvhd2'='aGVHD grade 2','agvhd34'='aGVHD grade 3-4',
                      'sex'='recipient sex')
CIBMTR_cox <- CIBMTR_cox %>%
  mutate(Variable = recode(Variable, !!!replacement_dict))

CIBMTR_cox <- CIBMTR_cox %>%
  mutate(across(c(`coef`, `lower .95`, `upper .95`), exp))

CIBMTR_cox<-CIBMTR_cox[CIBMTR_cox$Variable=='Class-I\nautoimmune allele',]

# Plotting using ggplot2 with conditional bold labels using plotmath and correct color application
ggplot(data = CIBMTR_cox, aes(x = Variable, y = `coef`, ymin = `lower .95`, ymax = `upper .95`)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black", size=1) +
  geom_errorbar(aes(color = "green"), size=1.2, linetype=1, width=0.2) + # Apply color directly in geom_errorbar
  geom_point(aes(color = "green"), size=5, shape=15) + # Apply color directly in geom_point
  geom_label(aes(label = label_text), parse = TRUE, nudge_x = 0.29, nudge_y = 0) +  
  geom_text(aes(y = `coef`, label = signif), nudge_x = 0.47, nudge_y = 0.0,size=7) + coord_flip() + 
  theme_bw(base_size = 22) +
  labs(y = "Hazard Ratio (Relapse Free Survival)", x = "Clinical Variable") +
  scale_color_manual(values = c("green" = "#0C7C59")) + # Ensure this matches with assigned colors
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black")) + ggtitle(label="Relapse-Free Survival")
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig2/RelHazzy_AI_only.pdf",width = 8,height = 3)

#########################
## TRIANGLE CORR PLOTS ##
#########################

library(corrplot)

# for each larger cohort:
all_variables <- c('cgvhd',"condint","ipssrpr","donorgp","AI.I","dnrsex","dnrage","ImprovedPHBR.I","ImprovedPHBR.II","agvhd34","agvhd2",'AI.II','sex','graftype',"Age.at.diagnosis")

#CIBMTR corplot
plot_df<-test_df_coxPH
plot_df<-plot_df[,c(all_variables)]

plot_df$graftype<-ifelse(plot_df$graftype==22,1,0)

# Calculate p-values for each pair of variables
p_values <- psych::corr.test(plot_df)$p
cor_matrix <- cor(plot_df)

replacement_dict <- c("AI.I" = "Class-I autoimmune allele","AI.II" = "Class-II autoimmune allele", "graftype22" = "stem cell source", "Age.at.diagnosis" = "age at diagnosis",
                      'dnrage'='donor age','ipssrpr'='status at transplant','dnrsex'='donor sex',
                      'condint'='conditioning intensity','donorgp'='donor type','ImprovedPHBR.II'='Improved PHBR-II',
                      'ImprovedPHBR.I'='Improved PHBR-I','cgvhd'='cGVHD','agvhd2'='aGVHD grade 2','agvhd34'='aGVHD grade 3-4',
                      'sex'='recipient sex')
colnames(cor_matrix) <- recode(colnames(cor_matrix), !!!replacement_dict)
rownames(cor_matrix) <- recode(rownames(cor_matrix), !!!replacement_dict)

colnames(p_values) <- recode(colnames(p_values), !!!replacement_dict)
rownames(p_values) <- recode(rownames(p_values), !!!replacement_dict)

#cor_matrix[cor_matrix==1]<-0

pdf("/Users/tjsears/Code/CIBMTR/FigsJun/FigS5/CorrPlot.pdf",width = 8,height = 8)

# Create a triangle correlation plot with p-values
corrplot(
  cor_matrix,diag=F,
  method = "color",
  type = "lower",col=rev(COL2('RdBu', 200)),
  #order = "hclust",
  tl.col = "black",mar=c(0,0,2,0),
  tl.srt = 45,tl.cex = 1.2,
  sig.level = c(0.001,0.01,0.05), # Set the significance level
  insig = "label_sig",outline=T,title = "Variable Correlation",
  p.mat = p_values,# Display p-values for all cells
  pch.cex = 1.7 # Adjust the size of p-value labels
  #addCoef.col = "black" # Set the color of p-value labels
)
dev.off()

####################################################
# plot pairwise interaction Pvalue of interactions #
####################################################

# for every combination of variables, set a coxph run where the interaction is included?
# use df_final from rel analysis

# Iterate through all possible combinations to find the best model based on AIC
get_cox_interaction <- function(data, variables) {
  
  #formula <- as.formula(paste("Surv(intxrel, rel) ~",paste(variables, collapse = "+"),"+",paste(variables, collapse = ":")))
  #formula <- as.formula(paste("Surv(intxrel, rel) ~",paste(variables, collapse = ":")))
  #formula <- as.formula(paste("Surv(tstart, tstop, event) ~",paste(variables, collapse = " + ","+",paste(variables, collapse = ":"))))
  formula <- as.formula(paste("Surv(tstart, tstop, event) ~",paste(variables, collapse = ":")))
  
  model <- coxph(formula, data = data)
  summary_model <- summary(model)
  worst_var <- rownames(summary_model$coefficients)[which.max(summary_model$coefficients[, "Pr(>|z|)"])]
  aic <- AIC(model)
  return(summary_model)
}

coxpval_df<-data.frame(matrix(nrow=length(best_variable_combination_rel),ncol=length(best_variable_combination_rel)))
rownames(coxpval_df)<-best_variable_combination_rel
colnames(coxpval_df)<-best_variable_combination_rel
coxpval_df[is.na(coxpval_df)]<-1

coxcoef_df<-data.frame(matrix(nrow=length(best_variable_combination_rel),ncol=length(best_variable_combination_rel)))
rownames(coxcoef_df)<-best_variable_combination_rel
colnames(coxcoef_df)<-best_variable_combination_rel
coxcoef_df[is.na(coxcoef_df)]<-0

for (i in best_variable_combination_rel){
  for (j in best_variable_combination_rel) {
    
    #summary<-get_cox_interaction(test_df_coxPH,c(i,j))
    summary<-get_cox_interaction(df_final,c(i,j))
    
    pval<-summary$coefficients[nrow(summary$coefficients),ncol(summary$coefficients)]
    coef<-summary$coefficients[nrow(summary$coefficients),1]
    
    coxpval_df[i,j]<-pval
    coxcoef_df[i,j]<-coef
  }
}

#triangle plot of direction and pval

coxcoef_df[is.na(coxcoef_df)]<-0
coxcoef_df<-data.matrix(coxcoef_df)

coxpval_df[is.na(coxpval_df)]<-0
coxpval_df<-data.matrix(coxpval_df)

corrplot(
  coxcoef_df,diag=F,
  method = "color",is.corr = FALSE,
  type = "lower",
  #order = "hclust",
  tl.col = "black",mar=c(0,0,2,0),
  tl.srt = 45,tl.cex = 1.2,
  sig.level = c(0.001,0.01,0.05), # Set the significance level
  insig = "label_sig",outline=T,title = "Variable Correlation",
  p.mat = coxpval_df,# Display p-values for all cells
  pch.cex = 1.7 # Adjust the size of p-value labels
  #addCoef.col = "black" # Set the color of p-value labels
)

# do a barplot version of this probably...


# Feature importance
cols=brewer.pal(8,"Dark2")
#cols=cols[c(1,8,2,4)]
cols=c("#5AB1BB","#A5C882","#F7DD72")
library(reshape2)

replacement_dict <- c("AI.I" = "class-I autoimmune","AI.II" = "class-II autoimmune", "graftype22" = "stem cell source", "Age.at.diagnosis" = "age at diagnosis",
                      'dnrage'='donor age','ipssrpr'='status at transplant','dnrsex'='donor sex',
                      'condint'='conditioning intensity','donorgp'='donor type','ImprovedPHBR.II'='Improved PHBR-II',
                      'ImprovedPHBR.I'='Improved PHBR-I','cgvhd'='cGVHD','agvhd2'='aGVHD grade 2','agvhd34'='aGVHD grade 3-4',
                      'sex'='recipient sex')
colnames(coxpval_df) <- recode(colnames(coxpval_df), !!!replacement_dict)
rownames(coxpval_df) <- recode(rownames(coxpval_df), !!!replacement_dict)

feet_long<-melt(coxpval_df)
feet_long$AltVar<-colnames(coxpval_df)

# switch var 1 and var 2

feet_long$FinalVar<-paste(feet_long$Var1,"/",feet_long$Var2)
feet_long<-feet_long[order(feet_long$value,decreasing = F),]
feet_long<-feet_long[!(feet_long$Var1==feet_long$Var2),]
feet_long = feet_long[seq(1, nrow(feet_long), 2), ]

feet_long$value<-p.adjust(feet_long$value,method="bonferroni")
feet_long$value<-abs(log2(feet_long$value))

feet_long<-feet_long[c(1:13),]

feet_long$FinalVar<-factor(feet_long$FinalVar,levels=c(feet_long$FinalVar))

feet_long$InteractionType<-c("MHC Autoimmune / Autoimmune","Clinical / Autoimmune","Clinical / Autoimmune",
                             "Clinical / Autoimmune","Clinical / Clinical","Clinical / Autoimmune",
                             "Clinical / Autoimmune","Clinical / Clinical","Clinical / Autoimmune",
                             "Clinical / Clinical", "Clinical / Clinical","Clinical / Clinical","Clinical / Clinical")

ggplot(feet_long, aes(x = FinalVar, y = value, fill = InteractionType)) + ylab(label="Log2 Interaction Pvalue")+
  geom_bar(stat = "identity") + theme_minimal(base_size = 16) + xlab(label=NULL)+
  scale_fill_manual(values = cols) + theme(axis.text.x = element_text(angle = 55, vjust = 1.07, hjust=1),
                                           plot.margin = margin(t = 10, r = 20, b = 10, l = 50, unit = "pt"), # Adjust plot margins
                                           panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(yintercept=4.321928, linetype="dashed", color = "red") + geom_label(label="Bonferroni < 0.05",x=11,y=5.4,show.legend = F,inherit.aes=F)
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig3/SimplifiedInteraction.pdf",width = 10,height = 5)





#########################
# Relapse-only analysis #
#########################


time_dep_df<-test_df_coxPH#[test_df_coxPH$intxsurv>=2,]

#convert relapse-free survival into just relapse

# if patient dies at the exact time of relapse, change relapse death to simple censorship
time_dep_df$rel<-ifelse((time_dep_df$intxrel==time_dep_df$intxsurv)&time_dep_df$dead==1,0,time_dep_df$rel)

df_time_dep <-
  tibble(Individual = time_dep_df$Individual,
         time = time_dep_df$intxcgvhd,
         cgvhd = time_dep_df$cgvhd)

df_time_dep_24 <-
  tibble(Individual = time_dep_df$Individual,
         time24 = time_dep_df$intxagvhd24,
         agvhd2 = time_dep_df$agvhd2)

df_time_dep_34 <-
  tibble(Individual = time_dep_df$Individual,
         time34 = time_dep_df$intxagvhd34,
         agvhd34 = time_dep_df$agvhd34)

df_ind <- tmerge(data1=time_dep_df[, c("intxsurv","intxrel","rel","Individual","condint","ipssrpr","donorgp","AI.I","dnrsex","dnrage","WorsePHBR.I","WorsePHBR.II","ImprovedPHBR.I","ImprovedPHBR.II",'AI.II','sex','graftype',"Age.at.diagnosis")], 
                 data2=time_dep_df[, c("intxsurv","intxrel","rel","Individual","condint","ipssrpr","donorgp","AI.I","dnrsex","dnrage","WorsePHBR.I","WorsePHBR.II","ImprovedPHBR.I","ImprovedPHBR.II",'AI.II','sex','graftype',"Age.at.diagnosis")], 
                 id=Individual, event=event(intxrel,rel)) 
df_final <-
  tmerge(data1=df_ind,
         data2=df_time_dep,
         id=Individual,
         cgvhd=tdc(time, cgvhd))
df_final <-
  tmerge(data1=df_final,
         data2=df_time_dep_24,
         id=Individual,
         agvhd2=tdc(time24, agvhd2))
df_final <-
  tmerge(data1=df_final,
         data2=df_time_dep_34,
         id=Individual,
         agvhd34=tdc(time34, agvhd34))

df_final$cgvhd[is.na(df_final$cgvhd)]<-0
df_final$agvhd2[is.na(df_final$agvhd2)]<-0
df_final$agvhd34[is.na(df_final$agvhd34)]<-0

# Define function to fit Cox model and find variable with the highest p-value
fit_cox_model <- function(data, variables) {
  print(variables)
  formula <- as.formula(paste("Surv(tstart, tstop, event) ~",paste(variables, collapse = " + ")))
  print(formula)
  model <- coxph(formula, data = data)
  summary_model <- summary(model)
  #print(summary_model)
  worst_var <- rownames(summary_model$coefficients)[which.max(summary_model$coefficients[, "Pr(>|z|)"])]
  aic <- AIC(model)
  return(list(aic = aic, worst_var = worst_var))
}

all_variables <- c('cgvhd',"condint","ipssrpr","donorgp","AI.I","dnrsex","dnrage","ImprovedPHBR.I","ImprovedPHBR.II","agvhd34","agvhd2",'AI.II','sex','graftype',"Age.at.diagnosis")

best_aic <- Inf
best_variable_combination <- NULL

# Iterate through all possible combinations to find the best model based on AIC
rerun <- TRUE
if (rerun) {
  repeat {
    res <- fit_cox_model(df_final, c(all_variables))
    aic <- res$aic
    worst_var <- res$worst_var
    if (aic < best_aic) {
      best_aic <- aic
      all_variables <- setdiff(all_variables, worst_var)
      best_variable_combination_rel <- all_variables
      print(aic)
      print(worst_var)
    } else {
      break
    }
  }
}

cat("Best variable combination:", best_variable_combination_rel, "\n")
cat("Best AIC:", best_aic, "\n")

best_variable_combination_rel<-c(best_variable_combination_rel)

# Fitting the best model
best_formula <- as.formula(paste("Surv(tstart, tstop, event) ~",paste(best_variable_combination_rel, collapse = " + ")))
final_model_REL <- coxph(best_formula, data = df_final)
rel_summary<-summary(final_model_REL)

best_formula2 <- as.formula(paste("Surv(tstart, tstop, event) ~ cgvhd:AI.I +",paste(best_variable_combination_rel, collapse = " + ")))
final_model2 <- coxph(best_formula2, data = df_final)
rel_summary2<-summary(final_model2)

#list of variables for later

CIBMTR_cox <- as.data.frame(cbind(rel_summary$coefficients[,c(1,3,4,5)],rel_summary$conf.int))
CIBMTR_cox$Variable<-rownames(CIBMTR_cox)
CIBMTR_cox$`lower .95`<-CIBMTR_cox$coef-2*CIBMTR_cox$`se(coef)`
CIBMTR_cox$`upper .95`<-CIBMTR_cox$coef+2*CIBMTR_cox$`se(coef)`
CIBMTR_cox<-CIBMTR_cox[order(CIBMTR_cox$coef),]
CIBMTR_cox$Variable<-factor(CIBMTR_cox$Variable,levels=CIBMTR_cox$Variable)

# Hazard plot of CIBMTR cox
# Assuming CIBMTR_cox is your existing DataFrame
CIBMTR_cox <- CIBMTR_cox %>%
  mutate(
    label_text = ifelse(`Pr(>|z|)` < 0.05, 
                        sprintf("bold('p = %.4f')", `Pr(>|z|)`),
                        sprintf("'p = %.4f'", `Pr(>|z|)`)
    )
  )
CIBMTR_cox$signif <- ifelse(CIBMTR_cox$`Pr(>|z|)` < 0.001, "***",
                            ifelse(CIBMTR_cox$`Pr(>|z|)` < 0.01, "**",
                                   ifelse(CIBMTR_cox$`Pr(>|z|)` < 0.05, "*", "")))

replacement_dict <- c("AI.I" = "Class-I\nautoimmune allele","AI.II" = "Class-II\nautoimmune allele", "graftype22" = "stem cell source", "Age.at.diagnosis" = "age at diagnosis",
                      'dnrage'='donor age','ipssrpr'='disease status at\ntime of transplant','dnrsex'='donor sex',
                      'condint'='conditioning intensity','donorgp'='donor type','ImprovedPHBR.II'='Improved PHBR-II',
                      'ImprovedPHBR.I'='Improved PHBR-I','cgvhd'='cGVHD','agvhd2'='aGVHD grade 2','agvhd34'='aGVHD grade 3-4',
                      'sex'='recipient sex')
CIBMTR_cox <- CIBMTR_cox %>%
  mutate(Variable = recode(Variable, !!!replacement_dict))

CIBMTR_cox <- CIBMTR_cox %>%
  mutate(across(c(`coef`, `lower .95`, `upper .95`), exp))

# Plotting using ggplot2 with conditional bold labels using plotmath and correct color application
ggplot(data = CIBMTR_cox, aes(x = Variable, y = `coef`, ymin = `lower .95`, ymax = `upper .95`)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black", size=1) +
  geom_errorbar(aes(color = "green"), size=1.2, linetype=1, width=0.2) + # Apply color directly in geom_errorbar
  geom_point(aes(color = "green"), size=5, shape=15) + # Apply color directly in geom_point
  geom_label(aes(label = label_text), parse = TRUE, nudge_x = 0.27, nudge_y = 0) +  
  geom_text(aes(y = `coef`, label = signif), nudge_x = 0.42, nudge_y = 0,size=7) + coord_flip() + 
  theme_bw(base_size = 22) +
  labs(y = "Hazard Ratio (Relapse)", x = "Clinical Variable") +
  scale_color_manual(values = c("green" = "#0C7C59")) + # Ensure this matches with assigned colors
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black")) + ggtitle(label="Relapse")
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig2/RelapseHazzy.pdf",width = 8,height = 12)



#################
# TRM  analysis #
#################

# 

time_dep_df<-test_df_coxPH[test_df_coxPH$intxsurv>2,] # 
time_dep_df$intxsurv<-ifelse(time_dep_df$rel==1,time_dep_df$intxrel,time_dep_df$intxsurv)
time_dep_df$dead<-ifelse(time_dep_df$rel==1,0,time_dep_df$dead)

df_time_dep <-
  tibble(Individual = time_dep_df$Individual,
         time = time_dep_df$intxcgvhd,
         cgvhd = time_dep_df$cgvhd)

df_time_dep_24 <-
  tibble(Individual = time_dep_df$Individual,
         time24 = time_dep_df$intxagvhd24,
         agvhd2 = time_dep_df$agvhd2)

df_time_dep_34 <-
  tibble(Individual = time_dep_df$Individual,
         time34 = time_dep_df$intxagvhd34,
         agvhd34 = time_dep_df$agvhd34)

df_ind <- tmerge(data1=time_dep_df[, c("intxsurv","dead","Individual","condint","ipssrpr","donorgp","AI.I","dnrsex","dnrage","WorsePHBR.I","WorsePHBR.II","ImprovedPHBR.I","ImprovedPHBR.II",'AI.II','sex','graftype',"Age.at.diagnosis")], 
                 data2=time_dep_df[, c("intxsurv","dead","Individual","condint","ipssrpr","donorgp","AI.I","dnrsex","dnrage","WorsePHBR.I","WorsePHBR.II","ImprovedPHBR.I","ImprovedPHBR.II",'AI.II','sex','graftype',"Age.at.diagnosis")], 
                 id=Individual, event=event(intxsurv,dead)) 
df_final <-
  tmerge(data1=df_ind,
         data2=df_time_dep,
         id=Individual,
         cgvhd=tdc(time, cgvhd))
df_final <-
  tmerge(data1=df_final,
         data2=df_time_dep_24,
         id=Individual,
         agvhd2=tdc(time24, agvhd2))
df_final <-
  tmerge(data1=df_final,
         data2=df_time_dep_34,
         id=Individual,
         agvhd34=tdc(time34, agvhd34))

df_final$cgvhd[is.na(df_final$cgvhd)]<-0
df_final$agvhd2[is.na(df_final$agvhd2)]<-0
df_final$agvhd34[is.na(df_final$agvhd34)]<-0

# Define function to fit Cox model and find variable with the highest p-value
fit_cox_model <- function(data, variables) {
  print(variables)
  formula <- as.formula(paste("Surv(tstart, tstop, event) ~",paste(variables, collapse = " + ")))
  print(formula)
  model <- coxph(formula, data = data)
  summary_model <- summary(model)
  #print(summary_model)
  worst_var <- rownames(summary_model$coefficients)[which.max(summary_model$coefficients[, "Pr(>|z|)"])]
  aic <- AIC(model)
  return(list(aic = aic, worst_var = worst_var))
}

all_variables <- c("condint","ipssrpr","donorgp","AI.I","dnrsex","dnrage","ImprovedPHBR.I","ImprovedPHBR.II","agvhd34","agvhd2",'AI.II','sex','graftype',"Age.at.diagnosis")

best_aic <- Inf
best_variable_combination <- NULL

# Iterate through all possible combinations to find the best model based on AIC
rerun <- TRUE
if (rerun) {
  repeat {
    res <- fit_cox_model(df_final, c(all_variables))
    aic <- res$aic
    worst_var <- res$worst_var
    if (aic < best_aic) {
      best_aic <- aic
      all_variables <- setdiff(all_variables, worst_var)
      best_variable_combination <- all_variables
      print(aic)
      print(worst_var)
    } else {
      break
    }
  }
}

cat("Best variable combination:", best_variable_combination, "\n")
cat("Best AIC:", best_aic, "\n")

# Fitting the best model
best_formula <- as.formula(paste("Surv(tstart, tstop, event) ~",paste(best_variable_combination, collapse = " + ")))
final_model <- coxph(best_formula, data = df_final)
os_summary<-summary(final_model)

best_formula2 <- as.formula(paste("Surv(tstart, tstop, event) ~ cgvhd:AI.I +",paste(best_variable_combination, collapse = " + ")))
final_model2 <- coxph(best_formula2, data = df_final)
os_summary2<-summary(final_model2)

CIBMTR_cox <- as.data.frame(cbind(os_summary$coefficients[,c(1,3,4,5)],os_summary$conf.int))
CIBMTR_cox$Variable<-rownames(CIBMTR_cox)
CIBMTR_cox$`lower .95`<-CIBMTR_cox$coef-2*CIBMTR_cox$`se(coef)`
CIBMTR_cox$`upper .95`<-CIBMTR_cox$coef+2*CIBMTR_cox$`se(coef)`
CIBMTR_cox<-CIBMTR_cox[order(CIBMTR_cox$coef),]
CIBMTR_cox$Variable<-factor(CIBMTR_cox$Variable,levels=CIBMTR_cox$Variable)

# Hazard plot of CIBMTR cox
# Assuming CIBMTR_cox is your existing DataFrame
CIBMTR_cox <- CIBMTR_cox %>%
  mutate(
    label_text = ifelse(`Pr(>|z|)` < 0.05, 
                        sprintf("bold('p = %.4f')", `Pr(>|z|)`),
                        sprintf("'p = %.4f'", `Pr(>|z|)`)
    )
  )

CIBMTR_cox$signif <- ifelse(CIBMTR_cox$`Pr(>|z|)` < 0.001, "***",
                            ifelse(CIBMTR_cox$`Pr(>|z|)` < 0.01, "**",
                                   ifelse(CIBMTR_cox$`Pr(>|z|)` < 0.05, "*", "")))

replacement_dict <- c("AI.I" = "Class-I\nautoimmune allele","AI.II" = "Class-II\nautoimmune allele", "graftype22" = "stem cell source", "Age.at.diagnosis" = "age at diagnosis",
                      'dnrage'='donor age','ipssrpr'='disease status at\ntime of transplant','dnrsex'='donor sex',
                      'condint'='conditioning intensity','donorgp'='donor type','ImprovedPHBR.II'='Improved PHBR-II',
                      'ImprovedPHBR.I'='Improved PHBR-I','cgvhd'='cGVHD','agvhd2'='aGVHD grade 2','agvhd34'='aGVHD grade 3-4',
                      'sex'='recipient sex')
CIBMTR_cox <- CIBMTR_cox %>%
  mutate(Variable = recode(Variable, !!!replacement_dict))

CIBMTR_cox <- CIBMTR_cox %>%
  mutate(across(c(`coef`, `lower .95`, `upper .95`), exp))

# Plotting using ggplot2 with conditional bold labels using plotmath and correct color application
ggplot(data = CIBMTR_cox, aes(x = Variable, y = `coef`, ymin = `lower .95`, ymax = `upper .95`)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black", size=1) +
  geom_errorbar(aes(color = "green"), size=1.2, linetype=1, width=0.2) + # Apply color directly in geom_errorbar
  geom_point(aes(color = "green"), size=5, shape=15) + # Apply color directly in geom_point
  geom_label(aes(label = label_text), parse = TRUE, nudge_x = 0.27, nudge_y = 0) +  
  geom_text(aes(y = `coef`, label = signif), nudge_x = 0.42, nudge_y = 0,size=7) + coord_flip() + 
  theme_bw(base_size = 22) +
  labs(y = "Hazard Ratio (TRM)", x = "Clinical Variable") +
  scale_color_manual(values = c("green" = "#0C7C59")) + # Ensure this matches with assigned colors
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black")) + ggtitle(label="TRM")
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig2/TRMHazzy.pdf",width = 8,height = 12)



###############
# Alt metrics #
###############

# quick univariate analysis of alternative PHBR metrics
install.packages("plotrix")
library(plotrix)

# load in alt metrics
alt_phbr<-read.table('/Users/tjsears/Code/CIBMTR/Tables/PHBR_alt_metrics.txt',sep='\t',header=T)
alt_phbr_plot<-merge(alt_phbr,test_df_coxPH,by.x='DID',by.y='Individual')
alt_phbr_plot$PHBR.II_diff_mean<-1-alt_phbr_plot$PHBR.II_diff_mean
alt_phbr_plot$PHBR.II_diff_sum<-1-alt_phbr_plot$PHBR.II_diff_sum

# change names of alt metrics

alt_metrics<-c('PHBR.I_best','I_diff_2_sum','I_diff_2_mean','TotalPresentable_I',
               'PHBR.II_best','PHBR.II_diff_mean','PHBR.II_diff_sum','TotalPresentable_II')

columns<-c("coef"     ,  "se(coef)",   "z"      ,    "Pr(>|z|)" ,  "exp(coef)",  "exp(-coef)", "lower .95" , "upper .95" , "Variable"  ,
          "label_text" ,"signif"    )

cox_df <- setNames(data.frame(matrix(ncol = length(columns), nrow = 0)), columns)

for (alt_var in alt_metrics){

  # run loop of relevant alt metrics
  formula<-as.formula(paste("Surv(intxsurv.x,dead.x) ~",alt_var))
  final_model2 <- coxph(formula,data = alt_phbr_plot)
  os_summary2<-summary(final_model2)
  
  CIBMTR_cox <- as.data.frame(cbind(t(as.data.frame(os_summary2$coefficients[,c(1,3,4,5)])),os_summary2$conf.int))
  CIBMTR_cox$Variable<-alt_var
  CIBMTR_cox$`lower .95`<-CIBMTR_cox$coef-2*CIBMTR_cox$`se(coef)`
  CIBMTR_cox$`upper .95`<-CIBMTR_cox$coef+2*CIBMTR_cox$`se(coef)`
  CIBMTR_cox<-CIBMTR_cox[order(CIBMTR_cox$coef),]
  CIBMTR_cox$Variable<-factor(CIBMTR_cox$Variable,levels=CIBMTR_cox$Variable)
  
  # Hazard plot of CIBMTR cox
  # Assuming CIBMTR_cox is your existing DataFrame
  CIBMTR_cox <- CIBMTR_cox %>%
    mutate(
      label_text = ifelse(`Pr(>|z|)` < 0.05, 
                          sprintf("bold('p = %.4f')", `Pr(>|z|)`),
                          sprintf("'p = %.4f'", `Pr(>|z|)`)
      )
    )
  
  CIBMTR_cox$signif <- ifelse(CIBMTR_cox$`Pr(>|z|)` < 0.001, "***",
                              ifelse(CIBMTR_cox$`Pr(>|z|)` < 0.01, "**",
                                     ifelse(CIBMTR_cox$`Pr(>|z|)` < 0.05, "*", "")))
  
  replacement_dict <- c('PHBR.I_best'='Best presented PHBR-I',
                        #'PHBR.I_diff_sum'='Sum of PHBR-I difference',
                        'I_diff_2_sum'='Sum of total PHBR-I',
                        'I_diff_2_mean'='Mean of total PHBR-I',
                        'TotalPresentable_I'='Number of PHBR-I\nneoantigens',
                        'PHBR.II_best'='Best presented PHBR-II',
                        #'PHBR.I_diff_sum'='Sum of PHBR-I difference',
                        'PHBR.II_diff_sum'='Sum of total PHBR-II',
                        'PHBR.II_diff_mean'='Mean of total PHBR-II',
                        'TotalPresentable_II'='Number of PHBR-II\nneoantigens')
  
  CIBMTR_cox <- CIBMTR_cox %>%
    mutate(Variable = recode(Variable, !!!replacement_dict))
  
  CIBMTR_cox <- CIBMTR_cox %>%
    mutate(across(c(`coef`, `lower .95`, `upper .95`), exp))
  
  #print(CIBMTR_cox)
  
  cox_df<-rbind(cox_df,CIBMTR_cox)
  
}

# Plotting using ggplot2 with conditional bold labels using plotmath and correct color application
ggplot(data = cox_df, aes(x = Variable, y = `coef`, ymin = `lower .95`, ymax = `upper .95`)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black", size=1) +
  geom_errorbar(aes(color = "green"), size=1.2, linetype=1, width=0.2) + # Apply color directly in geom_errorbar
  geom_point(aes(color = "green"), size=5, shape=15) + # Apply color directly in geom_point
  geom_label(aes(label = label_text), parse = TRUE, nudge_x = 0.27, nudge_y = 0) +  
  geom_text(aes(y = `coef`, label = signif), nudge_x = 0.42, nudge_y = 0,size=7) + coord_flip() + ylim(0,3.5) +
  theme_bw(base_size = 22) +
  labs(y = "Hazard Ratio (OS)", x = "") +
  scale_color_manual(values = c("green" = "#0C7C59")) + # Ensure this matches with assigned colors
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black")) + ggtitle(label="Overall Survival") +
  annotate("segment", x= 3, xend = 3, y=(cox_df$`lower .95`[cox_df$Variable == "Mean of total PHBR-I"]), 
           yend = 3.5,#max(cox_df$`upper .95`[cox_df$Variable == "Mean of total PHBR-I"]), 
           arrow = arrow(length = unit(0.5, "cm")),size=1, color = "#0C7C59")
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig2/overallSurvival_alt_phbr_Hazzy.pdf",width = 8,height = 12)

cox_df <- setNames(data.frame(matrix(ncol = length(columns), nrow = 0)), columns)

for (alt_var in alt_metrics){
  
  # run loop of relevant alt metrics
  formula<-as.formula(paste("Surv(intxrel.x,rel.x) ~",alt_var))
  final_model2 <- coxph(formula,data = alt_phbr_plot)
  os_summary2<-summary(final_model2)
  
  CIBMTR_cox <- as.data.frame(cbind(t(as.data.frame(os_summary2$coefficients[,c(1,3,4,5)])),os_summary2$conf.int))
  CIBMTR_cox$Variable<-alt_var
  CIBMTR_cox$`lower .95`<-CIBMTR_cox$coef-2*CIBMTR_cox$`se(coef)`
  CIBMTR_cox$`upper .95`<-CIBMTR_cox$coef+2*CIBMTR_cox$`se(coef)`
  CIBMTR_cox<-CIBMTR_cox[order(CIBMTR_cox$coef),]
  CIBMTR_cox$Variable<-factor(CIBMTR_cox$Variable,levels=CIBMTR_cox$Variable)
  
  # Hazard plot of CIBMTR cox
  # Assuming CIBMTR_cox is your existing DataFrame
  CIBMTR_cox <- CIBMTR_cox %>%
    mutate(
      label_text = ifelse(`Pr(>|z|)` < 0.05, 
                          sprintf("bold('p = %.4f')", `Pr(>|z|)`),
                          sprintf("'p = %.4f'", `Pr(>|z|)`)
      )
    )
  
  CIBMTR_cox$signif <- ifelse(CIBMTR_cox$`Pr(>|z|)` < 0.001, "***",
                              ifelse(CIBMTR_cox$`Pr(>|z|)` < 0.01, "**",
                                     ifelse(CIBMTR_cox$`Pr(>|z|)` < 0.05, "*", "")))
  
  replacement_dict <- c('PHBR.I_best'='Best presented PHBR-I',
                        #'PHBR.I_diff_sum'='Sum of PHBR-I difference',
                        'I_diff_2_sum'='Sum of total PHBR-I',
                        'I_diff_2_mean'='Mean of total PHBR-I',
                        'TotalPresentable_I'='Number of PHBR-I\nneoantigens',
                        'PHBR.II_best'='Best presented PHBR-II',
                        #'PHBR.I_diff_sum'='Sum of PHBR-I difference',
                        'PHBR.II_diff_sum'='Sum of total PHBR-II',
                        'PHBR.II_diff_mean'='Mean of total PHBR-II',
                        'TotalPresentable_II'='Number of PHBR-II\nneoantigens')
  
  
  CIBMTR_cox <- CIBMTR_cox %>%
    mutate(Variable = recode(Variable, !!!replacement_dict))
  
  CIBMTR_cox <- CIBMTR_cox %>%
    mutate(across(c(`coef`, `lower .95`, `upper .95`), exp))
  
  #print(CIBMTR_cox)
  
  cox_df<-rbind(cox_df,CIBMTR_cox)
  
}

# Plotting using ggplot2 with conditional bold labels using plotmath and correct color application
ggplot(data = cox_df, aes(x = Variable, y = `coef`, ymin = `lower .95`, ymax = `upper .95`)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black", size=1) +
  geom_errorbar(aes(color = "green"), size=1.2, linetype=1, width=0.2) + # Apply color directly in geom_errorbar
  geom_point(aes(color = "green"), size=5, shape=15) + # Apply color directly in geom_point
  geom_label(aes(label = label_text), parse = TRUE, nudge_x = 0.27, nudge_y = 0) +  
  geom_text(aes(y = `coef`, label = signif), nudge_x = 0.42, nudge_y = 0,size=7) + coord_flip() + ylim(0,3.5) +
  theme_bw(base_size = 22) +
  labs(y = "Hazard Ratio (RFS)",x='') +
  scale_color_manual(values = c("green" = "#0C7C59")) + # Ensure this matches with assigned colors
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black")) + ggtitle(label="Relapse-Free Survival") +
  annotate("segment", x= 3, xend = 3, y=(cox_df$`lower .95`[cox_df$Variable == "Mean of total PHBR-I"]), 
           yend = 3.5,#max(cox_df$`upper .95`[cox_df$Variable == "Mean of total PHBR-I"]), 
           arrow = arrow(length = unit(0.5, "cm")),size=1, color = "#0C7C59")
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig2/relapse_free_alt_phbr_Hazzy.pdf",width = 8,height = 12)

