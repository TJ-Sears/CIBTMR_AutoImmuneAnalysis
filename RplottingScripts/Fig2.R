# Fig 2 CIBMTR

library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(dplyr)

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

cols<-brewer.pal(n=3,"Dark2")
# COHORT CHARACTERISTICS

CIBMTR_df<-read.table("/Users/tjsears/Code/CIBMTR/Tables/complete_dataset_feb25.txt",sep="\t",header=T,row.names = 1)
CIBMTR_df[is.na(CIBMTR_df)] <-0# fill all other na with 0
#CIBMTR_df = CIBMTR_df('.', 0, regex=False)
CIBMTR_df$dnrsex[CIBMTR_df$dnrsex=="."]<-1
CIBMTR_df$dnrsex<-as.numeric(CIBMTR_df$dnrsex)
#CIBMTR_df=CIBMTR_df[CIBMTR_df$condint!=99,] # keep or drop?
#CIBMTR_df=CIBMTR_df[CIBMTR_df$donorgp!=99,]
CIBMTR_df$AutoImmune_Allele<-ifelse(CIBMTR_df$AI==1,"Present","Absent")
# Fraction of AI alleles

########################
# HLA Allele Bar plots #
########################

# read in summary table from python
class_I<-c('HLA-B27:05','HLA-B35:01',
           'HLA-B50:01','HLA-B51:01',#'HLA-B08:01',
           'HLA-B07:02',
           #'HLA-B07:04',
           'HLA-C05:01',
           'HLA-B13:02',
           'HLA-C12:03','HLA-A29:02')
class_II<-c(#'DRB1_1501',# >30%
            'DRB1_0801', 'DRB1_1302',#'DRB1_1601',#
            'DRB1_0301',
            'DRB1_0202',
            'DRB1_0803',#'DRB1_1001',
            'DRB1_0302'
            ,#'DRB1_0303',#'DRB1_0304',
            'DRB1_1402', 'DRB1_1401',
            'DQB10201', 'DQB10202')

HLA_table<-read.table("/Users/tjsears/Code/CIBMTR/Tables/HLA_table.txt",sep="\t",header=T)

HLA_A<-HLA_table[,c("Individual","HLA.A1")]
colnames(HLA_A)<-c("Individual","HLA.A2")
HLA_A<-rbind(HLA_A,HLA_table[,c("Individual","HLA.A2")])
colnames(HLA_A)<-c("Individual","HLA.A1")

table_hla<-table(HLA_A$HLA.A1)
table_hla_a_freq<-table_hla/(sum(table_hla))
HLA_A$HLA.A1<-factor(HLA_A$HLA.A1,levels=c(names(table_hla)[order(table_hla,decreasing = T)]))

HLA_A$AI<-ifelse(HLA_A$HLA.A1%in%class_I,"AI+","AI-")
ggplot(HLA_A, aes(x=HLA.A1,fill=AI)) +
  geom_bar() + 
  ggtitle("HLA A freq. and AI status") +
  xlab("Values") + theme_bw(base_size = 14) +
  ylab("Frequency") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
#ggsave("/Users/tjsears/Code/CIBMTR/Figs_Feb/Fig2/HLA_A.pdf",width = 9,height = 7)

HLA_A<-HLA_table[,c("Individual","HLA.B1")]
colnames(HLA_A)<-c("Individual","HLA.B2")
HLA_A<-rbind(HLA_A,HLA_table[,c("Individual","HLA.B2")])
colnames(HLA_A)<-c("Individual","HLA.A1")

table_hla<-table(HLA_A$HLA.A1)
table_hla_b_freq<-table_hla/(sum(table_hla))

HLA_A$HLA.A1<-factor(HLA_A$HLA.A1,levels=c(names(table_hla)[order(table_hla,decreasing = T)]))

HLA_A$AI<-ifelse(HLA_A$HLA.A1%in%class_I,"AI+","AI-")
ggplot(HLA_A, aes(x=HLA.A1,fill=AI)) +
  geom_bar() + 
  ggtitle("HLA B freq. and AI status") +
  xlab("Values") + theme_bw(base_size = 14) +
  ylab("Frequency") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
ggsave("/Users/tjsears/Code/CIBMTR/Figs_Feb/Fig2/HLA_B.pdf",width = 11,height = 7)

HLA_A<-HLA_table[,c("Individual","HLA.C1")]
colnames(HLA_A)<-c("Individual","HLA.C2")
HLA_A<-rbind(HLA_A,HLA_table[,c("Individual","HLA.C2")])
colnames(HLA_A)<-c("Individual","HLA.A1")

table_hla<-table(HLA_A$HLA.A1)
table_hla_c_freq<-table_hla/(sum(table_hla))

HLA_A$HLA.A1<-factor(HLA_A$HLA.A1,levels=c(names(table_hla)[order(table_hla,decreasing = T)]))

HLA_A$AI<-ifelse(HLA_A$HLA.A1%in%class_I,"AI+","AI-")
ggplot(HLA_A, aes(x=HLA.A1,fill=AI)) +
  geom_bar() + 
  ggtitle("HLA C freq. and AI status") +
  xlab("Values") + theme_bw(base_size = 14) +
  ylab("Frequency") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
ggsave("/Users/tjsears/Code/CIBMTR/Figs_Feb/Fig2/HLA_C.pdf",width = 9,height = 7)


# MHC-II

HLA_A<-HLA_table[,c("Individual","HLA.DQA")]
colnames(HLA_A)<-c("Individual","HLA.DQB")
#HLA_A<-rbind(HLA_A,HLA_table[,c("Individual","HLA.DRBB")])
#HLA_A<-rbind(HLA_A,HLA_table[,c("Individual","HLA.DRBB")])
HLA_A<-rbind(HLA_A,HLA_table[,c("Individual","HLA.DQB")])

colnames(HLA_A)<-c("Individual","HLA.A1")

table_hla<-table(HLA_A$HLA.A1)
table_hla_II_freq<-table_hla/(sum(table_hla))

HLA_A$HLA.A1<-factor(HLA_A$HLA.A1,levels=c(names(table_hla)[order(table_hla,decreasing = T)]))

HLA_A$AI<-ifelse(HLA_A$HLA.A1%in%class_II,"AI+","AI-")
ggplot(HLA_A, aes(x=HLA.A1,fill=AI)) +
  geom_bar() + 
  ggtitle("HLA Class II freq. and AI status") +
  xlab("Values") + theme_bw(base_size = 14) +
  ylab("Frequency") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels
ggsave("/Users/tjsears/Code/CIBMTR/Figs_Feb/Fig2/HLA_II.pdf",width = 9,height = 7)


# Look at AI freq in TCGA?

# all I want is simple bar plot showing comparison b/w
#load TCGA HLA types
class_i<-read.table("/Users/tjsears/Code/CIBMTR/Tables/all_mhc_i_types.txt",sep="\t",header=T)
class_ii<-read.table("/Users/tjsears/Code/CIBMTR/Tables/all_mhc_ii_types.csv",sep=",",header=T)


#load TCGA clinical data
clin_data<-read.table("/Users/tjsears/Code/CIBMTR/Tables/Liu2018.TCGA_survival.csv",sep=",",header=T)

#load in addtl TCGA data suggesting more severe disease
clin_data_mas<-read.table("/Users/tjsears/Code/CIBMTR/Tables/nationwidechildrens.org_clinical_patient_laml.txt",sep="\t",header=T)
colnames(clin_data_mas)<-clin_data_mas[1,]
clin_data_mas<-clin_data_mas[3:nrow(clin_data_mas),]

clin_data_i<-merge(clin_data,class_i,by.x="bcr_patient_barcode",by.y="X")
clin_data_i_ii<-merge(clin_data_i,class_ii,by.x="bcr_patient_barcode",by.y="X")
#clin_data_i_ii<-merge(clin_data_i_ii,clin_data_mas,by="bcr_patient_barcode")

clin_data_i_ii_AML<-clin_data_i_ii#[clin_data_i_ii$type=="READ",]
clin_data_i_ii_AML$AI.I<-rowSums(data.frame(clin_data_i_ii_AML$HLA.A1%in%class_I,
                                  clin_data_i_ii_AML$HLA.A2%in%class_I,
                                  clin_data_i_ii_AML$HLA.B1%in%class_I,
                                  clin_data_i_ii_AML$HLA.B2%in%class_I,
                                  clin_data_i_ii_AML$HLA.C1%in%class_I,
                                  clin_data_i_ii_AML$HLA.C2%in%class_I))

clin_data_i_ii_AML$AI.II<-rowSums(data.frame(clin_data_i_ii_AML$DRB1_1%in%class_I,
                                   clin_data_i_ii_AML$DRB1_2%in%class_II))

#bar chart TCGA
# Transforming Data to Long Format
plot_df<-clin_data_i_ii_AML[,c("AI.I","AI.II")]
data_long <- plot_df %>% 
  mutate(id = row_number()) %>%
  pivot_longer(-id, names_to = "class", values_to = "value") %>%
  group_by(class, value) %>%
  summarise(count = n(), .groups = 'drop')

data_long$value[data_long$value>=3]<-"3+"

data_long$class[data_long$class=='AI.I']<-'Class-I'
data_long$class[data_long$class=='AI.II']<-'Class-II'

# Plotting
ggplot(data_long, aes(x = class, y = count, fill = factor(value))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("#ADFFE8","#00CC8F","#007A56","#00523A")) +
  labs(x = "HLA Class", y = "Fraction of Pts", fill = "No. of Alleles") +
  theme_bw(base_size = 18) + theme(panel.grid.major = element_blank(), 
                                   panel.grid.minor = element_blank())+
  ggtitle("TCGA\nAutoimmune Alleles")
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig1/TCGA_AI_alleles.pdf",width = 6.5,height = 6)


# Transforming Data to Long Format
plot_df<-CIBMTR_df[,c("classI_AI","classII_AI")]
data_long <- plot_df %>% 
  mutate(id = row_number()) %>%
  pivot_longer(-id, names_to = "class", values_to = "value") %>%
  group_by(class, value) %>%
  summarise(count = n(), .groups = 'drop')

data_long$value[data_long$value>=3]<-"3+"

data_long$class[data_long$class=='classI_AI']<-'Class-I'
data_long$class[data_long$class=='classII_AI']<-'Class-II'

# Plotting
ggplot(data_long, aes(x = class, y = count, fill = factor(value))) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values = c("#ADFFE8","#00CC8F","#007A56","#00523A")) +
  labs(x = "HLA Class", y = "Fraction of Pts", fill = "No. of Alleles") +
  theme_bw(base_size = 18) + theme(panel.grid.major = element_blank(), 
                                   panel.grid.minor = element_blank())+
  ggtitle("CIBMTR Autoimmune Alleles")
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig1/AI_alleles.pdf",width = 6.5,height = 6)

# qucik maf
sum(0.034*0.69*0.86, # A
    0.965*0.31*0.86, # B
    0.965*0.69*0.14, # C
    0.034*0.31*0.86, #A B
    0.965*0.31*0.14, #B C
    0.034*0.31*0.14, # A B C
    0.034*0.69*0.14) # A C

#0.965*0.965*

#km curve of TCGA AI allele
#remove pts that are unlikely to undergo bone marrow transplant

clin_data_KM<-clin_data_i_ii_AML[!clin_data_i_ii_AML$OS.time=="MISSING",]
clin_data_KM$OS.time<-as.numeric(clin_data_KM$OS.time)
clin_data_KM$OS<-as.numeric(clin_data_KM$OS)

#clin_data_KM<-clin_data_KM[clin_data_KM$age_at_initial_pathologic_diagnosis.x<70,]
#clin_data_KM<-clin_data_KM[clin_data_KM$acute_myeloid_leukemia_calgb_cytogenetics_risk_category=="Poor"
#                           |clin_data_KM$acute_myeloid_leukemia_calgb_cytogenetics_risk_category=="Intermediate/Normal",]
#
#clin_data_KM<-clin_data_KM[!(clin_data_KM$leukemia_french_american_british_morphology_code=="M0 Undifferentiated"
#                             |clin_data_KM$leukemia_french_american_british_morphology_code=="M1"
#),]

clin_data_KM

clin_data_KM$AI<-ifelse(clin_data_KM$AI.I+clin_data_KM$AI.II>0,1,0)

km_trt_fit <- survfit(Surv(OS.time, OS) ~ AI, data=clin_data_KM)
res<- pairwise_survdiff(Surv(OS.time, OS) ~ AI, data=clin_data_KM,p.adjust.method = "none",rho=0) #, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

p1<-ggsurvplot(km_trt_fit,data=clin_data_KM,size=1.2,censor.shape="|", censor.size = 2,legend=c(0.2,0.15),
               pval=F,xlim=c(0,2000),break.x.by=200,ncensor.plot.height=0.25,
               font.x=18,font.y=18,font.legend=18,font.title=18,
               legend.labs=c("AI-","AI+"),palette=cols[2:1],legend.title="",
               ,xlab="Time in Days",ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))+ggtitle(label="Autoimmune Alleles & Overall Survival")
p1$plot+ggplot2::annotate("label",x=600,y = 0.8,size=6, label = paste("P =",round(res$p.value,4)))

ggsave("/Users/tjsears/Code/CIBMTR/Figs_Feb/Fig2/TCGA_AI_split.pdf",width = 6,height = 6)


temp<-clin_data_KM[clin_data_KM$type=="COAD"|clin_data_KM$type=="READ",]
temp$age_at_initial_pathologic_diagnosis<-as.numeric(temp$age_at_initial_pathologic_diagnosis)
all_composite <- summary(coxph(Surv(OS.time, OS) ~ AI.II+gender+age_at_initial_pathologic_diagnosis, data = temp))
all_composite


km_trt_fit <- survfit(Surv(intxrel, rel) ~ AutoImmune_Allele, data=CIBMTR_df)
res<- pairwise_survdiff(Surv(intxrel, rel) ~ AutoImmune_Allele, data=CIBMTR_df,p.adjust.method = "none",rho=0) #, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

coxph(Surv(intxrel, rel) ~ AutoImmune_Allele, data=CIBMTR_df)

cols<-brewer.pal(n=3,"Dark2")

p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df,size=1.2,censor.shape="|", censor.size = 6,legend=c("none"), #changed from top for abstact
               pval=F,xlim=c(0,65),break.x.by=10,ncensor.plot.height=0.25,risk.table.title="Number at Risk",
               font.x=18,font.y=18,font.legend=14,font.title=18,#risk.table = c("absolute"),
               palette=cols[2:1],legend.title="",
               ,xlab="Time in Months",ylim=c(0,1),ylab="Relapse-Free Survival (%)",surv.scale = c("percent"))+ggtitle(label="Autoimmune Alleles &\nRelapse-Free Survival")
p1$plot+ggplot2::annotate("label",x=50,y = 0.8,size=6, label = paste("P =",round(res$p.value,4)))
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig1/AI_split.pdf",width = 5,height = 4) #changed from 7 and 6

cols<-brewer.pal(n=3,"Dark2")

p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),
               pval=F,xlim=c(0,65),break.x.by=10,ncensor.plot.height=0.25,risk.table.title="Number at Risk",
               font.x=18,font.y=18,font.legend=10,font.title=18,risk.table = c("absolute"),
               palette=cols[2:1],legend.title="",
               ,xlab="Time in Months",ylim=c(0,1),ylab="Relapse-Free Survival (%)",surv.scale = c("percent"))#+ggtitle(label="Autoimmune Alleles & Relapse-Free Survival")
p1#$plot+annotate(geom = "table", x = 140, y = 0.1, label = (as.data.frame(res$p.value)))

ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig1/AI_split_pval.pdf",width = 8,height = 2)


km_trt_fit <- survfit(Surv(intxsurv, dead) ~ AutoImmune_Allele, data=CIBMTR_df)
res<- pairwise_survdiff(Surv(intxsurv, dead) ~ AutoImmune_Allele, data=CIBMTR_df,p.adjust.method = "none",rho=0) #, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

coxph(Surv(intxsurv, dead) ~ AutoImmune_Allele, data=CIBMTR_df)

cols<-brewer.pal(n=3,"Dark2")

p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),
               pval=F,xlim=c(0,65),break.x.by=10,ncensor.plot.height=0.25,risk.table.title="Number at Risk",
               font.x=18,font.y=18,font.legend=10,font.title=18,risk.table = c("absolute"),
               palette=cols[2:1],legend.title="",
               ,xlab="Time in Months",ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))#+ggtitle(label="Autoimmune Alleles & Relapse-Free Survival")
p1#$plot+annotate(geom = "table", x = 140, y = 0.1, label = (as.data.frame(res$p.value)))

ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig1/AI_split_OS_pval.pdf",width = 8,height = 2)


p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),
               pval=F,xlim=c(0,65),break.x.by=10,ncensor.plot.height=0.25,risk.table.title="Number at Risk",
               font.x=18,font.y=18,font.legend=14,font.title=18,#risk.table = c("absolute"),
               palette=cols[2:1],legend.title="",
               ,xlab="Time in Months",ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))+ggtitle(label="Autoimmune Alleles & Overall Survival")
p1$plot+ggplot2::annotate("label",x=50,y = 0.8,size=6, label = paste("P =",round(res$p.value,4)))

ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig1/AI_split_OS.pdf",width = 7,height = 6)







km_trt_fit <- survfit(Surv(intxrel, rel) ~ AI.I, data=CIBMTR_df)
res<- pairwise_survdiff(Surv(intxrel, rel) ~ AI.I, data=CIBMTR_df,p.adjust.method = "none",rho=0) #, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

cols<-brewer.pal(n=3,"Dark2")

p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),
               pval=F,xlim=c(0,65),break.x.by=10,ncensor.plot.height=0.25,risk.table.title="Number at Risk",
               font.x=18,font.y=18,font.legend=14,font.title=18,#risk.table = c("absolute"),
               palette=cols[2:1],legend.title="",
               ,xlab="Time in Months",ylim=c(0,1),ylab="Relapse-Free Survival (%)",surv.scale = c("percent"))+ggtitle(label="Autoimmune Alleles & Relapse-Free Survival")
p1$plot+ggplot2::annotate("label",x=50,y = 0.8,size=6, label = paste("P =",round(res$p.value,4)))
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig1/AI_split_AI_I.pdf",width = 7,height = 6)

cols<-brewer.pal(n=3,"Dark2")

p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),
               pval=F,xlim=c(0,65),break.x.by=10,ncensor.plot.height=0.25,risk.table.title="Number at Risk",
               font.x=18,font.y=18,font.legend=10,font.title=18,risk.table = c("absolute"),
               palette=cols[2:1],legend.title="",
               ,xlab="Time in Months",ylim=c(0,1),ylab="Relapse-Free Survival (%)",surv.scale = c("percent"))#+ggtitle(label="Autoimmune Alleles & Relapse-Free Survival")
p1#$plot+annotate(geom = "table", x = 140, y = 0.1, label = (as.data.frame(res$p.value)))

ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig1/AI_split_pval_AI_I.pdf",width = 8,height = 2)


km_trt_fit <- survfit(Surv(intxsurv, dead) ~ AI.I, data=CIBMTR_df)
res<- pairwise_survdiff(Surv(intxsurv, dead) ~ AI.I, data=CIBMTR_df,p.adjust.method = "none",rho=0) #, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

cols<-brewer.pal(n=3,"Dark2")

p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),
               pval=F,xlim=c(0,65),break.x.by=10,ncensor.plot.height=0.25,risk.table.title="Number at Risk",
               font.x=18,font.y=18,font.legend=10,font.title=18,risk.table = c("absolute"),
               palette=cols[2:1],legend.title="",
               ,xlab="Time in Months",ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))#+ggtitle(label="Autoimmune Alleles & Relapse-Free Survival")
p1#$plot+annotate(geom = "table", x = 140, y = 0.1, label = (as.data.frame(res$p.value)))

ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig1/AI_split_OS_pval_AI_I.pdf",width = 8,height = 2)


p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),
               pval=F,xlim=c(0,65),break.x.by=10,ncensor.plot.height=0.25,risk.table.title="Number at Risk",
               font.x=18,font.y=18,font.legend=14,font.title=18,#risk.table = c("absolute"),
               palette=cols[2:1],legend.title="",
               ,xlab="Time in Months",ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))+ggtitle(label="Autoimmune Alleles & Overall Survival")
p1$plot+ggplot2::annotate("label",x=50,y = 0.8,size=6, label = paste("P =",round(res$p.value,4)))

ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig1/AI_split_OS_AI_I.pdf",width = 7,height = 6)






km_trt_fit <- survfit(Surv(intxrel, rel) ~ AI.II, data=CIBMTR_df)
res<- pairwise_survdiff(Surv(intxrel, rel) ~ AI.II, data=CIBMTR_df,p.adjust.method = "none",rho=0) #, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res


cols<-brewer.pal(n=3,"Dark2")

p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),
               pval=F,xlim=c(0,65),break.x.by=10,ncensor.plot.height=0.25,risk.table.title="Number at Risk",
               font.x=18,font.y=18,font.legend=14,font.title=18,#risk.table = c("absolute"),
               palette=cols[2:1],legend.title="",
               ,xlab="Time in Months",ylim=c(0,1),ylab="Relapse-Free Survival (%)",surv.scale = c("percent"))+ggtitle(label="Autoimmune Alleles & Relapse-Free Survival")
p1$plot+ggplot2::annotate("label",x=50,y = 0.8,size=6, label = paste("P =",round(res$p.value,4)))
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig1/AI_split_AI_II.pdf",width = 7,height = 6)

cols<-brewer.pal(n=3,"Dark2")

p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),
               pval=F,xlim=c(0,65),break.x.by=10,ncensor.plot.height=0.25,risk.table.title="Number at Risk",
               font.x=18,font.y=18,font.legend=10,font.title=18,risk.table = c("absolute"),
               palette=cols[2:1],legend.title="",
               ,xlab="Time in Months",ylim=c(0,1),ylab="Relapse-Free Survival (%)",surv.scale = c("percent"))#+ggtitle(label="Autoimmune Alleles & Relapse-Free Survival")
p1#$plot+annotate(geom = "table", x = 140, y = 0.1, label = (as.data.frame(res$p.value)))

ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig1/AI_split_pval_AI_II.pdf",width = 8,height = 2)


km_trt_fit <- survfit(Surv(intxsurv, dead) ~ AI.II, data=CIBMTR_df)
res<- pairwise_survdiff(Surv(intxsurv, dead) ~ AI.II, data=CIBMTR_df,p.adjust.method = "none",rho=0) #, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

cols<-brewer.pal(n=3,"Dark2")

p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),
               pval=F,xlim=c(0,65),break.x.by=10,ncensor.plot.height=0.25,risk.table.title="Number at Risk",
               font.x=18,font.y=18,font.legend=10,font.title=18,risk.table = c("absolute"),
               palette=cols[2:1],legend.title="",
               ,xlab="Time in Months",ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))#+ggtitle(label="Autoimmune Alleles & Relapse-Free Survival")
p1#$plot+annotate(geom = "table", x = 140, y = 0.1, label = (as.data.frame(res$p.value)))

ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig1/AI_split_OS_pval_AI_II.pdf",width = 8,height = 2)


p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),
               pval=F,xlim=c(0,65),break.x.by=10,ncensor.plot.height=0.25,risk.table.title="Number at Risk",
               font.x=18,font.y=18,font.legend=14,font.title=18,#risk.table = c("absolute"),
               palette=cols[2:1],legend.title="",
               ,xlab="Time in Months",ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"))+ggtitle(label="Autoimmune Alleles & Overall Survival")
p1$plot+ggplot2::annotate("label",x=50,y = 0.8,size=6, label = paste("P =",round(res$p.value,4)))

ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig1/AI_split_OS_AI_II.pdf",width = 7,height = 6)


