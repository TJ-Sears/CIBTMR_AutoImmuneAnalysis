# Venn diagram CIBMTR


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


# read in summary table from python
class_I<-c('HLA-B27:05','HLA-B35:01',
           'HLA-B50:01','HLA-B51:01',
           'HLA-B07:02',
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
#ggsave("/Users/tjsears/Code/CIBMTR/Figs_Feb/Fig2/HLA_II.pdf",width = 9,height = 7)


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

clin_data_i_ii_AML$AI.II<-rowSums(data.frame(clin_data_i_ii_AML$DRB1_1%in%class_II,
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
#ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig1/AI_alleles.pdf",width = 6.5,height = 6)


# convert HLA table into allele based yes/no categories
HLA_table_plot<-HLA_table
HLA_table_plot$A1<-ifelse(HLA_table$HLA.A1%in%class_I,'AI+','AI-')
HLA_table_plot$A2<-ifelse(HLA_table$HLA.A2%in%class_I,'AI+','AI-')
HLA_table_plot$B1<-ifelse(HLA_table$HLA.B1%in%class_I,'AI+','AI-')
HLA_table_plot$B2<-ifelse(HLA_table$HLA.B2%in%class_I,'AI+','AI-')
HLA_table_plot$C1<-ifelse(HLA_table$HLA.C1%in%class_I,'AI+','AI-')
HLA_table_plot$C2<-ifelse(HLA_table$HLA.C2%in%class_I,'AI+','AI-')

HLA_table_plot<-HLA_table_plot[,c('A1','A2','B1','B2','C1','C2')]
HLA_table_plot$Patient<-seq(1,494)
HLA_table_plot_long<-pivot_longer(HLA_table_plot,cols=c('A1','A2','B1','B2','C1','C2'))
colnames(HLA_table_plot_long)<-c('Patient','Allele','Status')
HLA_table_plot_long$Allele<-as.factor(HLA_table_plot_long$Allele)
HLA_table_plot_long$Status<-as.factor(HLA_table_plot_long$Status)

ggplot(HLA_table_plot_long, aes(x = Allele, fill = Status)) +
  geom_bar(stat = "count") +
  labs(title = "Class-I Autoimmune Allele by Locus",
       x = "HLA Gene",
       y = "Patients with Autoimmune Allele") +
  theme_bw(base_size = 16) + scale_fill_manual(values = cols[c(3,1)]) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig1/AI_alleles_loci.pdf",width = 6.5,height = 6)

# calculate probability of possession of each number of alleles


# draw red box around 1+ cutoff?
data_long_prob<-data_long
data_long_prob$prob<-data_long_prob$count/494
data_long_prob<-data_long_prob[data_long_prob$class=='Class-I',]

ggplot(data_long_prob, aes(x = value, y=prob,fill=class)) +
  geom_bar(stat = "identity") +
  labs(title = "Class-I Autoimmune Allele Probability",
       x = "Number of Class-I Autoimmune Alleles",
       y = "Proportion of Patients") +
  theme_bw(base_size = 16) + scale_fill_manual(values = cols[c(1)]) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position='none') +
  # Add rectangle annotation around the first three columns
  annotate("rect", xmin = 1.5, xmax = 4.5, ymin = -0.01, ymax = 0.55, 
               alpha = 0.2, color = "red", fill = NA, size = 1) +
  annotate("label", x = 3.4, y = 0.35, label = "71% of patients have 1+\nClass-I autoimmune allele", color = "black", size = 5)
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig1/AI_alleles_probability.pdf",width = 6.5,height = 6)










