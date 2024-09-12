# CIBMTR Discovery

################
# REL ANALYSIS #
################
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
library(ggplot2)
library(RColorBrewer)
library(ggpmisc)

cols<-brewer.pal(n = 4,name = "Dark2")

data<-read.table("/Users/tjsears/Code/CIBMTR/Tables/complete_dataset_Feb25.txt",sep="\t",header=T,row.names = 1)
CIBMTR_df<-data

#fix missingness
CIBMTR_df[is.na(CIBMTR_df)] <-0# fill all other na with 0
CIBMTR_df$dnrsex[CIBMTR_df$dnrsex=="."]<-1
CIBMTR_df$dnrsex<-as.numeric(CIBMTR_df$dnrsex)
CIBMTR_df=CIBMTR_df[CIBMTR_df$condint!=99,]
CIBMTR_df=CIBMTR_df[CIBMTR_df$donorgp!=99,]
CIBMTR_df$CGVHD<-CIBMTR_df$cgvhd
CIBMTR_df$intxrel<-CIBMTR_df$intxrel*30.437
CIBMTR_df$intxsurv<-CIBMTR_df$intxsurv*30.437

CIBMTR_df$WorsePHBR.I<-CIBMTR_df$PHBR.I.change
CIBMTR_df$WorsePHBR.II<-CIBMTR_df$PHBR.II.change
CIBMTR_df$ImprovedPHBR.I<--1*CIBMTR_df$Poor..High..PHBR.I
CIBMTR_df$ImprovedPHBR.II<--1*CIBMTR_df$Poor..High..PHBR.II

CIBMTR_df$CGVHD<-ifelse(CIBMTR_df$CGVHD>0,"Yes","No")
CIBMTR_df$AI.I<-ifelse(CIBMTR_df$AI.I>0,"Yes","No")

#CIBMTR_df<-CIBMTR_df[CIBMTR_df$intxsurv>=270,]
#CIBMTR_df<-CIBMTR_df[CIBMTR_df$AI.I=='No',]
quantile(CIBMTR_df$intxcgvhd)

CIBMTR_df$ImprovedPHBR.II_bin<-ifelse(CIBMTR_df$ImprovedPHBR.II>0,1,CIBMTR_df$ImprovedPHBR.II)
CIBMTR_df$ImprovedPHBR.II_bin<-ifelse(CIBMTR_df$ImprovedPHBR.II<0,2,CIBMTR_df$ImprovedPHBR.II_bin)

km_trt_fit <- survfit(Surv(intxsurv, dead) ~ ImprovedPHBR.II_bin, data=CIBMTR_df)
res<- pairwise_survdiff(Surv(intxsurv, dead) ~ ImprovedPHBR.II_bin, data=CIBMTR_df,p.adjust.method = "none",rho=0)#, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),
               pval=F,xlim=c(0,2000),break.x.by=200,ncensor.plot.height=0.25,#risk.table.title="",
               font.x=18,font.y=18,font.legend=13,font.title=18,#risk.table = c("absolute"),
               palette=cols,legend.title="",title="CIBMTR Class-I Autoimmune Alleles & cGVHD\nRelapse-Free Survival (9mo Overall Survival Landmark)"
               ,xlab="Time in Days",ylim=c(0,1),ylab="Overall Survival (%)",surv.scale = c("percent"),
               risk.table.title="Number at Risk")
p1#$plot+annotate(geom = "table", x = 140, y = 0.1, label = (as.data.frame(res$p.value)))
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/FigS6/phbr_test_os.pdf",width = 10,height = 7)


p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),
               pval=F,xlim=c(0,2000),break.x.by=200,ncensor.plot.height=0.25,#risk.table.title="",
               font.x=18,font.y=18,font.legend=10,font.title=18,risk.table = c("absolute"),
               palette=cols,legend.title="",title="CIBMTR Class-I Autoimmune Alleles & cGVHD\nRelapse-Free Survival"
               ,xlab="Time in Days",ylim=c(0,1),ylab="Relapse-Free Survival (%)",surv.scale = c("percent"),
                risk.table.title="Number at Risk")
p1#$plot+annotate(geom = "table", x = 140, y = 0.1, label = (as.data.frame(res$p.value)))
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/FigS6/discc_km_risktable.pdf",width =11,height = 2)


library(gt)

# 4.160000e-10	 
# 4.6e-11
df<-data.frame(res$p.value)

df[3,1]<-"<0.0001"
df[3,3]<-"<0.0001"
df[3,2]<-"<0.0001"

colnames(df)<-c("AI.I=No, CGVHD=No ","AI.I=No, CGVHD=Yes ","AI.I=Yes, CGVHD=No ")
gt_table<-gt(df,rownames_to_stub = T,rowname_col = cols,auto_align =F) %>% # Color the header (column labels)
  tab_style(
    style = list(
      cell_fill(color = cols[1]),
      cell_text(color = "white")  # Adjust text color for better readability
    ),
    locations = cells_column_labels(columns = 2)
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = cols[2]),
      cell_text(color = "white")  # Adjust text color for better readability
    ),
    locations = cells_column_labels(columns = 3)
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = cols[3]),
      cell_text(color = "white")  # Adjust text color for better readability
    ),
    locations = cells_column_labels(columns = 4)
  ) %>%
  tab_style(
    style = list(cell_fill(color = cols[2]),cell_text(color = "white")),
    locations = cells_stub(rows = 1)
  )%>%
  tab_style(
    style = list(cell_fill(color = cols[3]),cell_text(color = "white")),
    locations = cells_stub(rows = 2)
  )%>%
  tab_style(
    style = list(cell_fill(color = cols[4]),cell_text(color = "white")),
    locations = cells_stub(rows = 3)
  )
 
gt_table

################
# BOOSTRAPPING #
################

set.seed(42)

# Bootstrap analysis of variables
CIBMTR_df$Stratification<-'Neither'
CIBMTR_df$Stratification[CIBMTR_df$CGVHD=='Yes']<-'CGVHD'
CIBMTR_df$Stratification[CIBMTR_df$AI.I=='Yes']<-'AI.I+'
CIBMTR_df$Stratification[CIBMTR_df$AI.I=='Yes'&CIBMTR_df$CGVHD=='Yes']<-'CGVHD and AI.I+'
CIBMTR_df$Stratification<-factor(CIBMTR_df$Stratification,levels=c('Neither','CGVHD','AI.I+','CGVHD and AI.I+'))

# Set the desired total number of rows
desired_size = 1000

# Sampling indices with replacement
sample_indices = sample(1:nrow(CIBMTR_df), size = desired_size, replace = TRUE)

# Creating the new DataFrame with duplicated rows
CIBMTR_df_boot = CIBMTR_df[sample_indices, ]

km_trt_fit <- survfit(Surv(intxrel, rel) ~ CGVHD+AI.I, data=CIBMTR_df_boot)
res<- pairwise_survdiff(Surv(intxrel, rel) ~ CGVHD+AI.I, data=CIBMTR_df_boot,p.adjust.method = "none",rho=0)#, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df_boot,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),
               pval=F,xlim=c(0,2000),break.x.by=200,ncensor.plot.height=0.25,#risk.table.title="",
               font.x=18,font.y=18,font.legend=13,font.title=18,#risk.table = c("absolute"),
               palette=cols,legend.title="",title="Bootstrapped (n=1000)\nCIBMTR Class-I Autoimmune Alleles & cGVHD\nRelapse-Free Survival (9mo Overall Survival Landmark)"
               ,xlab="Time in Days",ylim=c(0,1),ylab="Relapse-Free Survival (%)",surv.scale = c("percent"),
               risk.table.title="Number at Risk")
p1#$plot+annotate(geom = "table", x = 140, y = 0.1, label = (as.data.frame(res$p.value)))
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/FigS6/disc_km_split_boot.pdf",width = 10,height = 7)


p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df_boot,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),
               pval=F,xlim=c(0,2000),break.x.by=200,ncensor.plot.height=0.25,#risk.table.title="",
               font.x=18,font.y=18,font.legend=10,font.title=18,risk.table = c("absolute"),
               palette=cols,legend.title="",title="Bootstrapped (n=1000)\nCIBMTR Class-I Autoimmune Alleles & cGVHD\nRelapse-Free Survival"
               ,xlab="Time in Days",ylim=c(0,1),ylab="Relapse-Free Survival (%)",surv.scale = c("percent"),
               risk.table.title="Number at Risk")
p1#$plot+annotate(geom = "table", x = 140, y = 0.1, label = (as.data.frame(res$p.value)))
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/FigS6/discc_km_risktable_boot.pdf",width =11,height = 2)

library(gt)

# 4.160000e-10	 
# 4.6e-11
df<-data.frame(res$p.value)

df[3,1]<-"<0.0001"
df[3,3]<-"<0.0001"
df[3,2]<-"<0.0001"

colnames(df)<-c("AI.I=No, CGVHD=No ","AI.I=No, CGVHD=Yes ","AI.I=Yes, CGVHD=No ")
gt_table<-gt(df,rownames_to_stub = T,rowname_col = cols,auto_align =F) %>% # Color the header (column labels)
  tab_style(
    style = list(
      cell_fill(color = cols[1]),
      cell_text(color = "white")  # Adjust text color for better readability
    ),
    locations = cells_column_labels(columns = 2)
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = cols[2]),
      cell_text(color = "white")  # Adjust text color for better readability
    ),
    locations = cells_column_labels(columns = 3)
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = cols[3]),
      cell_text(color = "white")  # Adjust text color for better readability
    ),
    locations = cells_column_labels(columns = 4)
  ) %>%
  tab_style(
    style = list(cell_fill(color = cols[2]),cell_text(color = "white")),
    locations = cells_stub(rows = 1)
  )%>%
  tab_style(
    style = list(cell_fill(color = cols[3]),cell_text(color = "white")),
    locations = cells_stub(rows = 2)
  )%>%
  tab_style(
    style = list(cell_fill(color = cols[4]),cell_text(color = "white")),
    locations = cells_stub(rows = 3)
  )

gt_table

