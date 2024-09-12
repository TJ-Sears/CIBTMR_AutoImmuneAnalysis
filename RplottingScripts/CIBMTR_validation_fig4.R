# CIBMTR Validation Andrea

data<-read.table("/Users/tjsears/Code/CIBMTR/OriginalCohort/TJ_dat_test_feb25.tsv",sep="\t",header=T)
data$intxrel<-data$RelapseFree_days
data$rel<-data$Relapsed.1

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

CIBMTR_df<-data
CIBMTR_df$Glucksberg.Grading.of.aGVHD_cleaned[is.na(CIBMTR_df$Glucksberg.Grading.of.aGVHD_cleaned)]<-1
CIBMTR_df$AI.I<-CIBMTR_df$recipient
CIBMTR_df$CGVHD<-CIBMTR_df$Chronic.GVHD..Y.N.

CIBMTR_df$CGVHD<-ifelse(CIBMTR_df$CGVHD>0,"Yes","No")
CIBMTR_df$AI.I<-ifelse(CIBMTR_df$AI.I>0,"Yes","No")
CIBMTR_df<-CIBMTR_df[CIBMTR_df$OS_days>100,]

km_trt_fit <- survfit(Surv(intxrel, rel) ~ AI.I+CGVHD, data=CIBMTR_df)
res<- pairwise_survdiff(Surv(intxrel, rel) ~ AI.I+CGVHD, data=CIBMTR_df,p.adjust.method = "none",rho=0)#, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),risk.table=F,
               pval=F,xlim=c(0,2200),break.x.by=200,ncensor.plot.height=0.25,risk.table.title="Number at Risk",
               font.x=18,font.y=18,font.legend=13,font.title=18,#risk.table = c("absolute"),
               #risk.table.fontsize=3, risk.table.y.text=3,
               palette=cols,legend.title="",title="UCSD Class-I Autoimmune Alleles & cGVHD\nRelapse-Free Survival",
               ,xlab="Time in Days",ylim=c(0,1),ylab="Relapse-Free Survival (%)",surv.scale = c("percent"))
p1
ggsave("/Users/tjsears/Code/CIBMTR/FigsMay/Fig4/Val_AI_REL_km_split.pdf",width = 10,height = 7)

p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),risk.table=F,
               pval=F,xlim=c(0,2200),break.x.by=200,ncensor.plot.height=0.25,risk.table.title="Number at Risk",
               font.x=18,font.y=18,font.legend=10,font.title=18,risk.table = c("absolute"),
               #risk.table.fontsize=3, risk.table.y.text=3,
               palette=cols,legend.title="",title="UCSD Class-I Autoimmune Alleles & cGVHD\nRelapse-Free Survival",
               ,xlab="Time in Days",ylim=c(0,1),ylab="Relapse-Free Survival (%)",surv.scale = c("percent"))
p1

ggsave("/Users/tjsears/Code/CIBMTR/FigsMay/Fig4/Val_AI_REL_km_split_riskTable.pdf",width = 9,height = 2)


library(gt)

# 4.160000e-10	 
# 4.6e-11
df<-data.frame(res$p.value)

#df[3,1]<-"<0.0001"
#df[3,3]<-"<0.0001"

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
  


