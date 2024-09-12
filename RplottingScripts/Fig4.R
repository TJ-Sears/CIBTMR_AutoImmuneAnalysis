# CIBMTR Fig 3 


library(ggplot2)
library(RColorBrewer)

# read in CIBMTR data and explore time dependent covariate

CIBMTR_df<-read.table("/Users/tjsears/Code/CIBMTR/Tables/complete_dataset_feb5.txt",sep="\t",header=T,row.names = 1)
CIBMTR_df[is.na(CIBMTR_df)] <-0# fill all other na with 0
#CIBMTR_df = CIBMTR_df('.', 0, regex=False)
CIBMTR_df$dnrsex[CIBMTR_df$dnrsex=="."]<-1
CIBMTR_df$dnrsex<-as.numeric(CIBMTR_df$dnrsex)
#CIBMTR_df=CIBMTR_df[CIBMTR_df$condint!=99,] # keep or drop?
#CIBMTR_df=CIBMTR_df[CIBMTR_df$donorgp!=99,]

#rename PHBR variables
CIBMTR_df$Poor.Donor.PHBRI<-CIBMTR_df$Poor..High..PHBR.I
CIBMTR_df$Poor.Donor.PHBRII<-CIBMTR_df$Poor..High..PHBR.II

CIBMTR_df$Improved.Donor.PHBRI<-CIBMTR_df$PHBR.I.change
CIBMTR_df$Improved.Donor.PHBRII<-CIBMTR_df$PHBR.II.change


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

# taking top variables from initial python analysis 
regression_vars<-c("GVHD","ipssrpr","donorgp","AI","dnrsex","dnrage",'AI_plus','ClassI_difference','ClassII_difference',
                   "Poor..High..PHBR.I","Poor..High..PHBR.II","PHBR.I.change","PHBR.II.change")

#hazzy plot ICI
CIBMTR_df$AI_plus[CIBMTR_df$AI_plus>=3]<-3
#CIBMTR_df$AI_plus<-as.factor(CIBMTR_df$AI_plus)

plot_df<-CIBMTR_df
#plot_df$AI<-ifelse(plot_df$AI_plus>=1,1,0)

cols2<-brewer.pal(3,"Dark2")

#plot_df<-plot_df[plot_df$intxsurv>=2,]
plot_df_high<-plot_df[plot_df$cgvhd==1,]
plot_df_mid<-plot_df[plot_df$cgvhd==0,]

all_composite <- summary(coxph(Surv(intxrel, rel) ~ Poor..High..PHBR.I, data = plot_df_high))
all_mid <- summary(coxph(Surv(intxrel, rel) ~ Poor..High..PHBR.I, data = plot_df_mid))
#all_somatic <- summary(coxph(Surv(OS.time, OS) ~ Age_y+Gender+tissue+LAG3+CTLA4+CD274, data = plot_df_low))

big_cox<-rbind(all_composite$coefficients,all_mid$coefficients)
big_cox<-as.data.frame(big_cox[,c(1,3,5)])
big_cox$`lower .95`<-big_cox$coef-2*big_cox$`se(coef)`
big_cox$`upper .95`<-big_cox$coef+2*big_cox$`se(coef)`

big_cox$group<-factor(c(rep("CGVHD",nrow(all_composite$coefficients)),rep("NoCGVHD",nrow(all_composite$coefficients))),levels = c("CGVHD","NoCGVHD"))
big_cox$checkpoint_group<-as.factor(paste(big_cox$group,rownames(big_cox)))

ggplot(data = big_cox, aes(x = checkpoint_group, y = `coef`, ymin = `lower .95`, ymax = `upper .95`, colour = group)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", size=1) +
  geom_errorbar(size=1.2,color="grey40",linetype=1,width=0.2) +
  geom_point(size=10,shape=15) + geom_label(label=paste("p =",round(big_cox$`Pr(>|z|)`,4)),nudge_x = 0.3) +
  coord_flip() + theme_bw(base_size = 22) + labs(y="Hazard Ratio (OS)",x="Model") + scale_color_manual(values=cols2[c(1,3)])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        , panel.background = element_blank()
        , axis.line = element_line(colour = "black"),legend.position="none",
        ,axis.text = element_text(color="black"))  + ggtitle(label="Discovery (N=110)") 
ggsave("/Users/tjsears/Code/CIBMTR/plots/Jan31/test_hazzy.pdf",width = 10,height = 14)





