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
CIBMTR_df$WorsePHBR.I<-CIBMTR_df$PHBR.I.change
CIBMTR_df$WorsePHBR.II<-CIBMTR_df$PHBR.II.change
CIBMTR_df$ImprovedPHBR.I<--1*CIBMTR_df$Poor..High..PHBR.I
CIBMTR_df$ImprovedPHBR.II<--1*CIBMTR_df$Poor..High..PHBR.II

CIBMTR_df$CGVHD<-ifelse(CIBMTR_df$CGVHD>0,"Yes","No")
CIBMTR_df$AI.I<-ifelse(CIBMTR_df$AI.I>0,"Yes","No")

#library(Rcmdr)

#library(survival, pos=15)
#library(aod, pos=16)

rerun_simon=F
if (rerun_simon==T){
  
  
  TempDF <- CIBMTR_df
  TempDF <- TempDF[TempDF$AI.I=='No',]
  TempTD <- stsplit(TempDF, TempDF$intxrel, TempDF$rel, 
                    TempDF$intxcgvhd, TempDF$cgvhd, TempDF$intxrel)
  CoxModel.1 <- coxph(Surv(start_td, stop_td, endpoint_td==1) ~ covariate_td, 
                      data=TempTD, method="breslow")
  res <- NULL
  (res <- summary(CoxModel.1))
  cox.table <- NULL
  cox.table <- signif(cbind(t(res$conf.int[,c(1,3,4)]), 
                            p.value=res$coefficients[,5]), digits=4)
  rownames(cox.table) <- rownames(res$coefficients)
  colnames(cox.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", 
                                                              "Lower 95%CI", "Upper 95%CI", "p.value"))
  cox.table
  waldtest(CoxModel.1)
  print(cox.zph(CoxModel.1))
  
  
  # YOU NEED TO ACTIVATE THE PLUGIN ONCE RCMDR OPENS
  
  #####Simon-Makuch plot#####
  Mantel.Byar(Group = "covariate_td", Event = TempTD$endpoint_td, 
              StartTime = TempTD$start_td, StopTime = TempTD$stop_td, 
              method = c("Tominaga"), plot=1, landmark=0)
  
  
  TempDF <- CIBMTR_df
  TempDF <- TempDF[TempDF$AI.I=='Yes',]
  TempTD <- stsplit(TempDF, TempDF$intxrel, TempDF$rel, 
                    TempDF$intxcgvhd, TempDF$cgvhd, TempDF$intxrel)
  CoxModel.1 <- coxph(Surv(start_td, stop_td, endpoint_td==1) ~ covariate_td, 
                      data=TempTD, method="breslow")
  res <- NULL
  (res <- summary(CoxModel.1))
  cox.table <- NULL
  cox.table <- signif(cbind(t(res$conf.int[,c(1,3,4)]), 
                            p.value=res$coefficients[,5]), digits=4)
  rownames(cox.table) <- rownames(res$coefficients)
  colnames(cox.table) <- gettext(domain="R-RcmdrPlugin.EZR",c("Hazard ratio", 
                                                              "Lower 95%CI", "Upper 95%CI", "p.value"))
  cox.table
  waldtest(CoxModel.1)
  print(cox.zph(CoxModel.1))
  
  #####Simon-Makuch plot#####
  Mantel.Byar(Group = "covariate_td", Event = TempTD$endpoint_td, 
              StartTime = TempTD$start_td, StopTime = TempTD$stop_td, 
              method = c("Tominaga"), plot=1, landmark=0)
  
  
  # Now, we are going to take the outputs from the simon tests, copy paste them into tables
  AI_0<-read.table('/Users/tjsears/Code/CIBMTR/Tables/SimonCurve_AI0.txt',sep='\t',header=T)
  AI_0$CGVHD<-as.factor(AI_0$CGVHD)
  
  # re-load those tables into a new dataframe
  
  # then manually plot these as better looking curves
  
  # then finally we will make a hazard ratio comparison test betweent he two settings?
  ggplot(AI_0, aes(x = time, y = survival, group = CGVHD, color = CGVHD)) +
    geom_step() + # Use geom_step for stepwise curves
    labs(title = "Kaplan-Meier Curve",
         x = "Time",
         y = "Survival Probability") +
    theme_minimal() +
    scale_color_manual(values = c(`0` = "#00BFC4", `1` = "#F8766D"))
  
  
  
  # Now, we are going to take the outputs from the simon tests, copy paste them into tables
  AI_1<-read.table('/Users/tjsears/Code/CIBMTR/Tables/SimonCurve_AI1.txt',sep='\t',header=T)
  AI_1$CGVHD<-as.factor(AI_1$CGVHD)
  
  # re-load those tables into a new dataframe
  
  # then manually plot these as better looking curves
  
  # then finally we will make a hazard ratio comparison test betweent he two settings?
  ggplot(AI_1, aes(x = time, y = survival, group = CGVHD, color = CGVHD)) +
    geom_step() + # Use geom_step for stepwise curves
    labs(title = "Kaplan-Meier Curve",
         x = "Time",
         y = "Survival Probability") +
    theme_minimal() +
    scale_color_manual(values = c(`0` = "#00BFC4", `1` = "#F8766D"))
}

# compare and plot hazzies


# Example data
logHR1 <- log(0.5744)  # log of HR for group 1
logHR2 <- log(0.289)  # log of HR for group 2
SE1 <- 0.2567         # SE of log(HR) for group 1
SE2 <- 0.2011          # SE of log(HR) for group 2

# Calculate variances
var1 <- SE1^2
var2 <- SE2^2

# Calculate the Z-score
Z <- (logHR1 - logHR2) / sqrt(var1 + var2)

# Calculate p-value
p_value <- 2 * pnorm(-abs(Z))  # Two-tailed test

# Output the results
cat("Z-score:", Z, "\n")
cat("P-value:", p_value, "\n")

# data
# se 0.2567

#Hazard ratio Lower 95%CI Upper 95%CI p.value
#covariate_td       0.5744      0.3473        0.95 0.03078

# se 0.2011

#Hazard ratio Lower 95%CI Upper 95%CI   p.value
#covariate_td        0.289      0.1948      0.4286 6.762e-10


# plot hazard ratios
big_cox<-data.frame(AI_group=c("cGVHD\nClass-I autoimmune allele (-)",'cGVHD\nClass-I autoimmune allele (+)'),`Hazard ratio`=c(0.5744,0.289),
        `lower .95`=c(0.3473,0.1948),`upper .95`=c(0.95,0.4286),p.value=c(0.03078,6.762e-10))

big_cox <- big_cox %>%
  mutate(
    label_text = ifelse(p.value < 0.05, 
                        sprintf("bold('p = %g')", p.value),
                        sprintf("'p = %g'", p.value)
    )
  )

#big_cox$Hazard.ratio<-log(big_cox$Hazard.ratio)
#big_cox$lower..95<-log(big_cox$lower..95)
#big_cox$upper..95<-log(big_cox$upper..95)

pval_df<-data.frame(group1 = c("cGVHD\nClass-I autoimmune allele (-)"),
group2 = c("cGVHD\nClass-I autoimmune allele (+)"),
p.value = c("*")) #0.03516462

ggplot(data = big_cox, aes(x = AI_group, y = `Hazard.ratio`, ymin = `lower..95`, ymax = `upper..95`, colour = AI_group)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black", size=1) + ylim(c(0.1,1.25)) +
  geom_signif(y_position = c(1.075), comparisons=list(c("cGVHD\nClass-I autoimmune allele (-)","cGVHD\nClass-I autoimmune allele (+)")),vjust=0.5,
                     annotations = pval_df$p.value, map_signif_level = TRUE,size = 1,textsize = 12,colour = 'black',) +
  geom_errorbar(aes(color = "green"), size=1.2, linetype=1, width=0.2) + # Apply color directly in geom_errorbar
  geom_point(aes(color = "black"), size=5, shape=15) + # Apply color directly in geom_point
  geom_label(aes(label = label_text,color='black'), parse = TRUE, nudge_x = c(0.35,-0.4), nudge_y = 0) +  
  coord_flip() +
  theme_bw(base_size = 22) + 
  labs(y = "Hazard Ratio (Relapse-Free Survival)", x = NULL) +
  scale_color_manual(values = c("green" = "#0C7C59",black='gray20')) + # Ensure this matches with assigned colors
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black")) + ggtitle(label="cGVHD and Class-I Autoimmune Alleles")

ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig3/simon_hazzy.pdf",width = 10.5,height = 3.2)


##############################
# replotting simon km curves #
##############################

# using the HR and pval info from the simon plots, we can recreate these km curves properly

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

cols<-brewer.pal(n = 4,name = "Dark2")[3:4]

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

CIBMTR_df$cGVHD<-ifelse(CIBMTR_df$CGVHD>0,"Yes","No")
CIBMTR_df$AI.I<-ifelse(CIBMTR_df$AI.I>0,"Yes","No")

CIBMTR_df#<-CIBMTR_df[CIBMTR_df$intxsurv>=270,]
CIBMTR_df_0<-CIBMTR_df[CIBMTR_df$AI.I=='No',]
CIBMTR_df_1<-CIBMTR_df[CIBMTR_df$AI.I=='Yes',]

quantile(CIBMTR_df$intxcgvhd)

CIBMTR_df$ImprovedPHBR.II_bin<-ifelse(CIBMTR_df$ImprovedPHBR.II>0,1,CIBMTR_df$ImprovedPHBR.II)
CIBMTR_df$ImprovedPHBR.II_bin<-ifelse(CIBMTR_df$ImprovedPHBR.II<0,2,CIBMTR_df$ImprovedPHBR.II_bin)

km_trt_fit <- survfit(Surv(intxrel, rel) ~ cGVHD, data=CIBMTR_df_1)
res<- pairwise_survdiff(Surv(intxrel, rel) ~ cGVHD, data=CIBMTR_df_1,p.adjust.method = "none",rho=0)#, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

# annotate HR and pval
# pval = 1.052e-10 AI 1 simon curve
# pval = 0.03006 AI 0 simon curve
p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df_1,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),
               pval=F,xlim=c(0,2000),break.x.by=500,ncensor.plot.height=0.25,#risk.table.title="",
               font.x=18,font.y=18,font.legend=13,font.title=18,#risk.table = c("absolute"),
               palette=cols,legend.title="",title="cGVHD and Class-I Autoimmune Allele (+)\nRelapse-Free Survival"
               ,xlab="Time in Days",ylim=c(0,1),ylab="Relapse-Free Survival (%)",surv.scale = c("percent"),
               risk.table.title="Number at Risk")
p1$plot+ggplot2::annotate("label",x=1200,y = 0.9,size=6, label = paste("P =",1.052e-10))

ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig3/Disc_AI_REL_km_split_AI_1.pdf",width = 7,height = 6)

km_trt_fit <- survfit(Surv(intxrel, rel) ~ cGVHD, data=CIBMTR_df_0)
res<- pairwise_survdiff(Surv(intxrel, rel) ~ cGVHD, data=CIBMTR_df_0,p.adjust.method = "none",rho=0)#, p.adjust.method = "none", rho = 0)
res$p.value<-round(res$p.value, digits = 4)
res

p1<-ggsurvplot(km_trt_fit,data=CIBMTR_df_0,size=1.2,censor.shape="|", censor.size = 6,legend=c("top"),
               pval=F,xlim=c(0,2000),break.x.by=500,ncensor.plot.height=0.25,#risk.table.title="",
               font.x=18,font.y=18,font.legend=13,font.title=18,#risk.table = c("absolute"),
               palette=cols,legend.title="",title="cGVHD and Class-I Autoimmune Allele (-)\nRelapse-Free Survival"
               ,xlab="Time in Days",ylim=c(0,1),ylab="Relapse-Free Survival (%)",surv.scale = c("percent"),
               risk.table.title="Number at Risk",linetype = 'D3')
p1$plot+ggplot2::annotate("label",x=1200,y = 0.9,size=6, label = paste("P =",round(0.03006,4)))

ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig3/Disc_AI_REL_km_split_AI_0.pdf",width = 7,height = 6)






