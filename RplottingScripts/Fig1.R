# Fig 1 CIBMTR

library(ggplot2)
library(RColorBrewer)
library(survival)

# COHORT CHARACTERISTICS

CIBMTR_df<-read.table("/Users/tjsears/Code/CIBMTR/Tables/complete_dataset_feb25.txt",sep="\t",header=T,row.names = 1)
CIBMTR_df[is.na(CIBMTR_df)] <-0# fill all other na with 0
#CIBMTR_df = CIBMTR_df('.', 0, regex=False)
CIBMTR_df$dnrsex[CIBMTR_df$dnrsex=="."]<-1
CIBMTR_df$dnrsex<-as.numeric(CIBMTR_df$dnrsex)
#CIBMTR_df=CIBMTR_df[CIBMTR_df$condint!=99,] # keep or drop?
#CIBMTR_df=CIBMTR_df[CIBMTR_df$donorgp!=99,]

# Class I vs Class II fraction matching

# just take fraction of pts with a missmatch

plot_df<-CIBMTR_df[,c("ClassI_difference","ClassII_difference")]
plot_df$ClassI_difference<-ifelse(plot_df$ClassI_difference==0,1,0)
plot_df$ClassII_difference<-ifelse(plot_df$ClassII_difference==0,1,0)

plot_df<-data.frame(c(sum(plot_df$ClassI_difference)/length(plot_df$ClassI_difference),
                    sum(plot_df$ClassII_difference)/length(plot_df$ClassII_difference)))
plot_df$Class<-c("Class I","Class II")
colnames(plot_df)<-c("Fraction","Class")

cols<-brewer.pal(n = 3,"Dark2")

# Plotting
ggplot(plot_df, aes(x = Class, y = Fraction, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cols[1:2]) +
  labs(x = "HLA Class", y = "Fraction of Pts", fill = "No. Autoimmune Alleles") +
  theme_bw(base_size = 18) + theme(panel.grid.major = element_blank(), 
                                   panel.grid.minor = element_blank())+
  ggtitle("Fraction Matching HLA Alleles")
ggsave("/Users/tjsears/Code/CIBMTR/FigsJun/Fig1/FractionMatching.pdf",width = 6,height = 6)

# Gender
plot_df<-CIBMTR_df[,c("sex"),drop=F]
plot_df$Female<-ifelse(plot_df$sex==2,1,0)
plot_df$Male<-ifelse(plot_df$sex==1,1,0)

plot_df<-data.frame(c(sum(plot_df$Female)/length(plot_df$Female),
                      sum(plot_df$Male)/length(plot_df$Male)))
plot_df$Class<-c("Female","Male")
colnames(plot_df)<-c("Fraction","Sex")

cols<-brewer.pal(n = 3,"Dark2")

# Plotting
ggplot(plot_df, aes(x = Sex, y = Fraction, fill = Sex)) +
  geom_bar(stat = "identity", position = "stack") +
  #scale_fill_manual(values = cols[1:2]) +
  labs(x = "Patient Sex", y = "Fraction of Pts", fill = "Sex") +
  theme_bw(base_size = 18) + theme(panel.grid.major = element_blank(), 
                                   panel.grid.minor = element_blank())+
  ggtitle("Patient Sex")
ggsave("/Users/tjsears/Code/CIBMTR/Figs_Feb/Fig1/Sex.pdf",width = 6,height = 6)


# cgvhd and acute gvhd?
plot_df<-CIBMTR_df[,c("agvhd24",'agvhd34','cgvhd'),drop=F]
plot_df$AGVHD24<-ifelse(plot_df$agvhd24==1,1,0)
plot_df$AGVHD34<-ifelse(plot_df$agvhd34==1,1,0)
plot_df$CGVHD<-ifelse(plot_df$cgvhd==1,1,0)

plot_df<-data.frame(c(sum(plot_df$AGVHD24)/length(plot_df$AGVHD24),
                      sum(plot_df$AGVHD34)/length(plot_df$AGVHD34),
                      sum(plot_df$CGVHD)/length(plot_df$CGVHD)))
plot_df$Class<-c("acute GVHD 2-4","acute GVHD 3-4","chronic GVHD")
colnames(plot_df)<-c("Fraction","gvhd status")

cols<-brewer.pal(n = 3,"Paired")

# Plotting
ggplot(plot_df, aes(x = `gvhd status`, y = Fraction, fill = `gvhd status`)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cols[1:3]) +
  labs(x = "gvhd status", y = "Fraction of Pts", fill = "gvhd status") +
  theme_bw(base_size = 18) + theme(panel.grid.major = element_blank(), 
                                   panel.grid.minor = element_blank(),axis.text.x = element_blank())+
  ggtitle("Graft vs Host Disease status")
ggsave("/Users/tjsears/Code/CIBMTR/Figs_Feb/Fig1/gvhd_status.pdf",width = 6,height = 6)



# represent as fraction of pts. similar to age and gender

# survival statistics

### rel ###

# Calculate the median
median_value <- round(median(CIBMTR_df$intxrel),2)

# Create the histogram
ggplot(CIBMTR_df, aes(x = intxrel)) +
  geom_histogram(binwidth = 1, fill = "gray", color = "black") + # Adjust binwidth as needed
  geom_vline(xintercept = median_value, linetype = "dashed", color = "red") +
  annotate("label", x = median_value, y = 50, label = paste("Median =", median_value), hjust = -0.25, color = "black") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  labs(title = "Time to Relapse", x = "Months", y = "Frequency")
ggsave("/Users/tjsears/Code/CIBMTR/Figs_Feb/Fig1/intxrel.pdf",width = 6,height = 4)


### os ###

# Calculate the median
median_value <- round(median(CIBMTR_df$intxsurv),3)

# Create the histogram
ggplot(CIBMTR_df, aes(x = intxsurv)) +
  geom_histogram(binwidth = 1, fill = "gray", color = "black") + # Adjust binwidth as needed
  geom_vline(xintercept = median_value, linetype = "dashed", color = "red") +
  annotate("label", x = median_value, y = 45, label = paste("Median =", median_value), hjust = -0.25, color = "black") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  labs(title = "Time to Death", x = "Months", y = "Frequency")
ggsave("/Users/tjsears/Code/CIBMTR/Figs_Feb/Fig1/intxdead.pdf",width = 6,height = 4)


# Age 

# Calculate the median
median_value <- round(median(CIBMTR_df$Age.at.diagnosis),3)

# Create the histogram
ggplot(CIBMTR_df, aes(x = Age.at.diagnosis)) +
  geom_histogram(binwidth = 1, fill = "gray", color = "black") + # Adjust binwidth as needed
  geom_vline(xintercept = median_value, linetype = "dashed", color = "red") +
  annotate("label", x = median_value, y = 45, label = paste("Median =", median_value), hjust = -0.25, color = "black") +
  theme_bw(base_size = 18) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  labs(title = "Age at Diagnosis", x = "Years Old", y = "Frequency")
ggsave("/Users/tjsears/Code/CIBMTR/Figs_Feb/Fig1/Age.pdf",width = 6,height = 4)



# get new stats for T1 06/21/24

#agvhd24
test_df<-CIBMTR_df
test_df<-test_df[test_df$agvhd24==1,]
# report median time to 2-4 and 3-4 with confidence interval in perentheses
# Fit the Kaplan-Meier survival curve
fit <- survfit(Surv(intxagvhd24,agvhd24) ~ 1, data = test_df)

# Extract median time to event and confidence interval
median_time <- summary(fit)$table['median']
lower_ci <- summary(fit)$table['0.95LCL']
upper_ci <- summary(fit)$table['0.95UCL']

# Print the results
cat("Median time to event:", median_time, "\n")
cat("95% CI:", lower_ci, "-", upper_ci, "\n")


#agvhd34
test_df<-CIBMTR_df
test_df<-test_df[test_df$agvhd34==1,]
# report median time to 2-4 and 3-4 with confidence interval in perentheses
# Fit the Kaplan-Meier survival curve
fit <- survfit(Surv(intxagvhd34,agvhd34) ~ 1, data = test_df)

# Extract median time to event and confidence interval
median_time <- summary(fit)$table['median']
lower_ci <- summary(fit)$table['0.95LCL']
upper_ci <- summary(fit)$table['0.95UCL']

# Print the results
cat("Median time to event:", median_time, "\n")
cat("95% CI:", lower_ci, "-", upper_ci, "\n")



#cgvhd
test_df<-CIBMTR_df
test_df<-test_df[test_df$cgvhd==1,]
# report median time to 2-4 and 3-4 with confidence interval in perentheses
# Fit the Kaplan-Meier survival curve
fit <- survfit(Surv(intxcgvhd,cgvhd) ~ 1, data = test_df)

# Extract median time to event and confidence interval
median_time <- summary(fit)$table['median']
lower_ci <- summary(fit)$table['0.95LCL']
upper_ci <- summary(fit)$table['0.95UCL']

# Print the results
cat("Median time to event:", median_time, "\n")
cat("95% CI:", lower_ci, "-", upper_ci, "\n")
