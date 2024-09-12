# oncoprint for Fig 2 CIBMTR

library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(reshape2)

remove_Whitespace=T
# read in clinical data
clin_data<-read.table("/Users/tjsears/Code/CIBMTR/Tables/complete_dataset_feb25.txt",sep="\t",header=T)

# import phbr data, collapse into gene by mut type grid
phbr<-read.table("/Users/tjsears/Code/CIBMTR/Tables/donor_phbr.txt",sep="\t",header=T)

phbr$Gene<-sapply(strsplit(phbr$Mutation_ID, "_"), function(x) x[2])
phbr$Mut_Pt_Id<-paste(phbr$Gene,phbr$Individual)
phbr$AF<-phbr$t_alt_count/phbr$t_depth
phbr<-phbr[order(phbr$AF,decreasing = T),]
phbr<-phbr[!duplicated(phbr$Mut_Pt_Id,),]

phbr<-merge(phbr,clin_data,by="Individual",all=T)

matrix_df <- dcast(phbr, Gene ~ Individual, value.var = "Variant_Classification")
matrix_df[is.na(matrix_df)]<-"no mutation"
rownames(matrix_df)<-matrix_df$Gene
matrix_df<-matrix_df[,-1]

# order rows of matrix by rowsum
gene_list<-data.frame(table(phbr$Gene))[order(data.frame(table(phbr$Gene))[,2],decreasing = T),1]
ID_list<-data.frame(table(phbr$Individual))[order(data.frame(table(phbr$Individual))[,2],decreasing = T),1]

matrix_df<-matrix_df[gene_list,]

# take only top 50 genes
matrix_df<-matrix_df[1:30,]
gene_list<-gene_list[1:30]

numeric_matrix <- matrix_df
numeric_matrix[numeric_matrix=="no mutation"]<-0
numeric_matrix[numeric_matrix=="Missense_Mutation"]<-1
numeric_matrix[numeric_matrix=="In_Frame_Ins"]<-1
numeric_matrix[numeric_matrix=="In_Frame_Del"]<-1
numeric_matrix[numeric_matrix=="Frame_Shift_Del"]<-1
numeric_matrix[numeric_matrix=="Frame_Shift_Ins"]<-1
numeric_matrix <- apply(numeric_matrix, 2, as.numeric)
ID_list<-colnames(numeric_matrix)[order(colSums(numeric_matrix),decreasing = T)]
  
# have a toggle where we remove whitespace
if (remove_Whitespace==T){
numeric_matrix<-numeric_matrix[,colSums(numeric_matrix)>0]
ID_list<-colnames(numeric_matrix)[order(colSums(numeric_matrix),decreasing = T)]
}

matrix_df<-matrix_df[,ID_list]

# Clinical annotations
clinical_data <- data.frame(
  donor_HLA_mismatch = ifelse(clin_data$ClassI_difference>0|clin_data$ClassII_difference>0,"Yes","No"),
  cgvhd = ifelse(clin_data$cgvhd>0,"Yes","No"),
  agvhd = ifelse(clin_data$agvhd34>0,"Grade 3-4",ifelse(clin_data$agvhd24>0,"Grade 2","None")),
  Individual = clin_data$Individual,
  AI_I = ifelse(clin_data$AI.I>0,"class-I autoimmune allele (+)","class-I autoimmune allele (-)"),
  AI_II = ifelse(clin_data$AI.II>0,"class-II autoimmune allele (+)","class-II autoimmune allele (-)"),
  PHBR_I = clin_data$Poor..High..PHBR.I,
  PHBR_II = clin_data$Poor..High..PHBR.II
)

rownames(clinical_data) <- clinical_data$Individual
clinical_data[is.na(clinical_data)]<-0
clinical_data$agvhd[clinical_data$agvhd==0]<-'None'

clinical_data$donor_HLA_mismatch[is.na(clinical_data$donor_HLA_mismatch)] <- "No"
clinical_data$cgvhd[is.na(clinical_data$cgvhd)] <- "No"
clinical_data$agvhd[is.na(clinical_data$agvhd)] <- "None"
clinical_data$AI_I[is.na(clinical_data$AI_I)] <- "class-I autoimmune allele (-)"
clinical_data$AI_II[is.na(clinical_data$AI_II)] <- "class-II autoimmune allele (-)"

clinical_data <- clinical_data[clinical_data$Individual %in% colnames(matrix_df),]
rownames(clinical_data) <- clinical_data$Individual
clinical_data <- clinical_data[colnames(matrix_df),]
clinical_data <- clinical_data[,!colnames(clinical_data) %in% "Individual"]

# Perform T-tests comparing the frequency of mutations across the AI.I condition
gene_list <- rownames(matrix_df)

t_test_results_ai_i <- data.frame(Gene = gene_list, P_Value = NA)
t_test_results_ai_ii <- data.frame(Gene = gene_list, P_Value = NA)
t_test_results_ai <- data.frame(Gene = gene_list, P_Value = NA)
t_test_results_agvhd2 <- data.frame(Gene = gene_list, P_Value = NA)
t_test_results_agvhd34 <- data.frame(Gene = gene_list, P_Value = NA)
t_test_results_cgvhd <- data.frame(Gene = gene_list, P_Value = NA)
t_test_results_donormismatch <- data.frame(Gene = gene_list, P_Value = NA)
t_test_results_PHBRI <- data.frame(Gene = gene_list, P_Value = NA)
t_test_results_PHBRII <- data.frame(Gene = gene_list, P_Value = NA)

temp_matrix<-ifelse(matrix_df=='no mutation',0,1)

# correct pvals
# store t-stat too

for (gene in gene_list) {
  mutation_freq <- as.numeric(temp_matrix[gene, ])
  
  ai_i_groups <- clinical_data$AI_I
  t_test <- t.test(mutation_freq ~ ai_i_groups)
  t_test_results_ai_i[t_test_results_ai_i$Gene == gene, "P_Value"] <- t_test$p.value
  t_test_results_ai_i[t_test_results_ai_i$Gene == gene, "T_Stat"] <- t_test$statistic
  
  ai_ii_groups <- clinical_data$AI_II
  t_test <- t.test(mutation_freq ~ ai_ii_groups)
  t_test_results_ai_ii[t_test_results_ai_ii$Gene == gene, "P_Value"] <- t_test$p.value
  t_test_results_ai_ii[t_test_results_ai_ii$Gene == gene, "T_Stat"] <- t_test$statistic
  
  ai_ii_groups <- ifelse(clinical_data$AI_I=='class-I autoimmune allele (+)'|clinical_data$AI_II=='class-II autoimmune allele (+)',1,0)
  t_test <- t.test(mutation_freq ~ ai_ii_groups)
  t_test_results_ai[t_test_results_ai$Gene == gene, "P_Value"] <- t_test$p.value
  t_test_results_ai[t_test_results_ai$Gene == gene, "T_Stat"] <- t_test$statistic
  
  ai_i_groups <- ifelse(clinical_data$agvhd=='Grade 2',1,0)
  t_test <- t.test(mutation_freq ~ ai_i_groups)
  t_test_results_agvhd2[t_test_results_agvhd2$Gene == gene, "P_Value"] <- t_test$p.value
  t_test_results_agvhd2[t_test_results_agvhd2$Gene == gene, "T_Stat"] <- t_test$statistic
  
  ai_i_groups <- ifelse(clinical_data$agvhd=='Grade 3-4',1,0)
  t_test <- t.test(mutation_freq ~ ai_i_groups)
  t_test_results_agvhd34[t_test_results_agvhd34$Gene == gene, "P_Value"] <- t_test$p.value
  t_test_results_agvhd34[t_test_results_agvhd34$Gene == gene, "T_Stat"] <- t_test$statistic
  
  ai_i_groups <- clinical_data$cgvhd
  t_test <- t.test(mutation_freq ~ ai_i_groups)
  t_test_results_cgvhd[t_test_results_cgvhd$Gene == gene, "P_Value"] <- t_test$p.value
  t_test_results_cgvhd[t_test_results_cgvhd$Gene == gene, "T_Stat"] <- t_test$statistic
  
  ai_i_groups <- clinical_data$PHBR_I
  t_test <- glm(mutation_freq ~ ai_i_groups)
  t_test_results_PHBRI[t_test_results_PHBRI$Gene == gene, "P_Value"] <- summary(t_test)$coefficient[2,4]
  t_test_results_PHBRI[t_test_results_PHBRI$Gene == gene, "T_Stat"] <- summary(t_test)$coefficient[2,3]
  
  ai_i_groups <- clinical_data$PHBR_II
  t_test <- glm(mutation_freq ~ ai_i_groups)
  t_test_results_PHBRII[t_test_results_PHBRII$Gene == gene, "P_Value"] <- summary(t_test)$coefficient[2,4]
  t_test_results_PHBRII[t_test_results_PHBRII$Gene == gene, "T_Stat"] <- summary(t_test)$coefficient[2,3]
}

# Average allele depth for each gene (row)
avg_allele_depth <- aggregate(phbr$AF, by = list(phbr$Gene), FUN = mean)
rownames(avg_allele_depth) <- avg_allele_depth$Group.1
avg_allele_depth <- avg_allele_depth[gene_list,]
avg_allele_depth <- avg_allele_depth[,2]

# Get % mutated
mutation_percent <- as.data.frame(table(phbr$Gene))
rownames(mutation_percent) <- mutation_percent$Var1
mutation_percent <- mutation_percent[gene_list,]
mutation_percent <- round(mutation_percent[,2] / (494 / 100), 1) # get % of cohort mutated
#mutation_percent <- paste("(", mutation_percent, ")%", sep = "")
#mutation_percent[1] <- paste(mutation_percent[1], "Percent Mutated")

# Rename mutations to look better
matrix_df[] <- lapply(matrix_df, function(x) gsub("_", " ", x))
matrix_df[is.na(matrix_df)] <- 0
# Create the oncoprint
matrix_df[matrix_df == 'In Frame Ins'] <- "In Frame Ins/Del"
matrix_df[matrix_df == 'In Frame Del'] <- "In Frame Ins/Del"

pdf(file="/Users/tjsears/Code/CIBMTR/FigsJun/FigS4/Oncoprint.pdf",width = 12,height = 6)

clinical_data<-clinical_data[,c(1:4)]
ht_list = Heatmap(as.matrix(matrix_df),
                  col = c("no mutation" = "white", 
                          "Missense Mutation" = "blue", 
                          "In Frame Ins/Del" = "#3BB273",
                          "Frame Shift Ins" = "#E1BC29",
                          "Frame Shift Del" = "black"),
                  top_annotation = HeatmapAnnotation(df = clinical_data, simple_anno_size = unit(0.3, "cm"),
                                                     col = list(donor_HLA_mismatch = c("Yes" = "#5CD3FF", "No" = "grey90"), 
                                                                cgvhd = c("Yes" = "#F2C078", "No" = "grey90"), 
                                                                agvhd = c("Grade 3-4" = "gold", "Grade 2" = "darkgrey", "None" = "grey90")),
                                                     gp = gpar(col = "black", lwd = 0.3), border = F),
                  right_annotation = rowAnnotation(PercentMutated = anno_barplot(mutation_percent,rot = 45,
                                                                                 border = TRUE, gp = gpar(fill = "#FF000080"))),
                  cluster_rows = F, show_column_names = F, column_split = clinical_data$AI_I,
                  cluster_columns = F, border_gp = gpar(col = "black", lwd = 0.3), rect_gp = gpar(col = "grey80", lwd = 0.1),
                  row_names_gp = gpar(cex = 0.7), row_names_side = "left", show_heatmap_legend = F
)

# Draw the heatmap
draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

##### matrix of statistical comparisons ##### 

# compare top mutations for AI, AI.I, AI.II, agvhd2, agvhd2-4, cgvhd

pval_mat<-cbind(t_test_results_ai_i[,c(2)],t_test_results_ai_ii[,2],t_test_results_ai[,2],
                t_test_results_agvhd2[,2],t_test_results_agvhd34[,2],
                t_test_results_cgvhd[,2])

rownames(pval_mat)<-t_test_results_ai_i$Gene
colnames(pval_mat)<-c('class-I autoimmune allele','class-II autoimmune allele','AI','agvhd2','agvhd34','cgvhd')

pval_mat_adj<-apply(pval_mat,2,function(x) p.adjust(x,method = "BH"))

t_mat<-cbind(t_test_results_ai_i[,c(3)],t_test_results_ai_ii[,3],t_test_results_ai[,3],
                t_test_results_agvhd2[,3],t_test_results_agvhd34[,3],
                t_test_results_cgvhd[,3])

rownames(t_mat)<-t_test_results_ai_i$Gene
colnames(t_mat)<-c('class-I autoimmune allele','class-II autoimmune allele','AI','agvhd2','agvhd34','cgvhd')

library(RColorBrewer)
rev.color <- colorRampPalette(c("slateblue", "whitesmoke", "firebrick1"))(200)

# compare TMB with AI, AI.I, AI.II, agvhd2, agvhd2-4, cgvhd

phbr_stat_test<-phbr[!is.na(phbr$t_depth),]
phbr_stat_df<-as.data.frame(table(phbr_stat_test$Individual))

# merge with clinical data

clin_TMB_test<-merge(clin_data,phbr_stat_df,by.x='DID',by.y='Var1',all.x = T)
clin_TMB_test$Freq[is.na(clin_TMB_test$Freq)]<-0
clin_TMB_test$TMB<-clin_TMB_test$Freq
clin_TMB_test$AI<-as.factor(clin_TMB_test$AI)
clin_TMB_test$AI.I<-as.factor(clin_TMB_test$AI.I)
clin_TMB_test$AI.II<-as.factor(clin_TMB_test$AI.II)
clin_TMB_test$cgvhd<-as.factor(clin_TMB_test$cgvhd)
clin_TMB_test$agvhd2<-ifelse(clin_TMB_test$agvhd24!=0&clin_TMB_test$agvhd34==0,1,0)
clin_TMB_test$agvhd2<-as.factor(clin_TMB_test$agvhd2)
clin_TMB_test$agvhd34<-as.factor(clin_TMB_test$agvhd34)
clin_TMB_test$HLA_mismatch<-ifelse(clin_TMB_test$ClassI_difference>0|clin_TMB_test$ClassII_difference>0,"Yes","No")

clin_TMB_test$PHBRI<-clin_TMB_test$Poor..High..PHBR.I
clin_TMB_test$PHBRI[is.na(clin_TMB_test$PHBRI)]<-0
clin_TMB_test$PHBRII<-clin_TMB_test$Poor..High..PHBR.II
clin_TMB_test$PHBRII[is.na(clin_TMB_test$PHBRII)]<-0

TMB_t<-c()
TMB_p<-c()

for (clinvar in c('AI.I','AI.II','AI','agvhd2','agvhd34','cgvhd','PHBRI','PHBRII')){
  
  if (clinvar %in% c('PHBRI','PHBRII')) {
    
    TMB<-clin_TMB_test$TMB
    ai_i_groups <- clin_TMB_test[,clinvar]
    test<-glm(TMB~ai_i_groups)
    TMB_t<-c(TMB_t,summary(test)$coefficient[2,3])
    TMB_p<-c(TMB_p,summary(test)$coefficient[2,4])
    
  } else {
    a<-clin_TMB_test[clin_TMB_test[,clinvar]==0,'TMB']
    b<-clin_TMB_test[clin_TMB_test[,clinvar]==1,'TMB']
    test<-t.test(a,b)
    TMB_t<-c(TMB_t,test$statistic)
    TMB_p<-c(TMB_p,test$p.value)
  }
  
}

pval_mat_plot<-rbind(pval_mat_adj,TMB_p)
t_mat_plot<-rbind(t_mat,TMB_t)
rownames(t_mat_plot)[nrow(t_mat_plot)]<-'TMB'
rownames(pval_mat_plot)[nrow(pval_mat_plot)]<-'TMB'

library(ggplot2)
library(corrplot)

pdf(file="/Users/tjsears/Code/CIBMTR/FigsJun/FigS4/SomaticClinicalPlot.pdf",width = 10,height = 6)

#colnames(t_mat_plot)[c(ncol(t_mat_plot)-1,ncol(t_mat_plot))]<-c('Improved PHBR-I','Improved PHBR-II')
colnames(t_mat_plot)[3]<-'any autoimmune allele'
colnames(pval_mat_plot)[3]<-'any autoimmune allele'

corrplot(
  t(t_mat_plot),diag=T,
  method = "color",is.corr = FALSE,
  type = "full",col = rev.color,#col.lim = c(-0.35,0.25),
  #order = "hclust",
  tl.col = "black",mar=c(0,0,2,2),cl.ratio=0.4,cl.cex=0.7,
  tl.srt = 45,tl.cex = 0.7,cl.pos = 'b',cl.offset = 0.5,
  sig.level = c(0.001,0.01,0.05), # Set the significance level
  insig = "label_sig",outline=T,title = "Clinical x somatic associations",
  p.mat = t(pval_mat_plot),# Display p-values for all cells
  pch.cex = 1.5 # Adjust the size of p-value labels
  #addCoef.col = "black" # Set the color of p-value labels
)
dev.off()
