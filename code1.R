rm(list = ls())
library(Maaslin2)
library(MicrobiotaProcess)
library(tidyverse)
library(compositions)
library(pheatmap)
library(phyloseq)
library(ggalluvial)
library(dplyr)
library(qvalue)
library(gridExtra)
library(ggsignif)
library(ggpubr)
library(data.table)
library(ggplot2)
library(ggExtra)
library(ggtext)
library(vegan)
library(caret)
library(data.table)
library(ggtext)
library(ggbreak)
meta_D<-read.table('D:/microbiome/micom/202404/paper/data/meta_Depression.txt',sep = '\t',header = T)
meta_H<-read.table('D:/microbiome/micom/202404/paper/data/meta_healthy.txt',sep = '\t',header = T)
colnames(meta_H)[8]<-'disease'
colnames(meta_D)[8]<-'disease'
meta_H$disease<-'Healthy'
meta_D$disease<-'Depression'
age<-rbind(meta_H[,c(6,8),drop=F],meta_D[,c(6,8),drop=F])
p1<-ggplot(age, aes(x = Host.age, fill = disease)) +
  geom_density(alpha = 0.5) + 
  labs(x = "Age(Before matching)") +
  theme_minimal() +
  scale_fill_manual(values=c('Depression'='#E7B800', 'Healthy'='#00AFBB'))+
  ggtitle('C')+
  theme(legend.position = 'none')+
  theme(plot.title = element_text(face = 'bold'))+
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())
p1##p.val:6.265e-15
BMI<-rbind(meta_H[,c(7,8),drop=F],meta_D[,c(7,8),drop=F])
p2<-ggplot(BMI, aes(x = BMI, fill = disease)) +
  geom_density(alpha = 0.5) +  
  labs(x = "BMI(Before matching)") +
  theme_minimal() +
  scale_fill_manual(values=c('Depression'='#E7B800', 'Healthy'='#00AFBB'))+
  xlim(c(0,55))+
  ggtitle('E')+
  theme(plot.title = element_text(face = 'bold'))+
  theme(legend.position = 'none')+
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())
p2##p.val:0.4529
meta_H<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/meta_H_used.csv',row.names = 1)
meta_D<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/meta_D_used.csv',row.names = 1)
colnames(meta_H)[5]<-'disease'
colnames(meta_D)[5]<-'disease'
meta_H$disease<-'Healthy'
meta_D$disease<-'Depression'
age<-rbind(meta_H[,c(3,5),drop=F],meta_D[,c(3,5),drop=F])
p3<-ggplot(age, aes(x = Host.age, fill = disease)) +
  geom_density(alpha = 0.5) + 
  labs(x = "Age(After matching)") +
  theme_minimal() +
  scale_fill_manual(values=c('Depression'='#E7B800', 'Healthy'='#00AFBB'))+
  ggtitle('D')+
  theme(plot.title = element_text(face = 'bold'))+
  theme(legend.position = 'none')+
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())
p3##p.val:0.5574
BMI<-rbind(meta_H[,c(4,5),drop=F],meta_D[,c(4,5),drop=F])
p4<-ggplot(BMI, aes(x = BMI, fill = disease)) +
  geom_density(alpha = 0.5) +  
  labs(x = "BMI(After matching)") +
  theme_minimal() +
  scale_fill_manual(values=c('Depression'='#E7B800', 'Healthy'='#00AFBB'))+
  ggtitle('F')+
  theme(plot.title = element_text(face = 'bold'))+
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())
p4##p.val:0.8053
pdf(file = 'D:/microbiome/micom/202404/paper/figures/fig1.pdf',width = 15,height = 5.82)
grid.arrange(p1, p3,p2,p4,nrow = 1)
dev.off()

abundance_D<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/depression_genus.csv',row.names = 1)
abundance_D<-abundance_D[match(rownames(meta_D),colnames(abundance_D))]
temp<-abundance_D
temp[temp!=0]<-1
abundance_D<-abundance_D[which(rowSums(temp)>ncol(abundance_D)*0.1)%>%as.numeric(),]
abundance_H<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/healthy_genus.csv',row.names = 1)
abundance_H<-abundance_H[match(rownames(meta_H),colnames(abundance_H))]
temp<-abundance_H
temp[temp!=0]<-1
#abundance_H<-abundance_H[which(rowSums(temp)>ncol(abundance_H)*0.1)%>%as.numeric(),]
abundance<-merge(abundance_D,abundance_H,by = "row.names",all = TRUE)
abundance[is.na(abundance)]<-0
rownames(abundance)<-abundance$Row.names
abundance$Row.names<-NULL
#abundance_clr<-t(abundance/100)%>%clr()%>%as.matrix()%>%t()
meta<-rbind(meta_D,meta_H)
files<-list.files('D:/microbiome/micom/genus_models/genus_models/')
for(i in 1:length(files)){
  files[i]<-substr(files[i],1,nchar(files[i])-4)
}
abundance1<-abundance_H[files,]
abundance1<-na.omit(abundance1)
colSums(abundance1)%>%hist()
##########alpha
OTU<-phyloseq::otu_table(abundance,taxa_are_rows = T)
meta1<-phyloseq::sample_data(meta)
physeq<-phyloseq(OTU,meta1)
mpse<-as.MPSE(physeq)
mpse@assays@data@listData[["Abundance"]]<-mpse@assays@data@listData[["Abundance"]]*1e10
mpse %<>% 
  mp_decostand(.abundance=Abundance)
cols<-c('Depression'='#E7B800', 'Healthy'='#00AFBB')
mpse %<>% mp_cal_alpha(.abundance = Abundance,force = T)
p5 <- mpse %>% 
  mp_plot_alpha(
    .alpha = c(Observe, Shannon,Chao1,Simpson, Pielou),
    .group = disease,
  ) +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) +
  theme(
    legend.position="none",
    strip.background = element_rect(colour=NA, fill="grey")
  )+
  ggtitle('A')+
  theme(plot.title = element_text(face = 'bold'))
p5
##########anosim
mpse %>% 
  mp_anosim(.abundance=Abundance, .group=disease, action="get")
mpse %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")
p6 <- mpse %>% mp_plot_dist(.distmethod = bray, .group = disease, group.test=TRUE, textsize=2)+
  ggtitle('C')+
  theme(plot.title = element_text(face = 'bold'))
p6 
##########pcoa
dist_bray <- phyloseq::distance(physeq, method = "bray")
pcoa <- ordinate(physeq, method = "PCoA", distance = dist_bray)
pcoa_df <- data.frame(pcoa$vectors)
pcoa_df$SampleID <- rownames(pcoa_df)
sample_data<-cbind(rownames(meta),meta$disease)%>%as.data.frame()
colnames(sample_data)<-c('SampleID','disease')
rownames(sample_data)<-sample_data$SampleID
pcoa_df <- merge(pcoa_df, sample_data, by.x = "SampleID", by.y = "SampleID")
p7<-ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = disease)) +
  geom_point(alpha = 0.4) +
  labs(
    x = paste0("PCOA1 (", round(pcoa[["values"]][["Eigenvalues"]][1], 2), "%)"),
    y = paste0("PCOA2 (", round(pcoa[["values"]][["Eigenvalues"]][2], 2), "%)"),
    color = "label"
  ) +
  theme_minimal()+
  theme(legend.position = 'none')+
  stat_ellipse(aes(color = disease),type = "t", level = 0.95)+
  scale_color_manual(values=c('Depression'='#E7B800', 'Healthy'='#00AFBB'))+
  annotate("text", x = -0.3, y = -0.3, label = "ANOSIM R: 0.007962 
           P: 0.004")+
  ggtitle('B')+
  theme(plot.title = element_text(face = 'bold'))
p7
p7<-ggMarginal(p7, type = "boxplot", margins = "both", size = 5, groupColour = TRUE, groupFill = TRUE)
p7
#abundance1<-t(abundance)%>%as.data.frame()
#env<-cbind(meta$Host.age,meta$BMI)%>%as.data.frame()
#colnames(env)<-c('Age','BMI')
#rownames(env)<-rownames(abundance1)
#cca <- cca(abundance1 ~ ., data = env)
#cca_summary<-summary(cca)
#cca_df<-cca_summary$sites%>%as.data.frame()
#cca_df$disease<-meta$disease[match(rownames(cca_df),rownames(meta))]
#scores <- as.data.frame(scores(cca)$biplot)
#p7<-ggplot(cca_df, aes(x = CCA1, y = CCA2, color = disease)) +
#  geom_point(alpha = 0.4) +
#  geom_segment(data = scores, aes(x = 0, xend = CCA1*20, y = 0, yend = CCA2*20), 
#               arrow = arrow(length = unit(0.2, "cm")), color = "red")+
#  labs(
#    x = paste0("CCA1 (", round(cca$CCA$eig[1] * 100, 2), "%)"),
#    y = paste0("CCA2 (", round(cca$CCA$eig[2] * 100, 2), "%)"),
#    color = "label"
#  ) +
#  theme_minimal()+
#  theme(legend.position = 'none')+
#  stat_ellipse(aes(color = disease),type = "t", level = 0.95)+
#  scale_color_manual(values=c('Depression'='#E7B800', 'Healthy'='#00AFBB'))+
#  annotate("text", x = -10, y = -20, label = "ANOSIM R: 0.007962 
#           P: 0.004")+
#  ggtitle('B')
#p7<-ggMarginal(p7, type = "boxplot", margins = "both", size = 5, groupColour = TRUE, groupFill = TRUE)
#p7

pdf(file = 'D:/microbiome/micom/202404/paper/figures/fig2_1.pdf',width = 15,height = 4)
grid.arrange(p5,p7,ncol=2)
dev.off()
#########
#top5_D<-rownames(abundance_D)[order(rowSums(abundance_D),decreasing = T)[1:5]]
#df_D<-rowSums(abundance_D[top5_D,])%>%as.data.frame()
#colnames(df_D)<-'Depression'
#df_D<-rbind(df_D,abundance_D[-order(rowSums(abundance_D),decreasing = T)[1:5],]%>%rowSums()%>%sum())
#rownames(df_D)[6]<-'Others'
#top5_H<-rownames(abundance_H)[order(rowSums(abundance_H),decreasing = T)[1:5]]
#df_H<-rowSums(abundance_H[top5_H,])%>%as.data.frame()
#colnames(df_H)<-'Healthy'
#df_H<-rbind(df_H,abundance_H[-order(rowSums(abundance_H),decreasing = T)[1:5],]%>%rowSums()%>%sum())
#rownames(df_H)[6]<-'Others'
#df<-merge(df_D,df_H,by = "row.names",all = TRUE)
#rownames(df)<-df$Row.names
#df$Row.names<-NULL
#df[is.na(df)]<-0
#df$genus<-rownames(df)
#temp<-colSums(df[,1:2])%>%as.numeric()
#df[,1]<-df[,1]/temp[1]
#df[,2]<-df[,2]/temp[2]
#data_long <- df %>%
#  tidyr::gather(key = "group", value = "value", -genus)
#data_alluvial <- data_long %>%
#  dplyr::mutate(flow = rep(1:length(unique(data_long$genus)),2))
#desired_order<-c('Healthy','Depression')
#colnames(data_alluvial)[1]<-'Genus'
#p8<-ggplot(data_alluvial, aes(x = group, y = value, fill = Genus, stratum = Genus, alluvium = flow)) +
#  geom_stratum(width = 0.4) + 
#  geom_flow(stat = "alluvium", lode.guidance = "forward", aes.flow = "backward") + 
#  theme_minimal() +
#  labs(x = "", y = "abundance") +
#  scale_fill_brewer(type = "qual", palette = "Set3") +
#    scale_x_discrete(limits=desired_order)+
#  theme(legend.text=element_text(face="italic"))+
#  ggtitle('D')+
#  theme(plot.title = element_text(face = 'bold'))
#p8
############Wilcoxon signed rank test
abundance_D<-abundance[,match(rownames(meta_D),colnames(abundance))]
abundance_H<-abundance[,match(rownames(meta_H),colnames(abundance))]
res<-matrix(0,nrow(abundance),2)%>%as.data.frame()
rownames(res)<-rownames(abundance)
colnames(res)<-c('p.val','fdr')
for(i in 1:nrow(abundance)){
  res[i,1]<-wilcox.test(abundance_D[i,]%>%as.numeric(),abundance_H[i,]%>%as.numeric(),paired = T)$p.val
}
res[,2]<-p.adjust(res[,1],method = 'BH')
diff_top50<-rownames(res)[order(res[,2],decreasing = F)[1:50]]
diff_top50
diff<-rownames(res)[which(res[,2]<0.05)]
#diff1<-rownames(res)[which(res[,2]>1e-5&res[,2]<0.05)]
#diff1
#exp<-abundance_clr[diff_top50,]
#label<-meta[,5,drop=F]
#ann_colors = list(disease=c(Healthy="#00AFBB",Depression="#E7B800"))
#col <- c(colorRampPalette(c("#4DBBD5FF", "white"))(abs(min(exp))*100),
#         colorRampPalette('white')(1),
#         colorRampPalette(c("white", "#F39B7FFF"))(max(exp)*100))
#p9<-pheatmap(exp,show_colnames = F,show_rownames = T,
#             annotation_col = label,
#             annotation_colors = ann_colors,
#             color = col,
#             cluster_cols = F,
#             fontfamily="serif")
#p9
#pdf(file = 'diff_genus.pdf',width = 8,height = 15)
#pheatmap(exp,show_colnames = F,show_rownames = T,
#         annotation_col = label,
#         annotation_colors = ann_colors,
#         color = col,
#         cluster_cols = F,
#         fontfamily="serif")
#dev.off()
#diff<-rownames(res)[which(res[,2]<0.005)]
#a<-c('Bifidobacterium','Blautia','Ruminococcus','Oscillibacter','Faecalibacterium','Corprococcus','Alistipes','Dialister')
#intersect(a,diff)
#diff_abundance<-abundance[diff,]
#top5_abundance<-rowSums(diff_abundance)%>%as.data.frame()
#top5_genus<-diff_abundance[rownames(top5_abundance)[order(top5_abundance[,1],decreasing = T)[1:5]],]
#top5_genus<-t(top5_genus)%>%as.data.frame()
#top5_genus$id<-rownames(top5_genus)
#top5_genus<-pivot_longer(top5_genus, 
#                  cols = -id,
#                  names_to = "Genus",
#                  values_to = "Relative Abundance" 
#)
#top5_genus$label<-meta$disease[match(top5_genus$id,rownames(meta))]
#p10<-ggplot(top5_genus, aes(x = `Relative Abundance`, y = Genus, fill = label)) +
#  geom_boxplot(outlier.shape = NA) +
#  theme_minimal() +
#  scale_fill_manual(values=c('Depression'='#E7B800', 'Healthy'='#00AFBB'))+
#  theme(axis.text.y =  element_text(face = "italic"))+
#  ggtitle('C')
#p10
#pdf(file = 'D:/microbiome/micom/202404/paper/figures/fig2.pdf',width = 12,height = 10)
#grid.arrange(arrangeGrob(p5, p7,ncol=2),p10, nrow = 2)
#dev.off()
#exp1<-abundance_clr[diff1,]
#col <- c(colorRampPalette(c("#4DBBD5FF", "white"))(abs(min(exp1))*100),
#         colorRampPalette('white')(1),
#         colorRampPalette(c("white", "#F39B7FFF"))(max(exp1)*100))
#p5<-pheatmap(exp1,show_colnames = F,show_rownames = F,
#             annotation_col = label,
#             annotation_colors = ann_colors,
#             color = col)
##########Maaslin2
#meta_D1<-meta_D
#meta_D1$disease[1:100]<-'Healthy'
#fit_data = Maaslin2(input_data     = t(abundance_D), 
#                    input_metadata = meta_D1, 
#                    min_prevalence = 0,
#                    normalization  = "NONE",
#                    max_significance = 1,
#                    output         = "D:/microbiome/micom/202404/paper/result/1-1/Maaslin2_result/Depression", 
#                    fixed_effects  = c("Host.age",'BMI'))
#
#meta_H1<-meta_H
#meta_H1$disease[1:100]<-'Depression'
#fit_data = Maaslin2(input_data     = t(abundance_H), 
#                    input_metadata = meta_H1, 
#                    min_prevalence = 0,
#                    normalization  = "NONE",
#                    max_significance = 1,
#                    output         = "D:/microbiome/micom/202404/paper/result/1-1/Maaslin2_result/Healthy", 
#                    fixed_effects  = c("Host.age",'BMI'))
############
#coef_D<-read.table('D:/microbiome/micom/202404/paper/result/1-1/Maaslin2_result/Depression/all_results.tsv',header = T)
#coef_H<-read.table('D:/microbiome/micom/202404/paper/result/1-1/Maaslin2_result/Healthy/all_results.tsv',header = T)
#coef_D<-coef_D[which(coef_D$metadata=='Host.age'),]
#coef_H<-coef_H[which(coef_H$metadata=='Host.age'),]
#coef_D<-coef_D[match(intersect(coef_D$feature,coef_H$feature),coef_D$feature),]
#coef_H<-coef_H[match(intersect(coef_D$feature,coef_H$feature),coef_H$feature),]
#coef<-cbind(coef_D$coef,coef_H$coef)%>%as.data.frame()
#colnames(coef)<-c('Depression','Healthy')
#rownames(coef)<-coef_D$feature
#coef$index<-abs(coef[,1]-coef[,2])
#genus<-c('Ruminococcus','Lactonifactor')
######
#a<-cor(t(abundance_D),meta_D$Host.age)
#b<-cor(t(abundance_H),meta_H$Host.age)
###########
Reactionabundance<-read.csv('D:/microbiome/micom/202404/paper/result/1-1/Results/ReactionAbundance.csv',row.names = 1)
#Reactionpresence<-Reactionabundance
#Reactionpresence[Reactionpresence!=0]<-1
#Subsystemabundance<-read.csv('D:/microbiome/micom/202404/paper/result/1-1/Results/SubsystemAbundance.csv',row.names = 1)
#Subsystempresence<-Subsystemabundance
#Subsystempresence[Subsystempresence!=0]<-1
#genuspresence<-abundance
#genuspresence[genuspresence!=0]<-1
#df<-rbind(genuspresence,Reactionpresence,Subsystempresence)
#df<-t(df)
#label_col<-c(rep('Genus',nrow(genuspresence)),
#         rep('Reactions',nrow(Reactionpresence)),
#         rep('Subsystem',nrow(Subsystempresence)))%>%as.data.frame()
#colnames(label_col)<-'label'
#rownames(label_col)<-colnames(df)
#label_row<-meta[,5,drop=F]
#col <-c("white", "#4DBBD5FF")
#p11<-pheatmap(df,show_colnames = F,show_rownames = F,
#             annotation_row = label_row,
#             annotation_col = label_col,
#             annotation_colors = ann_colors,
#             color = col,
#             cluster_rows = F,
#             cluster_cols = F,
#             legend = F)
#p11
########Wilcoxon signed rank test (Reaction)
res_reaction<-matrix(0,nrow(Reactionabundance),2)
colnames(res_reaction)<-c('p.val','fdr')
rownames(res_reaction)<-rownames(Reactionabundance)
Reactionabundance_D<-Reactionabundance[,match(rownames(meta_D),colnames(Reactionabundance))]
Reactionabundance_H<-Reactionabundance[,match(rownames(meta_H),colnames(Reactionabundance))]
for(i in 1:nrow(Reactionabundance)){
  print(i)
  res_reaction[i,1]<-wilcox.test(Reactionabundance_D[i,]%>%as.numeric(),Reactionabundance_H[i,]%>%as.numeric(),paired = T)$p.val
}
res_reaction[,2]<-p.adjust(res_reaction[,1],method = 'BH')
diff_reaction<-rownames(res_reaction)[which(res_reaction[,2]<0.05)]
diff_reaction
df<-t(Reactionabundance[rownames(res_reaction)[order(res_reaction[,1],decreasing = F)[1:20]],,drop=F])%>%as.data.frame()
df$label<-meta$disease[match(rownames(df),rownames(meta))]
df1<-pivot_longer(df, 
                  cols = -label,
                  names_to = "reaction",
                  values_to = "value" 
)
p12<-ggplot(df1, aes(x = value, y = reaction, fill = label)) +
      geom_boxplot(outlier.shape = NA) +
      theme_minimal() +
      scale_x_break(c(0.1,0.99), scales = 0.5,ticklabels=c(0.99,0.995,1))+
      scale_fill_manual(values=c('Depression'='#E7B800', 'Healthy'='#00AFBB'))+
      scale_x_continuous(labels = function(x) sprintf("%.2f", x))+
  ggtitle('D')+
  theme(plot.title = element_text(face = 'bold'))
p12


#pdf(file = 'diff_reactions.pdf',width = 8,height = 15)
#pheatmap(exp,show_colnames = F,show_rownames = T,
#         annotation_col = label,
#         annotation_colors = ann_colors,
#         color = col,
#         cluster_cols = F)
#dev.off()
########Wilcoxon signed rank test (Subsystem)
#res_Subsystem<-matrix(0,nrow(Subsystemabundance),2)
#colnames(res_Subsystem)<-c('p.val','fdr')
#rownames(res_Subsystem)<-rownames(Subsystemabundance)
#Subsystemabundance_D<-Subsystemabundance[,match(rownames(meta_D),colnames(Subsystemabundance))]
#Subsystemabundance_H<-Subsystemabundance[,match(rownames(meta_H),colnames(Subsystemabundance))]
#for(i in 1:nrow(Subsystemabundance)){
#  print(i)
#  res_Subsystem[i,1]<-wilcox.test(Subsystemabundance_D[i,]%>%as.numeric(),Subsystemabundance_H[i,]%>%as.numeric(),paired = T)$p.val
#}
#res_Subsystem[,2]<-p.adjust(res_Subsystem[,1],method = 'BH')
#diff_Subsystem<-rownames(res_Subsystem)[which(res_Subsystem[,1]<0.05)]
#diff_Subsystem
#exp<-Subsystemabundance[diff_Subsystem,]
#col <- c(colorRampPalette('white')(1),
#         colorRampPalette(c("white", "#F39B7FFF"))(max(exp)*100))
#p13<-pheatmap(exp,show_colnames = F,show_rownames = T,
#             annotation_col = label_row,
#             annotation_colors = ann_colors,
#             color = col,
#             cluster_cols = F)
#p13
############alpha/beta(reaction)
temp<-colSums(Reactionabundance)%>%as.numeric()
for(i in 1:ncol(Reactionabundance)){
  Reactionabundance[,i]<-Reactionabundance[,i]/temp[i]
}
OTU<-phyloseq::otu_table(Reactionabundance,taxa_are_rows = T)
#rownames(meta)<-meta$sampleid
meta1<-phyloseq::sample_data(meta)
physeq<-phyloseq(OTU,meta1)
mpse<-as.MPSE(physeq)
mpse@assays@data@listData[["Abundance"]]<-mpse@assays@data@listData[["Abundance"]]*1e20
mpse %<>% 
  mp_decostand(.abundance=Abundance)
cols<-c('Depression'='#E7B800', 'Healthy'='#00AFBB')
mpse %<>% mp_cal_alpha(.abundance = Abundance,force = T)
p14 <- mpse %>% 
  mp_plot_alpha(
    .alpha = c(Observe, Shannon,Chao1,Simpson, Pielou),
    .group = disease,
  ) +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) +
  theme(
    legend.position="none",
    strip.background = element_rect(colour=NA, fill="grey")
  )+
  ggtitle('A')+
  theme(plot.title = element_text(face = 'bold'))

mpse %>% 
  mp_anosim(.abundance=Abundance, .group=disease, action="get")

p14

dist_bray <- phyloseq::distance(physeq, method = "bray")
pcoa <- ordinate(physeq, method = "PCoA", distance = dist_bray)
pcoa_df <- data.frame(pcoa$vectors)
pcoa_df$SampleID <- rownames(pcoa_df)
sample_data<-cbind(rownames(meta),meta$disease)%>%as.data.frame()
colnames(sample_data)<-c('SampleID','disease')
rownames(sample_data)<-sample_data$SampleID
pcoa_df <- merge(pcoa_df, sample_data, by.x = "SampleID", by.y = "SampleID")
p15<-ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = disease)) +
  geom_point(alpha = 0.4) +
  labs(
    x = paste0("PCOA1 (", round(pcoa[["values"]][["Eigenvalues"]][1], 2), "%)"),
    y = paste0("PCOA2 (", round(pcoa[["values"]][["Eigenvalues"]][2], 2), "%)"),
    color = "label"
  ) +
  theme_minimal()+
  theme(legend.position = 'none')+
  stat_ellipse(aes(color = disease),type = "t", level = 0.95)+
  scale_color_manual(values=c('Depression'='#E7B800', 'Healthy'='#00AFBB'))+
  annotate("text", x = -0.1, y = -0.2, label = "ANOSIM R: -0.0003651 
           P: 0.493")+
  ggtitle('B')+
  theme(plot.title = element_text(face = 'bold'))
p15
p15<-ggMarginal(p15, type = "boxplot", margins = "both", size = 5, groupColour = TRUE, groupFill = TRUE)
p15
#abundance1<-t(Reactionabundance)%>%as.data.frame()
#env<-cbind(meta$Host.age,meta$BMI)%>%as.data.frame()
#colnames(env)<-c('Age','BMI')
#rownames(env)<-rownames(abundance1)
#cca <- cca(abundance1 ~ ., data = env)
#cca_summary<-summary(cca)
#cca_df<-cca_summary$sites%>%as.data.frame()
#cca_df$disease<-meta$disease[match(rownames(cca_df),rownames(meta))]
#scores <- as.data.frame(scores(cca)$biplot)
#p15<-ggplot(cca_df, aes(x = CCA1, y = CCA2, color = disease)) +
#  geom_point(alpha = 0.4) +
#  geom_segment(data = scores, aes(x = 0, xend = CCA1*20, y = 0, yend = CCA2*20), 
#               arrow = arrow(length = unit(0.2, "cm")), color = "red")+
#  labs(
#    x = paste0("CCA1 (", round(cca$CCA$eig[1] * 100, 2), "%)"),
#    y = paste0("CCA2 (", round(cca$CCA$eig[2] * 100, 2), "%)"),
#    color = "label"
#  ) +
#  theme_minimal()+
#  theme(legend.position = 'none')+
#  stat_ellipse(aes(color = disease),type = "t", level = 0.95)+
#  scale_color_manual(values=c('Depression'='#E7B800', 'Healthy'='#00AFBB'))+
#  annotate("text", x = -10, y = -20, label = "ANOSIM R: 0.007962 
#           P: 0.004")+
#  ggtitle('B')
#p15<-ggMarginal(p15, type = "boxplot", margins = "both", size = 5, groupColour = TRUE, groupFill = TRUE)
#p15
#########
fc_reaction<-rowMeans(Reactionabundance_D)/rowMeans(Reactionabundance_H)%>%as.data.frame()
res_reaction<-cbind(res_reaction,fc_reaction)
colnames(res_reaction)[3]<-'log2fc'
res_reaction[which(res_reaction$log2fc!=0),3]<-log2(res_reaction[which(res_reaction$log2fc!=0),3])
res_reaction$`-log10p.val`<-(-log10(res_reaction$p.val))

cut_off_log2FC<-1
cut_off_p<-0.05
res_reaction$sig<-ifelse(res_reaction$p.val < cut_off_p &
         abs(res_reaction$log2fc) >= cut_off_log2FC,  
       ifelse(res_reaction$log2fc > cut_off_log2FC ,'Up','Down'),'No significance')
p16<-ggplot(res_reaction, aes(x =log2fc, y=`-log10p.val`,colour = sig)) + 
  geom_point(alpha=0.65, size=2) +  
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) +  
  geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=4,col="black",lwd=0.8) + 
  geom_hline(yintercept = -log10(cut_off_p), lty=4,col="black",lwd=0.8) +  
  labs(x="log2FC", y="-log10(p.val)") +  
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0),
        legend.position="right", 
        legend.title = element_blank()
  )+
  ggtitle('C')+
  theme(plot.title = element_text(face = 'bold'))
p16

############alpha/beta(Subsystemabundance)
#temp<-colSums(Subsystemabundance)%>%as.numeric()
#for(i in 1:ncol(Subsystemabundance)){
#  Subsystemabundance[,i]<-Subsystemabundance[,i]/temp[i]
#}
#OTU<-phyloseq::otu_table(Subsystemabundance,taxa_are_rows = T)
#meta1<-phyloseq::sample_data(meta)
#physeq<-phyloseq(OTU,meta1)
#mpse<-as.MPSE(physeq)
#mpse@assays@data@listData[["Abundance"]]<-mpse@assays@data@listData[["Abundance"]]*1e20
#mpse %<>% 
#  mp_decostand(.abundance=Abundance)
#cols<-c('Depression'='#E7B800', 'Healthy'='#00AFBB')
#mpse %<>% mp_cal_alpha(.abundance = Abundance,force = T)
#p16 <- mpse %>% 
#  mp_plot_alpha(
#    .alpha = c(Observe, Shannon,Chao1,Simpson, Pielou),
#    .group = disease,
#  ) +
#  scale_fill_manual(values=cols) +
#  scale_color_manual(values=cols) +
#  theme(
#    legend.position="none",
#    strip.background = element_rect(colour=NA, fill="grey")
#  )
#p16
#dist_bray <- phyloseq::distance(physeq, method = "bray")
#pcoa <- ordinate(physeq, method = "PCoA", distance = dist_bray)
#pcoa_df <- data.frame(pcoa$vectors)
#pcoa_df$SampleID <- rownames(pcoa_df)
#sample_data<-cbind(rownames(meta),meta$disease)%>%as.data.frame()
#colnames(sample_data)<-c('SampleID','disease')
#rownames(sample_data)<-sample_data$SampleID
#pcoa_df <- merge(pcoa_df, sample_data, by.x = "SampleID", by.y = "SampleID")
#p17<-ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = disease)) +
#  geom_point(size = 2) +
#  labs(x = paste0("PCoA1 (", round(pcoa$values$Relative_eig[1] * 100, 2), "%)"),
#       y = paste0("PCoA2 (", round(pcoa$values$Relative_eig[2] * 100, 2), "%)")) +
#  theme_minimal() +
#  stat_ellipse(aes(color = disease),type = "t", level = 0.95) +
#  scale_color_manual(values=c('Depression'='#E7B800', 'Healthy'='#00AFBB'))
#p17<-ggMarginal(p17, type = "boxplot", margins = "both", size = 5, groupColour = TRUE, groupFill = TRUE)
#p17
###############
#abundance_H<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/healthy_genus.csv',row.names = 1)
#files<-list.files('D:/microbiome/micom/genus_models/genus_models')
#for(i in 1:length(files)){
#  files[i]<-substr(files[i],1,nchar(files[i])-4)
#}
#abundance_H<-abundance_H[intersect(rownames(abundance_H),files),]
#sampleids<-list.files('D:/microbiome/micom/202404/paper/result/1-1/FBA_result/healthy')
#for(i in 1:length(sampleids)){
#  sampleids[i]<-substr(sampleids[i],1,nchar(sampleids[i])-4)
#}
#for(sampleid in sampleids){
#  print(sampleid)
#  flux_file<-paste0('D:/microbiome/micom/202404/paper/result/1-1/FBA_result/healthy/',sampleid,'.csv')
#  flux<-read.csv(flux_file,row.names = 1)
#  flux<-flux[which(flux[,1]!=0),,drop=F]
#  temp<-abundance_H[,match(colnames(flux),colnames(abundance_H)),drop=F]
#  temp<-temp[which(temp[,1]!=0),,drop=F]
#  flux$genus<-0
#  flux$metabolite<-0
#  for(i in 1:nrow(temp)){
#    temp1<-grep(rownames(temp)[i],rownames(flux))
#    flux[temp1,2]<-rownames(temp)[i]
#  }
#  for(i in 1:nrow(flux)){
#    if(flux$genus[i]!='0'){
#      flux$metabolite[i]<-substr(rownames(flux)[i],nchar(flux$genus[i])+1,nchar(rownames(flux)[i]))
#    }else{
#      flux$metabolite[i]<-rownames(flux)[i]
#    }
#  }
#  flux$genus[which(flux$genus=='0')]<-'thylakoid'
#  metabolite<-unique(flux$metabolite)
#  genus<-unique(flux$genus)
#  temp3<-matrix(0,length(metabolite),1)%>%as.data.frame()
#  rownames(temp3)<-metabolite
#  temp_a<-abundance_H[,sampleid,drop=F]
#  temp_a<-temp_a/colSums(temp_a)
#  for(i in 1:length(genus)){
#    if(genus[i]!='thylakoid'){
#      flux[which(flux$genus==genus[i]),1]<-flux[which(flux$genus==genus[i]),1]*temp_a[genus[i],]
#    }
#  }
#  for(i in 1:length(metabolite)){
#    temp3[i,1]<-flux[which(flux$metabolite==metabolite[i]),1]%>%sum()
#  }
#  colnames(temp3)<-sampleid
#  output_file<-paste0('D:/microbiome/micom/202404/paper/result/1-1//FBA_result_p/healthy/',sampleid,'.csv')
#  write.csv(temp3,output_file)
#}

#temp4<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/FBA_result_p/healthy/',sampleids[1],'.csv'),row.names = 1)
#temp5<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/FBA_result_p/healthy/',sampleids[2],'.csv'),row.names = 1)
#merged_df <- merge(temp4, temp5, by = "row.names", all = TRUE)
#merged_df[is.na(merged_df)] <- 0
#rownames(merged_df) <- merged_df$Row.names
#merged_df$Row.names <- NULL
#for(i in 3:length(sampleids)){
#  print(i)
#  temp6<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/FBA_result_p/healthy/',sampleids[i],'.csv'),row.names = 1)
#  merged_df <- merge(merged_df, temp6, by = "row.names", all = TRUE)
#  merged_df[is.na(merged_df)] <- 0
#  rownames(merged_df) <- merged_df$Row.names
#  merged_df$Row.names <- NULL
#}
#write.csv(merged_df,'D:/microbiome/micom/202404/paper/result/1-1/metabolite_H.csv')
########
#abundance_D<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/depression_genus.csv',row.names = 1)
#files<-list.files('D:/microbiome/micom/genus_models/genus_models')
#for(i in 1:length(files)){
#  files[i]<-substr(files[i],1,nchar(files[i])-4)
#}
#abundance_D<-abundance_D[intersect(rownames(abundance_D),files),]
#sampleids<-list.files('D:/microbiome/micom/202404/paper/result/1-1/FBA_result/depression')
#for(i in 1:length(sampleids)){
#  sampleids[i]<-substr(sampleids[i],1,nchar(sampleids[i])-4)
#}
#for(sampleid in sampleids){
#  print(sampleid)
#  flux_file<-paste0('D:/microbiome/micom/202404/paper/result/1-1/FBA_result/depression/',sampleid,'.csv')
#  flux<-read.csv(flux_file,row.names = 1)
#  flux<-flux[which(flux[,1]!=0),,drop=F]
#  temp<-abundance_D[,match(colnames(flux),colnames(abundance_D)),drop=F]
#  temp<-temp[which(temp[,1]!=0),,drop=F]
#  flux$genus<-0
#  flux$metabolite<-0
#  for(i in 1:nrow(temp)){
#    temp1<-grep(rownames(temp)[i],rownames(flux))
#    flux[temp1,2]<-rownames(temp)[i]
#  }
#  for(i in 1:nrow(flux)){
#    if(flux$genus[i]!='0'){
#      flux$metabolite[i]<-substr(rownames(flux)[i],nchar(flux$genus[i])+1,nchar(rownames(flux)[i]))
#    }else{
#      flux$metabolite[i]<-rownames(flux)[i]
#    }
#  }
#  flux$genus[which(flux$genus=='0')]<-'thylakoid'
#  metabolite<-unique(flux$metabolite)
#  genus<-unique(flux$genus)
#  temp3<-matrix(0,length(metabolite),1)%>%as.data.frame()
#  rownames(temp3)<-metabolite
#  temp_a<-abundance_D[,sampleid,drop=F]
#  temp_a<-temp_a/colSums(temp_a)
#  for(i in 1:length(genus)){
#    if(genus[i]!='thylakoid'){
#      flux[which(flux$genus==genus[i]),1]<-flux[which(flux$genus==genus[i]),1]*temp_a[genus[i],]
#    }
#  }
#  for(i in 1:length(metabolite)){
#    temp3[i,1]<-flux[which(flux$metabolite==metabolite[i]),1]%>%sum()
#  }
#  colnames(temp3)<-sampleid
#  output_file<-paste0('D:/microbiome/micom/202404/paper/result/1-1//FBA_result_p/depression/',sampleid,'.csv')
#  write.csv(temp3,output_file)
#}
#
#temp4<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/FBA_result_p/depression/',sampleids[1],'.csv'),row.names = 1)
#temp5<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/FBA_result_p/depression/',sampleids[2],'.csv'),row.names = 1)
#merged_df <- merge(temp4, temp5, by = "row.names", all = TRUE)
#merged_df[is.na(merged_df)] <- 0
#rownames(merged_df) <- merged_df$Row.names
#merged_df$Row.names <- NULL
#for(i in 3:length(sampleids)){
#  print(i)
#  temp6<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/FBA_result_p/depression/',sampleids[i],'.csv'),row.names = 1)
#  merged_df <- merge(merged_df, temp6, by = "row.names", all = TRUE)
#  merged_df[is.na(merged_df)] <- 0
#  rownames(merged_df) <- merged_df$Row.names
#  merged_df$Row.names <- NULL
#}
#write.csv(merged_df,'D:/microbiome/micom/202404/paper/result/1-1/metabolite_D.csv')


########
metabolite_D<-read.csv('D:/microbiome/micom/202404/paper/result/1-1/metabolite_D.csv',row.names = 1)
metabolite_H<-read.csv('D:/microbiome/micom/202404/paper/result/1-1/metabolite_H.csv',row.names = 1)
metabolite_D<-metabolite_D[intersect(rownames(metabolite_D),rownames(metabolite_H)),]
metabolite_H<-metabolite_H[intersect(rownames(metabolite_D),rownames(metabolite_H)),]
metabolite_D<-metabolite_D[,match(rownames(meta_D),colnames(metabolite_D))]
metabolite_H<-metabolite_H[,match(rownames(meta_H),colnames(metabolite_H))]
metabolite<-cbind(metabolite_D,metabolite_H)
meta_D<-read.table('D:/microbiome/micom/202404/paper/data/meta_Depression.txt',sep = '\t',header = T)
meta_H<-read.table('D:/microbiome/micom/202404/paper/data/meta_healthy.txt',sep = '\t',header = T)
colnames(meta_H)[8]<-'disease'
colnames(meta_D)[8]<-'disease'
meta_H$disease<-'Healthy'
meta_D$disease<-'Depression'
meta<-rbind(meta_D,meta_H)

res_metabolite<-matrix(0,nrow(metabolite_D),2)%>%as.data.frame()
rownames(res_metabolite)<-rownames(metabolite_D)
colnames(res_metabolite)<-c('p.val','fdr')
for(i in 1:nrow(metabolite_D)){
  print(i)
  res_metabolite[i,1]<-wilcox.test(metabolite_D[i,]%>%as.numeric(),metabolite_H[i,]%>%as.numeric(),paired = T)$p.val
}
res_metabolite[,2]<-p.adjust(res_metabolite[,1],method = 'BH')
diff_metabolite<-rownames(res_metabolite)[which(res_metabolite[,2]<0.05)]
diff_metabolite
df<-t(metabolite[rownames(res_metabolite)[which(res_metabolite[,2]<0.05)],,drop=F])%>%as.data.frame()
df<-scale(df,scale = T,center = T)%>%as.data.frame()
metabolite_all<-fread('D:/microbiome/micom/202404/paper/data/MetaboliteDatabase.txt')%>%as.data.frame()
compartment<-c()
for(i in 1:ncol(df)){
  compartment<-c(compartment,substr(colnames(df)[i],nchar(colnames(df)[i])-2,nchar(colnames(df)[i])))
}
name<-c()
for(i in 1:ncol(df)){
  name<-c(name,substr(colnames(df)[i],1,nchar(colnames(df)[i])-3))
}
name1<-metabolite_all[match(name,metabolite_all[,1]),2]
name1<-paste0(name1,compartment)
colnames(df)<-name1
df$label<-meta$disease[match(rownames(df),meta$Run.ID)]


df1<-pivot_longer(df, 
                  cols = -label,
                  names_to = "metabolite",
                  values_to = "flux" 
)

#p<-list()
#for(i in 1:(ncol(df)-1)){
#  print(i)
#  temp<-df1[which(df1$metabolite==colnames(df)[i]),]
#  p[[i]]<-ggplot(temp, aes(x = label, y = flux, fill = label)) +
#    geom_boxplot(outlier.shape = NA) +
#    theme_minimal() +
#    scale_fill_manual(values=c('Depression'='#E7B800', 'Healthy'='#00AFBB'))+
#    coord_cartesian(ylim=c(quantile(temp$flux,0.1),quantile(temp$flux,0.9)+0.2))+
#    
#    ggtitle(colnames(df)[i])
#}
p18<-ggplot(df1, aes(x = flux, y = metabolite, fill = label)) +
  geom_boxplot(outlier.alpha=0.2) +
  theme_minimal() +
  scale_fill_manual(values=c('Depression'='#E7B800', 'Healthy'='#00AFBB'))+
  #xlim(c(-5,5))+
  labs(x='flux\n mmol/[gDW h](scaled)',y='')+
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = 'right',
        plot.title = element_text(hjust = -2.5,face = 'bold'))+
  ggtitle('A')
  
p18
###############
allmetabolite<-fread('D:/microbiome/micom/202404/paper/data/MetaboliteDatabase.txt')%>%as.data.frame()
diff_metabolite<-strsplit(diff_metabolite,'\\[')%>%unlist
diff_metabolite<-diff_metabolite[-grep('\\]',diff_metabolite)]%>%unique
id<-allmetabolite$V6[match(diff_metabolite,allmetabolite$V1)]
id<-id[id!='']
id
write.table(id,'D:/microbiome/micom/202404/paper/result/1-1/keggid_diff.txt',quote = F,sep = '\t',row.names = F,col.names = F)

stitch<-fread('D:/microbiome/micom/202404/paper/result/1-1/stitch_interactions (1).tsv')%>%as.data.frame()
protein<-c(stitch$`#node1`,
           stitch$node2)
names(protein)<-c(stitch$node1_external_id,stitch$node2_external_id)
protein<-protein[grep('ENSP',names(protein))]
protein<-unique(protein)
write.table(protein,'D:/microbiome/micom/202404/paper/result/1-1/protein.txt',quote = F,sep = '\t',row.names = F,col.names = F)

protein<-read.table('D:/microbiome/micom/202404/paper/result/1-1/protein.txt',sep='\t')
gene_list<-protein$V1
library(httr)
library(jsonlite)
library(stringr)
get_uniprot_id <- function(gene) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/search?query=gene_exact:",
                gene, "+AND+organism_id:9606+AND+reviewed:true&format=json&size=1")
  res <- GET(url)
  if (res$status_code != 200) return(NA)
  data <- fromJSON(content(res, "text", encoding = "UTF-8"))
  if (length(data$results) == 0) return(NA)
  return(data[["results"]][["primaryAccession"]])
}
get_fasta <- function(uniprot_id) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id, ".fasta")
  res <- GET(url)
  if (res$status_code != 200) return(NULL)
  return(content(res, "text", encoding = "UTF-8"))
}
fasta_all <- ""
for (gene in gene_list) {
  cat("Processing:", gene, "\n")
  uid <- get_uniprot_id(gene)
  if (is.na(uid)) {
    warning(paste("No UniProt ID found for", gene))
    next
  }
  fasta_seq <- get_fasta(uid)
  if (is.null(fasta_seq)) {
    warning(paste("Failed to retrieve FASTA for", uid))
    next
  }
  fasta_all <- paste0(fasta_all, fasta_seq, "\n")
}

# 5. ?????? FASTA ??????
writeLines(fasta_all, "output_sequences.fasta")
##############
#metabolite_used_D<-metabolite_D[diff_metabolite,]%>%t()%>%as.data.frame()
#metabolite_used_H<-metabolite_H[diff_metabolite,]%>%t()%>%as.data.frame()
#metabolite_used_D$label<-'Depression'
#metabolite_used_H$label<-'Healthy'
#set.seed(123)
#temp<-sample(1:nrow(metabolite_used_D),0.7*nrow(metabolite_used_D))
#train_D<-metabolite_used_D[temp,]
#train_H<-metabolite_used_H[temp,]
#data_train<-rbind(train_D,train_H)
#test_D<-metabolite_used_D[-temp,]
#test_H<-metabolite_used_H[-temp,]
#data_test<-rbind(test_D,test_H)
#train_index <- createFolds(data_train$label, k = 10)
#randomForestFit <- data_train %>% train(label ~ .,
#                                        method = "rf",
#                                        data = .,
#                                        tuneLength = 5,
#                                        trControl = trainControl(method = "cv", indexOut = train_index))
#randomForestFit
#pr_rf<-predict(randomForestFit,data_test,type = 'prob')
#library(pROC)
#roc_curve<-roc(data_test$label,pr_rf[,1])
#plot(roc_curve)
#sp.obj <- ci.sp(roc_curve, sensitivities=seq(0, 1, .01), boot.n=100) 
#plot(sp.obj, type="shape", col="#5bd1d7")
#text(0,0.5,paste0("AUC:",round(roc_curve[["auc"]],4)))










#pdf(file = 'diff_metabolite.pdf',width = 10,height = 5)
#ggplot(df1, aes(x = metabolite, y = value, fill = label)) +
#  geom_boxplot() +
#  theme_minimal() +
#  scale_fill_manual(values=c('Depression'='#E7B800', 'Healthy'='#00AFBB'))
#dev.off()
#######34dhphe[e]
#sampleids<-list.files('D:/microbiome/micom/202404/paper/result/1-1/FBA_result/depression')
#for(i in 1:length(sampleids)){
#  sampleids[i]<-substr(sampleids[i],1,nchar(sampleids[i])-4)
#}
#temp7<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/FBA_result/depression/',sampleids[1],'.csv'),row.names = 1)
#temp7<-temp7[grep('34dhphe\\[e\\]',rownames(temp7)),,drop=F]
#temp8<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/FBA_result/depression/',sampleids[2],'.csv'),row.names = 1)
#temp8<-temp8[grep('34dhphe\\[e\\]',rownames(temp8)),,drop=F]
#
#merged_df <- merge(temp7, temp8, by = "row.names", all = TRUE)
#merged_df[is.na(merged_df)] <- 0
#rownames(merged_df) <- merged_df$Row.names
#merged_df$Row.names <- NULL
#for(i in 3:length(sampleids)){
#  print(i)
#  temp9<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/FBA_result/depression/',sampleids[i],'.csv'),row.names = 1)
#  temp9<-temp9[grep('34dhphe\\[e\\]',rownames(temp9)),,drop=F]
#  merged_df <- merge(merged_df, temp9, by = "row.names", all = TRUE)
#  merged_df[is.na(merged_df)] <- 0
#  rownames(merged_df) <- merged_df$Row.names
#  merged_df$Row.names <- NULL
#}
#write.csv(merged_df,'D:/microbiome/micom/202404/paper/result/1-1/34dhphe[e]_D.csv')
#######
#L_Dopa_D<-read.csv('D:/microbiome/micom/202404/paper/result/1-1/34dhphe[e]_D.csv',row.names = 1)
#L_Dopa_H<-read.csv('D:/microbiome/micom/202404/paper/result/1-1/34dhphe[e]_H.csv',row.names = 1)
#for(i in 1:nrow(L_Dopa_D)){
#  rownames(L_Dopa_D)[i]<-substr(rownames(L_Dopa_D)[i],1,nchar(rownames(L_Dopa_D)[i])-10)
#}
#for(i in 1:nrow(L_Dopa_H)){
#  rownames(L_Dopa_H)[i]<-substr(rownames(L_Dopa_H)[i],1,nchar(rownames(L_Dopa_H)[i])-10)
#}
#L_Dopa_H<-abs(L_Dopa_H)
#L_Dopa_D<-abs(L_Dopa_D)
#top20_D<-rownames(L_Dopa_D)[order(rowSums(L_Dopa_D),decreasing = T)[1:20]]
#df_D<-rowSums(L_Dopa_D[top20_D,])%>%as.data.frame()
#colnames(df_D)<-'Depression'
#df_D<-rbind(df_D,L_Dopa_D[-order(rowSums(L_Dopa_D),decreasing = T)[1:20],]%>%rowSums()%>%sum())
#rownames(df_D)[21]<-'Others'
#top20_H<-rownames(L_Dopa_H)[order(rowSums(L_Dopa_H),decreasing = T)[1:20]]
#df_H<-rowSums(L_Dopa_H[top20_H,])%>%as.data.frame()
#colnames(df_H)<-'Healthy'
#df_H<-rbind(df_H,L_Dopa_H[-order(rowSums(L_Dopa_H),decreasing = T)[1:20],]%>%rowSums()%>%sum())
#rownames(df_H)[21]<-'Others'
#df<-merge(df_D,df_H,by = "row.names",all = TRUE)
#rownames(df)<-df$Row.names
#df$Row.names<-NULL
#df[is.na(df)]<-0
#df$genus<-rownames(df)
##temp<-colSums(df[,1:2])%>%as.numeric()
##df[,1]<-df[,1]/temp[1]
##df[,2]<-df[,2]/temp[2]
#data_long <- df %>%
#  tidyr::gather(key = "group", value = "value", -genus)
#data_alluvial <- data_long %>%
#  dplyr::mutate(flow = rep(1:length(unique(data_long$genus)),2))
#desired_order<-c('Healthy','Depression')
#col<-rgb(runif(21), runif(21), runif(21))
#p19<-ggplot(data_alluvial, aes(x = group, y = value, fill = genus, stratum = genus, alluvium = flow)) +
#  geom_stratum(width = 0.4) + 
#  geom_flow(stat = "alluvium", lode.guidance = "forward", aes.flow = "backward") + 
#  theme_minimal() +
#  labs(x = "", y = "L-Dopa[e]") +
#  scale_x_discrete(limits=desired_order)+
#  ggtitle('B')
#p19
#
#
#tmao_D<-read.csv('D:/microbiome/micom/202404/paper/result/1-1/tmao[e]_D.csv',row.names = 1)
#tmao_H<-read.csv('D:/microbiome/micom/202404/paper/result/1-1/tmao[e]_H.csv',row.names = 1)
#for(i in 1:nrow(tmao_D)){
#  rownames(tmao_D)[i]<-substr(rownames(tmao_D)[i],1,nchar(rownames(tmao_D)[i])-7)
#}
#for(i in 1:nrow(tmao_H)){
#  rownames(tmao_H)[i]<-substr(rownames(tmao_H)[i],1,nchar(rownames(tmao_H)[i])-7)
#}
#tmao_H<-abs(tmao_H)
#tmao_D<-abs(tmao_D)
#top20_D<-rownames(tmao_D)[order(rowSums(tmao_D),decreasing = T)[1:20]]
#df_D<-rowSums(tmao_D[top20_D,])%>%as.data.frame()
#colnames(df_D)<-'Depression'
#df_D<-rbind(df_D,tmao_D[-order(rowSums(tmao_D),decreasing = T)[1:20],]%>%rowSums()%>%sum())
#rownames(df_D)[21]<-'Others'
#top20_H<-rownames(tmao_H)[order(rowSums(tmao_H),decreasing = T)[1:20]]
#df_H<-rowSums(tmao_H[top20_H,])%>%as.data.frame()
#colnames(df_H)<-'Healthy'
#df_H<-rbind(df_H,tmao_H[-order(rowSums(tmao_H),decreasing = T)[1:20],]%>%rowSums()%>%sum())
#rownames(df_H)[21]<-'Others'
#df<-merge(df_D,df_H,by = "row.names",all = TRUE)
#rownames(df)<-df$Row.names
#df$Row.names<-NULL
#df[is.na(df)]<-0
#df$genus<-rownames(df)
##temp<-colSums(df[,1:2])%>%as.numeric()
##df[,1]<-df[,1]/temp[1]
##df[,2]<-df[,2]/temp[2]
#data_long <- df %>%
#  tidyr::gather(key = "group", value = "value", -genus)
#data_alluvial <- data_long %>%
#  dplyr::mutate(flow = rep(1:length(unique(data_long$genus)),2))
#desired_order<-c('Healthy','Depression')
#col<-rgb(runif(22), runif(22), runif(22))
#p20<-ggplot(data_alluvial, aes(x = group, y = value, fill = genus, stratum = genus, alluvium = flow)) +
#  geom_stratum(width = 0.4) + 
#  geom_flow(stat = "alluvium", lode.guidance = "forward", aes.flow = "backward") + 
#  theme_minimal() +
#  labs(x = "", y = "Trimethylamine N-oxide[e]") +
#  scale_x_discrete(limits=desired_order)+
#  ggtitle('C')
#p20
#
#pdf(file = 'D:/microbiome/micom/202404/paper/figures/fig4.pdf',width = 10,height = 7)
#grid.arrange(p18,arrangeGrob(p19, p20,ncol=2), nrow = 2)
#dev.off()
#
########
#L_Dopa_D<-read.csv('D:/microbiome/micom/202404/paper/result/1-1/34dhphe[e]_D.csv',row.names = 1)
#L_Dopa_H<-read.csv('D:/microbiome/micom/202404/paper/result/1-1/34dhphe[e]_H.csv',row.names = 1)
#L_Dopa <- merge(L_Dopa_D, L_Dopa_H, by = "row.names", all = TRUE)
#L_Dopa[is.na(L_Dopa)] <- 0
#rownames(L_Dopa) <- L_Dopa$Row.names
#L_Dopa$Row.names <- NULL
#for(i in 1:nrow(L_Dopa)){
#  rownames(L_Dopa)[i]<-substr(rownames(L_Dopa)[i],1,nchar(rownames(L_Dopa)[i])-10)
#}
#L_Dopa_D<-L_Dopa[,match(colnames(L_Dopa_D),colnames(L_Dopa))]
#L_Dopa_H<-L_Dopa[,match(colnames(L_Dopa_H),colnames(L_Dopa))]
#res_L_Dopa<-matrix(0,nrow(L_Dopa),2)%>%as.data.frame()
#rownames(res_L_Dopa)<-rownames(L_Dopa)
#colnames(res_L_Dopa)<-c('p.val','fdr')
#for(i in 1:nrow(L_Dopa_D)){
#  res_L_Dopa[i,1]<-wilcox.test(L_Dopa_D[i,]%>%as.numeric(),L_Dopa_H[i,]%>%as.numeric(),paired = T)$p.val  
#}
#res_L_Dopa[,2]<-p.adjust(res_L_Dopa[,1],method = 'BH')
#diff_L_Dopa<-rownames(res_L_Dopa)[which(res_L_Dopa[,2]<0.1)]
#setdiff(diff_L_Dopa,diff)
##########
#tmao_D<-read.csv('D:/microbiome/micom/202404/paper/result/1-1/tmao[e]_D.csv',row.names = 1)
#tmao_H<-read.csv('D:/microbiome/micom/202404/paper/result/1-1/tmao[e]_H.csv',row.names = 1)
#tmao <- merge(tmao_D, tmao_H, by = "row.names", all = TRUE)
#tmao[is.na(tmao)] <- 0
#rownames(tmao) <- tmao$Row.names
#tmao$Row.names <- NULL
#for(i in 1:nrow(tmao)){
#  rownames(tmao)[i]<-substr(rownames(tmao)[i],1,nchar(rownames(tmao)[i])-7)
#}
#tmao_D<-tmao[,match(colnames(tmao_D),colnames(tmao))]
#tmao_H<-tmao[,match(colnames(tmao_H),colnames(tmao))]
#res_tmao<-matrix(0,nrow(tmao),2)%>%as.data.frame()
#rownames(res_tmao)<-rownames(tmao)
#colnames(res_tmao)<-c('p.val','fdr')
#for(i in 1:nrow(tmao_D)){
#  res_tmao[i,1]<-wilcox.test(tmao_D[i,]%>%as.numeric(),tmao_H[i,]%>%as.numeric(),paired = T)$p.val  
#}
#res_tmao[,2]<-p.adjust(res_tmao[,1],method = 'BH')
#diff_tmao<-rownames(res_tmao)[which(res_tmao[,1]<0.05)]
#setdiff(diff_tmao,diff)



#######FEA
reaction_all<-read.csv('D:/microbiome/micom/202404/paper/data/ReactionDatabase.csv',row.names = 1)
reaction_used<-reaction_all[match(rownames(res_reaction),rownames(reaction_all)),,drop=F]
diff_reaction_used<-reaction_all[match(diff_reaction,rownames(reaction_all)),,drop=F]
reaction_used<-na.omit(reaction_used)
reaction_used$Subsystem[which(reaction_used$Subsystem=='')]<-'Other'
Subsystem_list <- split(rownames(reaction_used), reaction_used$Subsystem)
all_reactions<-rownames(reaction_used)
selected_reactions<-rownames(diff_reaction_used)
res_FEA <- data.frame(
  Subsystem = character(),
  PValue = numeric(),
  stringsAsFactors = FALSE
)
N <- length(all_reactions)
M <- length(selected_reactions)
for (Subsystem in names(Subsystem_list)) {
  K <- length(Subsystem_list[[Subsystem]])
  x <- sum(selected_reactions %in% Subsystem_list[[Subsystem]])
  p_value <- phyper(x - 1, K, N - K, M, lower.tail = FALSE)
  res_FEA <- rbind(res_FEA, data.frame(Subsystem = Subsystem, PValue = p_value))
}
res_FEA <- res_FEA %>%
  mutate(AdjustedPValue = p.adjust(PValue, method = "BH")) %>%
  arrange(AdjustedPValue)
res_FEA_sig<-res_FEA[which(res_FEA$AdjustedPValue<0.05),]
index<-match(res_FEA_sig$Subsystem,names(Subsystem_list))
ratio<-c()
for(i in 1:length(index)){
  ratio<-c(ratio,length(which(diff_reaction_used$Subsystem==names(Subsystem_list)[[index[i]]]))/length(Subsystem_list[[index[i]]]))
}
res_FEA_sig$neg_log_padj <- -log10(res_FEA_sig$AdjustedPValue)
res_FEA_sig$ratio<-ratio
p21<-ggplot(res_FEA_sig, aes(x = ratio, y = Subsystem, size = neg_log_padj, color = AdjustedPValue)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(name = "-log10(p.adjust)") +
  scale_color_gradient(low = "#4DBBD5FF", high = "#F39B7FFF", name = "fdr") +
  labs(title = "",
       x = "Reaction Ratio",
       y = "Subsystem") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(hjust = -1.05,face = 'bold'))+
  guides(color = guide_legend(ncol = 2),
         size = guide_legend(ncol = 2))+
  ggtitle('E')
  
p21
a<-ggplot()+theme_minimal()
pdf(file = 'D:/microbiome/micom/202404/paper/figures/fig3_12.pdf',width = 16,height = 14)
grid.arrange(arrangeGrob(p14, p15,nrow=2),arrangeGrob(p16,a,p21,nrow=3,heights = c(0.3,0.4,0.3)), ncol = 2)
dev.off()
pdf(file = 'D:/microbiome/micom/202404/paper/figures/p12.pdf',width = 9,height = 6)
p12
dev.off()

##########
#metabolite_D<-read.csv('D:/microbiome/micom/202404/paper/result/1-1/metabolite_D.csv',row.names = 1)
#metabolite_H<-read.csv('D:/microbiome/micom/202404/paper/result/1-1/metabolite_H.csv',row.names = 1)
#merged_df <- merge(metabolite_D, metabolite_H, by = "row.names", all = TRUE)
#merged_df[is.na(merged_df)] <- 0
#rownames(merged_df) <- merged_df$Row.names
#merged_df$Row.names <- NULL
#merged_df<-abs(merged_df)
#OTU<-phyloseq::otu_table(merged_df,taxa_are_rows = T)
#meta1<-phyloseq::sample_data(meta)
#physeq<-phyloseq(OTU,meta1)
#mpse<-as.MPSE(physeq)
#mpse@assays@data@listData[["Abundance"]]<-mpse@assays@data@listData[["Abundance"]]*1e10
#mpse %<>% 
#  mp_decostand(.abundance=Abundance)
#cols<-c('Depression'='#E7B800', 'Healthy'='#00AFBB')
#mpse %<>% mp_cal_alpha(.abundance = Abundance,force = T)
#p22 <- mpse %>% 
#  mp_plot_alpha(
#    .alpha = c(Observe, Shannon,Chao1,Simpson, Pielou),
#    .group = disease,
#  ) +
#  scale_fill_manual(values=cols) +
#  scale_color_manual(values=cols) +
#  theme(
#    legend.position="none",
#    strip.background = element_rect(colour=NA, fill="grey")
#  )+
#  ggtitle('A')
#p22
###########anosim
#
###########pcoa
#dist_bray <- phyloseq::distance(physeq, method = 'bray')
#pcoa <- ordinate(physeq, method = "PCoA", distance = dist_bray)
#pcoa_df <- data.frame(pcoa$vectors)
#pcoa_df$SampleID <- rownames(pcoa_df)
#sample_data<-cbind(rownames(meta),meta$disease)%>%as.data.frame()
#colnames(sample_data)<-c('SampleID','disease')
#rownames(sample_data)<-sample_data$SampleID
#pcoa_df <- merge(pcoa_df, sample_data, by.x = "SampleID", by.y = "SampleID")
#pcoa_df$age<-meta$Host.age[match(pcoa_df$SampleID,rownames(meta))]
#pcoa_df$BMI<-meta$BMI[match(pcoa_df$SampleID,rownames(meta))]
#pcoa_df$age[which(pcoa_df$age>=0&pcoa_df$age<=20)]<-'0-20'
#pcoa_df$age[which(pcoa_df$age>20&pcoa_df$age<=40)]<-'20-40'
#pcoa_df$age[which(pcoa_df$age>40&pcoa_df$age<=60)]<-'40-60'
#pcoa_df$age[which(pcoa_df$age>60&pcoa_df$age<=80)]<-'60-80'
#p23<-ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = disease, shape = age, size = BMI)) +
#  geom_point(alpha = 0.4) +
#  labs(
#    x = paste0("PCoA1 (", round(pcoa$values$Relative_eig[1] * 100, 2), "%)"),
#    y = paste0("PCoA2 (", round(pcoa$values$Relative_eig[2] * 100, 2), "%)"),
#    color = "label",
#    shape = "Age",
#    size = "BMI"
#  ) +
#  theme_minimal()+
#  scale_color_manual(values=c('Depression'='#E7B800', 'Healthy'='#00AFBB'))+
#  annotate("text", x = 0.4, y = -0.4, label = "ANOSIM R: 0.007962 
#           P: 0.004")+
#  ggtitle('B')
#p23<-ggMarginal(p23, type = "boxplot", margins = "both", size = 5, groupColour = TRUE, groupFill = TRUE)
#p23
########
#a<-res[rownames(res_L_Dopa),]
#a<-na.omit(a)
#b<-res_L_Dopa[match(rownames(a),rownames(res_L_Dopa)),]
#aa<-cbind(a[,1],b[1])
#colnames(aa)<-c('p.val1','p.val2')
#write.csv(aa,'D:/microbiome/micom/202404/paper/result/1-1/res_genus-ldopa.csv')
#
#a<-res[rownames(res_tmao),]
#a<-na.omit(a)
#b<-res_tmao[match(rownames(a),rownames(res_tmao)),]
#aa<-cbind(a[,1],b[1])
#colnames(aa)<-c('p.val1','p.val2')
#write.csv(aa,'D:/microbiome/micom/202404/paper/result/1-1/res_genus-tmao.csv')
########

metabolites<-c('15dap\\[e\\]','34dhphe\\[c\\]','3c3hmp\\[c\\]','5dglcn\\[c\\]','acACP\\[c\\]',
               'C02528\\[c\\]','cinnm\\[c\\]','cu2\\[e\\]',
               'gdpddman\\[c\\]','gdpfuc\\[c\\]','idon_L\\[c\\]')

#for(i in 1:length(metabolites)){
#  print(i)
#  sampleids_D<-list.files('D:/microbiome/micom/202404/paper/result/1-1/FBA_result/depression')
#  for(j in 1:length(sampleids_D)){
#    sampleids_D[j]<-substr(sampleids_D[j],1,nchar(sampleids_D[j])-4)
#  }
#  temp7<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/FBA_result/depression/',sampleids_D[1],'.csv'),row.names = 1)
#  temp7<-temp7[grep(metabolites[i],rownames(temp7)),,drop=F]
#  temp8<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/FBA_result/depression/',sampleids_D[2],'.csv'),row.names = 1)
#  temp8<-temp8[grep(metabolites[i],rownames(temp8)),,drop=F]
#  
#  merged_df <- merge(temp7, temp8, by = "row.names", all = TRUE)
#  merged_df[is.na(merged_df)] <- 0
#  rownames(merged_df) <- merged_df$Row.names
#  merged_df$Row.names <- NULL
#  for(k in 3:length(sampleids_D)){
#    temp9<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/FBA_result/depression/',sampleids_D[k],'.csv'),row.names = 1)
#    temp9<-temp9[grep(metabolites[i],rownames(temp9)),,drop=F]
#    merged_df <- merge(merged_df, temp9, by = "row.names", all = TRUE)
#    merged_df[is.na(merged_df)] <- 0
#    rownames(merged_df) <- merged_df$Row.names
#    merged_df$Row.names <- NULL
#  }
#  write.csv(merged_df,paste0('D:/microbiome/micom/202404/paper/result/1-1/',substr(metabolites[i],1,nchar(metabolites[i])-5),'_D.csv'))
#  
#  
#  sampleids_H<-list.files('D:/microbiome/micom/202404/paper/result/1-1/FBA_result/healthy')
#  for(j1 in 1:length(sampleids_H)){
#    sampleids_H[j1]<-substr(sampleids_H[j1],1,nchar(sampleids_H[j1])-4)
#  }
#  temp7<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/FBA_result/healthy/',sampleids_H[1],'.csv'),row.names = 1)
#  temp7<-temp7[grep(metabolites[i],rownames(temp7)),,drop=F]
#  temp8<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/FBA_result/healthy/',sampleids_H[2],'.csv'),row.names = 1)
#  temp8<-temp8[grep(metabolites[i],rownames(temp8)),,drop=F]
#  
#  merged_df <- merge(temp7, temp8, by = "row.names", all = TRUE)
#  merged_df[is.na(merged_df)] <- 0
#  rownames(merged_df) <- merged_df$Row.names
#  merged_df$Row.names <- NULL
#  for(k1 in 3:length(sampleids_H)){
#    temp9<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/FBA_result/healthy/',sampleids_H[k1],'.csv'),row.names = 1)
#    temp9<-temp9[grep(metabolites[i],rownames(temp9)),,drop=F]
#    merged_df <- merge(merged_df, temp9, by = "row.names", all = TRUE)
#    merged_df[is.na(merged_df)] <- 0
#    rownames(merged_df) <- merged_df$Row.names
#    merged_df$Row.names <- NULL
#  }
#  write.csv(merged_df,paste0('D:/microbiome/micom/202404/paper/result/1-1/',substr(metabolites[i],1,nchar(metabolites[i])-5),'_H.csv')) 
#} 

#########
i<-1
temp_D<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/',substr(metabolites[i],1,nchar(metabolites[i])-5),'_D.csv'),row.names = 1)
temp_H<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/',substr(metabolites[i],1,nchar(metabolites[i])-5),'_H.csv'),row.names = 1)
meta_H<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/meta_H_used.csv',row.names = 1)
meta_D<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/meta_D_used.csv',row.names = 1)
colnames(meta_H)[5]<-'disease'
colnames(meta_D)[5]<-'disease'
meta_H$disease<-'Healthy'
meta_D$disease<-'Depression'
meta<-rbind(meta_D,meta_H)
for(j in 1:nrow(temp_D)){
  rownames(temp_D)[j]<-substr(rownames(temp_D)[j],1,nchar(rownames(temp_D)[j])-(nchar(metabolites[i])-2))
}
for(j in 1:nrow(temp_H)){
  rownames(temp_H)[j]<-substr(rownames(temp_H)[j],1,nchar(rownames(temp_H)[j])-(nchar(metabolites[i])-2))
}

temp <- merge(temp_D, temp_H, by = "row.names", all = TRUE)
temp[is.na(temp)] <- 0
rownames(temp) <- temp$Row.names
temp$Row.names <- NULL
temp<-t(temp)%>%as.data.frame()
temp<-temp[,-which(colSums(temp)==0)]
temp$sampleid<-rownames(temp)

df1<-pivot_longer(temp, 
                  cols = -sampleid,
                  names_to = "genus",
                  values_to = "flux" 
)
df1$metabolite<-substr(metabolites[i],1,nchar(metabolites[i])-5)
df1$label<-meta$disease[match(df1$sampleid,rownames(meta))]
for(i in 2:length(metabolites)){
  print(i)
  temp_D<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/',substr(metabolites[i],1,nchar(metabolites[i])-5),'_D.csv'),row.names = 1)
  temp_H<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/',substr(metabolites[i],1,nchar(metabolites[i])-5),'_H.csv'),row.names = 1)
  for(j in 1:nrow(temp_D)){
    rownames(temp_D)[j]<-substr(rownames(temp_D)[j],1,nchar(rownames(temp_D)[j])-(nchar(metabolites[i])-2))
  }
  for(j in 1:nrow(temp_H)){
    rownames(temp_H)[j]<-substr(rownames(temp_H)[j],1,nchar(rownames(temp_H)[j])-(nchar(metabolites[i])-2))
  }
  
  temp <- merge(temp_D, temp_H, by = "row.names", all = TRUE)
  temp[is.na(temp)] <- 0
  rownames(temp) <- temp$Row.names
  temp$Row.names <- NULL
  temp<-t(temp)%>%as.data.frame()
  temp<-temp[,-which(colSums(temp)==0)]
  temp$sampleid<-rownames(temp)
  
  df3<-pivot_longer(temp, 
                    cols = -sampleid,
                    names_to = "genus",
                    values_to = "flux" 
  )
  df3$metabolite<-substr(metabolites[i],1,nchar(metabolites[i])-5)
  df3$label<-meta$disease[match(df3$sampleid,rownames(meta))]
  df1<-rbind(df1,df3)
}
#pdf(file = 'test.pdf',width = 10,height = 1200)
#ggplot(df1, aes(x = label, y = flux, fill = label)) +
#  geom_boxplot(outlier.shape = NA) +
#  theme(axis.text.x = element_blank(),     # ihxh=4ff,f g->
#        axis.ticks.x = element_blank())+
#  stat_compare_means(comparisons = list(c("Depression", "Healthy")),label.y = 1.5)+
#  scale_fill_manual(values=c('Depression'='#E7B800', 'Healthy'='#00AFBB'))+
#  coord_cartesian(ylim = c(-1, 1))+
#  facet_grid(genus ~ metabolite)
#dev.off()


#write.csv(df1,'data.csv')
#########
#diff<-read.csv('D:/microbiome/micom/202404/paper/result/1-1/diff_taxon.csv',header = F)%>%as.character()
library(tidyverse)
library(RColorBrewer)
library(effsize)
#df1$flux<-scale(df1$flux,scale = T,center = T)
metabolite_used<-unique(df1$metabolite)
metabolite_all<-fread('D:/microbiome/micom/202404/paper/data/MetaboliteDatabase.txt')%>%as.data.frame()
metabolite_used1<-metabolite_all[match(metabolite_used,metabolite_all[,1]),2]
df1$label<-meta$disease[match(df1$sampleid,rownames(meta))]
meta<-df1[!duplicated(df1$sampleid),c(1,5)]
res_all<-data.frame()

for(i in 1:length(metabolite_used)){
  print(i)
  temp_m<-df1[df1$metabolite==metabolite_used[i],]
  con<-meta$sampleid[meta$label=="Healthy"]
  ds<-meta$sampleid[meta$label=="Depression"]
  m1<-df1[df1$metabolite==metabolite_used[i],]
  mat1<-pivot_wider(m1[,1:3],names_from = sampleid,values_from = flux)
  mat1<-column_to_rownames(mat1,var = "genus")
  mat1<-mat1[,meta$sampleid]
  fc1<-rowMeans(mat1[,ds])/rowMeans(mat1[,con])
  p1<-c()
  cliff_delta<-c()
  a<-c()
  for (j in 1:nrow(mat1)) {
    ptmp<-wilcox.test(as.numeric(mat1[j,ds]),as.numeric(mat1[j,con]),paired = T)$p.val
    cliff_delta_temp<-cliff.delta(as.numeric(mat1[j,ds]),as.numeric(mat1[j,con]))$estimate
    p1<-c(p1,ptmp)
    cliff_delta<-c(cliff_delta,cliff_delta_temp)
  }
  res1<-data.frame(genus=rownames(mat1),fc=fc1,p=p1,metabolite=metabolite_used1[i],cliff_delta=cliff_delta,a=rowMeans(mat1[,ds]))
  res_all<-rbind(res_all,res1)
}
res_all$sig<-'no'
res_all$sig[which(res_all$p<0.05)]<-'yes'
res_all1<-res_all[res_all$sig=='yes',]
#setdiff(res_all1$genus%>%unique(),diff_top50)
cu_genus<-res_all1$genus[which(res_all1$metabolite=='Cu2+')]

res_all2<-res_all1[which(res_all1$metabolite!='Cu2+'),]
res_cu<-res_all1[which(res_all1$metabolite=='Cu2+'),]
res_all2$a[res_all2$a>0]<-'export'
res_all2$a[res_all2$a<0]<-'import'
colnames(res_all2)[6]<-'direction'
res_all2 <- res_all2 %>%
  mutate(significance = case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  ))

p22<-ggplot2::ggplot()+
  geom_point(data = res_all2,
             aes(x = genus, y = metabolite,size = fc,color = significance,shape=direction))+
  theme_minimal()+
  #scale_fill_manual(c('*'='red','**'='green'))+
  theme(axis.text.x = element_text(angle = 90,size = 12,face = "italic"),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = 'left',
        plot.title = element_text(hjust = -2,face = 'bold'))+
  scale_y_discrete(position = 'right')+
  labs(x='',y='')+
  guides(color = guide_legend(ncol = 1),
         size = guide_legend(ncol = 2)) +
  scale_x_discrete(labels = function(y) {
    ifelse(y %in% setdiff(res_all2$genus%>%unique(),diff), 
           paste0("<span style='color:red;'>", y, "</span>"), 
           y)
  })+
  scale_shape_manual(values = c("export" = 24, "import" = 25))+
  theme(axis.text.x=element_markdown())+
  ggtitle('B')
p22
pdf(file = 'D:/microbiome/micom/202404/paper/figures/fig4.pdf',width = 15,height = 12)
grid.arrange(p18,p22, ncol = 2)
dev.off()
#ggsave("dotp.pdf",height = 15,width = 6,limitsize = FALSE)
######

i<-8
cu_D<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/',substr(metabolites[i],1,nchar(metabolites[i])-5),'_D.csv'),row.names = 1)
cu_H<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/',substr(metabolites[i],1,nchar(metabolites[i])-5),'_H.csv'),row.names = 1)
taxon<-read.csv('D:/microbiome/micom/202404/paper/data/taxa.csv')
for(j in 1:nrow(cu_D)){
  rownames(cu_D)[j]<-substr(rownames(cu_D)[j],1,nchar(rownames(cu_D)[j])-(nchar(metabolites[i])-2))
}
for(j in 1:nrow(cu_H)){
  rownames(cu_H)[j]<-substr(rownames(cu_H)[j],1,nchar(rownames(cu_H)[j])-(nchar(metabolites[i])-2))
}

cu <- merge(cu_D, cu_H, by = "row.names", all = TRUE)
cu[is.na(cu)] <- 0
rownames(cu) <- cu$Row.names
cu$Row.names <- NULL
cu_data<-cu
cu<-t(cu)%>%as.data.frame()
cu<-cu[,-which(colSums(cu)==0)]
cu<-abs(cu)
for(i in 1:nrow(cu)){
  cu[i,]<-cu[i,]/rowSums(cu)[i]
}
cu$sampleid<-rownames(cu)

df<-pivot_longer(cu, 
                  cols = -sampleid,
                  names_to = "genus",
                  values_to = "flux" 
)
df$phyla<-taxon$p__[match(df$genus,taxon$g__)]
df$flux<-abs(df$flux)
df$flux<-log2(df$flux+1)
df$metabolic<-'Cu2+'
df$group<-meta$disease[match(df$sampleid,rownames(meta))]
set.seed(1)
col<-rgb(runif(length(unique(df$phyla))),runif(length(unique(df$phyla))),runif(length(unique(df$phyla))))
names(col)<-unique(df$phyla)
ggplot(df, aes(x = sampleid, y = flux, fill = phyla)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12))+
  facet_grid(group~.,scales = "free_x",as.table = F)




#######
p_phyla<-list()
p_genus<-list()
meta_H<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/meta_H_used.csv',row.names = 1)
meta_D<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/meta_D_used.csv',row.names = 1)
colnames(meta_H)[5]<-'disease'
colnames(meta_D)[5]<-'disease'
meta_H$disease<-'Healthy'
meta_D$disease<-'Depression'
meta<-rbind(meta_D,meta_H)
i=8
for(i in 1:length(metabolites)){
  cu_D<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/',substr(metabolites[i],1,nchar(metabolites[i])-5),'_D.csv'),row.names = 1)
  cu_H<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/',substr(metabolites[i],1,nchar(metabolites[i])-5),'_H.csv'),row.names = 1)
  cu_D<-cu_D[,-which(colnames(cu_H)=='ERR1090035')]
  cu_H<-cu_H[,-which(colnames(cu_H)=='ERR1090035')]
  taxon<-read.csv('D:/microbiome/micom/202404/paper/data/taxa.csv')
  for(j in 1:nrow(cu_D)){
    rownames(cu_D)[j]<-substr(rownames(cu_D)[j],1,nchar(rownames(cu_D)[j])-(nchar(metabolites[i])-2))
  }
  for(j in 1:nrow(cu_H)){
    rownames(cu_H)[j]<-substr(rownames(cu_H)[j],1,nchar(rownames(cu_H)[j])-(nchar(metabolites[i])-2))
  }
  
  cu <- merge(cu_D, cu_H, by = "row.names", all = TRUE)
  cu[is.na(cu)] <- 0
  rownames(cu) <- cu$Row.names
  cu$Row.names <- NULL
  cu<-t(cu)%>%as.data.frame()
  cu<-cu[,-which(colSums(cu)==0)]
  cu$sampleid<-rownames(cu)
  df<-pivot_longer(cu, 
                   cols = -sampleid,
                   names_to = "genus",
                   values_to = "flux" 
  )
  df$phyla<-taxon$p__[match(df$genus,taxon$g__)]
  df$flux<-abs(df$flux)
  df$flux<-log2(df$flux+1)
  df$metabolic<-metabolites[i]
  set.seed(1)
  col<-rgb(runif(length(unique(df$phyla))),runif(length(unique(df$phyla))),runif(length(unique(df$phyla))))
  names(col)<-unique(df$phyla)
  col_genus<-rgb(runif(length(unique(df$genus))),runif(length(unique(df$genus))),runif(length(unique(df$genus))))
  names(col_genus)<-unique(df$genus)
  p_genus[[i]]<-ggplot(df, aes(x = sampleid, y = flux, fill = genus)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = col_genus)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          strip.text = element_text(size = 12),
          legend.position = 'none')+
    facet_wrap(~ metabolic)
}
for(i in 1:length(metabolites)){
  cu_D<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/',substr(metabolites[i],1,nchar(metabolites[i])-5),'_D.csv'),row.names = 1)
  cu_H<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/',substr(metabolites[i],1,nchar(metabolites[i])-5),'_H.csv'),row.names = 1)
  cu_D<-cu_D[,-which(colnames(cu_H)=='ERR1090035')]
  cu_H<-cu_H[,-which(colnames(cu_H)=='ERR1090035')]
  taxon<-read.csv('D:/microbiome/micom/202404/paper/data/taxa.csv')
  for(j in 1:nrow(cu_D)){
    rownames(cu_D)[j]<-substr(rownames(cu_D)[j],1,nchar(rownames(cu_D)[j])-(nchar(metabolites[i])-2))
  }
  for(j in 1:nrow(cu_H)){
    rownames(cu_H)[j]<-substr(rownames(cu_H)[j],1,nchar(rownames(cu_H)[j])-(nchar(metabolites[i])-2))
  }
  
  cu <- merge(cu_D, cu_H, by = "row.names", all = TRUE)
  cu[is.na(cu)] <- 0
  rownames(cu) <- cu$Row.names
  cu$Row.names <- NULL
  cu<-t(cu)%>%as.data.frame()
  cu<-cu[,-which(colSums(cu)==0)]
  cu$sampleid<-rownames(cu)
  df<-pivot_longer(cu, 
                   cols = -sampleid,
                   names_to = "genus",
                   values_to = "flux" 
  )
  df$phyla<-taxon$p__[match(df$genus,taxon$g__)]
  df$flux<-abs(df$flux)
  df$flux<-log2(df$flux+1)
  df$metabolic<-metabolites[i]
  Firm<-df[which(df$phyla=='Proteobacteria'),]
  Firm_summarized <- Firm %>%
    group_by(sampleid) %>%
    summarize(
      flux = sum(flux),
    )
  set.seed(1)
  col_phyla<-rgb(runif(length(unique(df$phyla))),runif(length(unique(df$phyla))),runif(length(unique(df$phyla))))
  names(col_phyla)<-unique(df$phyla)
  col_phyla<-c(col_phyla,'#E7B800','#00AFBB')
  names(col_phyla)[19:20]<-c('Depression','Healthy')
  col_genus<-rgb(runif(length(unique(df$genus))),runif(length(unique(df$genus))),runif(length(unique(df$genus))))
  names(col_genus)<-unique(df$genus)
  Firm_summarized_D<-Firm_summarized[match(rownames(meta)[which(meta$disease=='Depression')],Firm_summarized$sampleid),]
  Firm_summarized_H<-Firm_summarized[match(rownames(meta)[which(meta$disease=='Healthy')],Firm_summarized$sampleid),]
  levels<-c(Firm_summarized_D$sampleid[order(Firm_summarized_D$flux)],rownames(meta_H)[match(Firm_summarized_D$sampleid[order(Firm_summarized_D$flux)],rownames(meta_D))])
  ##df$sampleid<-factor(df$sampleid,levels = Firm_summarized$sampleid[order(Firm_summarized$flux)])
  df$sampleid<-factor(df$sampleid,levels = levels)
  df$label<-meta$disease[match(df$sampleid,rownames(meta))]
  df<-na.omit(df)
  #tile_data <- df %>%
  #  mutate(xmin = as.numeric(factor(sampleid)) - 0.5,
  #         xmax = as.numeric(factor(sampleid)) + 0.5,
  #         ymin = -1,
  #         ymax = -0.5)
  p_phyla[[i]]<-ggplot(df, aes(x = sampleid, y = flux, fill = phyla)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          strip.text = element_text(size = 12)) +
    scale_fill_manual(values = col_phyla)+
    facet_wrap(~ label, scales = "free_x", ncol = 1)
  
}
p25<-ggplot(df, aes(x = sampleid, y = flux, fill = phyla)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12),
        plot.title = element_text(face = 'bold')) +
  scale_fill_manual(values = col_phyla)+
  labs(y='flux\n mmol/[gDW h](scaled)',x='Sample id')+
  facet_wrap(~ label, scales = "free_x", ncol = 1)+
  ggtitle('D')
p25
########

cu_data$genus<-rownames(cu_data)
cu_data<-pivot_longer(cu_data, 
                      cols = -genus,
                      names_to = "sampleid",
                      values_to = "flux" 
)
cu_data$label<-meta$label[match(cu_data$sampleid,meta$sampleid)]
p<-list()
for(i in 1:length(cu_genus)){
  temp<-cu_data[which(cu_data$genus==cu_genus[i]),]
  p[[i]]<-ggplot(temp, aes(x = label, y = flux, fill = label)) +
    geom_violin() +
    labs( x = "", y = '') +
    facet_wrap(~ genus)+
    stat_compare_means(aes(group = label),label = "p.signif",label.x = 1.5) +
    theme(axis.text.x = element_blank(),
          strip.text = element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),#top, right, bottom, and left
          legend.position = 'none')+
    scale_fill_manual(values=c('Depression'='#E7B800', 'Healthy'='#00AFBB'))
}
pdf(file = 'D:/microbiome/micom/202404/paper/figures/fig6.pdf',width = 36,height = 20)
grid.arrange(arrangeGrob(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],ncol=9),
             arrangeGrob(p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],ncol=9),
             arrangeGrob(p[[19]],p[[20]],p[[21]],p[[22]],p[[23]],p[[24]],p[[25]],p[[26]],p[[27]],ncol=9),
             arrangeGrob(p[[28]],p[[29]],p[[30]],p[[31]],p[[32]],p[[33]],p[[34]],p[[35]],p[[36]],ncol=9),
             arrangeGrob(p[[37]],p[[38]],p[[39]],p[[40]],p[[41]],p[[42]],p[[43]],p[[44]],p[[45]],ncol=9),
             nrow = 5)
dev.off()
#########
Firm<-df[which(df$phyla=='Firmicutes'),]
Firm_summarized <- Firm %>%
  group_by(sampleid) %>%
  summarize(
    flux = sum(flux),
  )
########
i=8
cu_D<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/',substr(metabolites[i],1,nchar(metabolites[i])-5),'_D.csv'),row.names = 1)
cu_H<-read.csv(paste0('D:/microbiome/micom/202404/paper/result/1-1/',substr(metabolites[i],1,nchar(metabolites[i])-5),'_H.csv'),row.names = 1)
cu <- merge(cu_D, cu_H, by = "row.names", all = TRUE)
cu[is.na(cu)] <- 0
rownames(cu) <- cu$Row.names
cu$Row.names <- NULL
cu<-t(cu)%>%as.data.frame()
cu<-cu[,-which(colSums(cu)==0)]
cu<-abs(cu)
for(i in 1:ncol(cu)){
  colnames(cu)[i]<-substr(colnames(cu)[i],1,nchar(colnames(cu)[i])-6)
}
cu<-cu[,match(cu_genus,colnames(cu))]
cu$group<-meta$disease[match(rownames(cu),rownames(meta))]
cu<-cbind(cu[,ncol(cu)],cu[,-ncol(cu)])
colnames(cu)[1]<-'label'
cu_D<-cu[which(cu$label=='Depression'),]
cu_H<-cu[which(cu$label=='Healthy'),]
temp1<-colMeans(cu_D[,-1])
temp2<-colMeans(cu_H[,-1])
temp<-rbind(temp1,temp2)%>%as.data.frame()
temp$label<-c('Depression','Healthy')
temp<-cbind(temp[,ncol(temp)],temp[,-ncol(temp)])

p23<-ggradar(temp,
        grid.line.width = 0,
        axis.label.size= 0,
        group.line.width = 0.2,
        group.point.size = 0.2,
        plot.extent.x.sf = 1.6,
        background.circle.colour = 'white',
        grid.label.size = 0,
        grid.min = 0,
        grid.mid = 0.2,
        grid.max = 0.3,
        group.colours = c('#E7B800','#00AFBB'),
        background.circle.transparency = 0,
        legend.position = 'none')
p23
p24<-ggplot()+theme_minimal()
pdf(file = 'D:/microbiome/micom/202404/paper/figures/fig4_12.pdf',width = 15,height = 18)
grid.arrange(arrangeGrob(p18,arrangeGrob(p22,p24,nrow = 2),ncol=2),p25,nrow=2)
dev.off()

pdf(file = 'D:/microbiome/micom/202404/paper/figures/cu_radar.pdf',width = 7,height = 4)
p23
dev.off()
###########
library(mediation)
data<-rbind(abundance,metabolite)
data<-t(data)%>%as.data.frame()
for(i in 1:ncol(data)){
  colnames(data)[i]<-paste0('f_',colnames(data)[i])
}
data$group<-meta$disease[match(rownames(data),rownames(meta))]
data$gender<-meta$Sex[match(rownames(data),rownames(meta))]
data$age<-meta$Host.age[match(rownames(data),rownames(meta))]
data$BMI<-meta$BMI[match(rownames(data),rownames(meta))]
genusnames<-diff
metabolitenames<-diff_metabolite
data$group<-(as.factor(data$group)%>%as.numeric()-1)
colnames(data)<-gsub("\\[|\\]","_",colnames(data))
f1<-paste0('group','~','f_',genusnames[i])
f1<-as.formula(f1)
Y1<-glm(f1,data=data,family = binomial("probit"))
f2<-paste0('group','~','f_',metabolitenames[j],'+gender+BMI+age')
f2<-as.formula(f2)
Y2<-lm(f2,data=data)
############
library(mediation)
data_a<-abundance[diff,]
data_m<-metabolite[diff_metabolite,]
data<-rbind(data_a,data_m)
data<-t(data)%>%as.data.frame()
allfeatures<-colnames(data)
for(i in 1:ncol(data)){
  colnames(data)[i]<-paste0('features',i)
}
data$group<-meta$disease[match(rownames(data),rownames(meta))]
data$gender<-meta$Sex[match(rownames(data),rownames(meta))]
data$age<-meta$Host.age[match(rownames(data),rownames(meta))]
data$BMI<-meta$BMI[match(rownames(data),rownames(meta))]
data$group<-(as.factor(data$group)%>%as.numeric()-1)
result_mediate<-c()
for(i in 1:nrow(data_a)){
  print(i)
  for(j in (nrow(data_a)+1):(ncol(data)-4)){
    print(j)
    f1<-paste0('features',j,'~','features',i,'+gender+BMI+age')
    f1<-as.formula(f1)
    Y1<-lm(f1,data=data)
    f2<-paste0('group','~','features',i,'+','features',j,'+gender+BMI+age')
    f2<-as.formula(f2)
    Y2<-glm(f2,data=data,family = binomial("probit"))
    r_causal2 <- mediate(Y1, Y2, treat=paste0('features',i), mediator=paste0('features',j), covariates=c('gender','BMI','age'), boot=TRUE, sims=1000)
    res_temp<-c(paste0('features',i),paste0('features',j),
                r_causal2$d1.p,r_causal2$d0.p,
                r_causal2$d1,r_causal2$d0)
    result_mediate <- append(result_mediate, list(res_temp))
  }
}
features1<-c()
features2<-c()
res_m<-c()
for(i in 1:length(result_mediate)){
  if(result_mediate[[i]][3]%>%as.numeric()<0.001&
     result_mediate[[i]][4]%>%as.numeric()<0.001){
    features1<-append(features1,allfeatures[match(result_mediate[[i]][1],colnames(data))])
    features2<-append(features2,allfeatures[match(result_mediate[[i]][2],colnames(data))])
    temp<-c(allfeatures[match(result_mediate[[i]][1],colnames(data))],
            allfeatures[match(result_mediate[[i]][2],colnames(data))],
            result_mediate[[i]][3]%>%as.numeric())%>%t()%>%as.data.frame()
    res_m<-rbind(res_m,temp)
  }
}

features1<-c()
features2<-c()

for(i in 1:length(result_mediate)){
  if(result_mediate[[i]][3]%>%as.numeric()<0.001&
     result_mediate[[i]][4]%>%as.numeric()<0.001){
    features1<-append(features1,result_mediate[[i]][1])
    features2<-append(features2,result_mediate[[i]][2])
  }
}

colnames(res_m)<-c('genus','metabolites','p.val')
res_m[,3]<-as.numeric(res_m[,3])
res_m1<-res_m
#########
res_m1<-read.csv('D:/microbiome/micom/202404/paper/result/1-1/res_m2.csv',row.names = 1)
features1<-res_m1$genus
features2<-res_m1$metabolites
features1<-match(features1,allfeatures)
features1<-paste0('features',features1)
features2<-match(features2,allfeatures)
features2<-paste0('features',features2)

result_mediate1<-c()

for(i in 1:nrow(res_m1)){
  print(i)
  f0<-paste0('group','~',features1[i],'+gender+BMI+age')#d8-d;ei~h*ei+f77fe g4 
  f0<-as.formula(f0)
  Y0<-glm(f0,data=data,family = binomial("probit"))
  f1<-paste0(features2[i],'~',features1[i],'+gender+BMI+age')#d8-d;ei~h*ei+f77fe g4 
  f1<-as.formula(f1)
  Y1<-lm(f1,data=data)
  print(summary(Y1))
  f2<-paste0('group','~',features1[i],'+',features2[i],'+gender+BMI+age')#e ei~h*ei+d8-d;ei+f77fe g4 
  f2<-as.formula(f2)
  Y2<-glm(f2,data=data,family = binomial("probit"))
  r_causal2 <- mediate(Y1, Y2, treat=features1[i], mediator=features2[i], covariates=c('gender','BMI','age'), boot=TRUE, sims=1000)
  model_summary<-summary(r_causal2)
  res_temp<-c(features1[i],
              features2[i],
              r_causal2$d.avg,
              r_causal2$d.avg.p,
              r_causal2$z0.p,
              r_causal2$z1.p,
              r_causal2$z.avg,
              r_causal2$z.avg.p,
              Y1$coefficients[2]%>%as.numeric(),
              Y2$coefficients[2]%>%as.numeric(),
              Y2$coefficients[3]%>%as.numeric(),
              r_causal2$n.avg,
              Y0$coefficients[2]
              )
  result_mediate1 <- append(result_mediate1, list(res_temp))
}
features1<-c()
features2<-c()
d.avg<-c()
z.avg<-c()
coefficientsy1<-c()
coefficientsy2<-c()
coefficientsy3<-c()
d.avg.p<-c()
z.avg.p<-c()
n.avg<-c()

res_m2<-c()
for(i in 1:length(result_mediate1)){
  if(result_mediate1[[i]][8]%>%as.numeric()<0.01&result_mediate1[[i]][4]%>%as.numeric()<0.01){
    print(i)
    features1<-result_mediate1[[i]][1]
    features2<-result_mediate1[[i]][2]
    d.avg<-result_mediate1[[i]][3]
    z.avg<-result_mediate1[[i]][7]
    d.avg.p<-result_mediate1[[i]][4]
    z.avg.p<-result_mediate1[[i]][8]
    coefficientsy1<-result_mediate1[[i]][9]
    coefficientsy2<-result_mediate1[[i]][10]
    coefficientsy3<-result_mediate1[[i]][13]
    n.avg<-result_mediate1[[i]][12]
    temp<-c(allfeatures[match(features1,colnames(data))],
            allfeatures[match(features2,colnames(data))],
            coefficientsy1,
            coefficientsy2,
            coefficientsy3,
            d.avg.p,
            z.avg.p,
            d.avg,
            z.avg,
            n.avg)%>%t()%>%as.data.frame()
    colnames(temp)<-c('genus','metabolites','coefficients1','coefficients2','coefficients3','d.avg.p','z.avg.p','d.avg','z.avg','n.avg')
    res_m2<-rbind(res_m2,temp)
  }
}
res_m2<-res_m2[res_m2$n.avg%>%as.numeric()>0.1,]
result_mediate2<-c()
for(i in 1:nrow(res_m1)){
  print(i)
  f1<-paste0('group','~',features1[i],'+gender+BMI+age')
  f1<-as.formula(f1)
  Y1<-glm(f1,data=data,family = binomial("probit"))
  f2<-paste0(features2[i],'~',features1[i],'+','group','+gender+BMI+age')
  f2<-as.formula(f2)
  Y2<-lm(f2,data=data)
  r_causal2 <- mediate(Y1, Y2, treat=features1[i], mediator=features2[i], covariates=c('gender','BMI','age'), boot=TRUE, sims=1000)
  res_temp<-c(features1[i],
              features2[i],
              r_causal2$d.avg,
              r_causal2$d.avg.p,
              r_causal2$z0.p,
              r_causal2$z1.p,
              r_causal2$z.avg,
              r_causal2$z.avg.p,
              #r_causal2$n.avg.p
              Y1$coefficients[2]%>%as.numeric(),
              Y2$coefficients[3]%>%as.numeric())
  result_mediate2 <- append(result_mediate2, list(res_temp))
}

res_m3<-c()
for(i in 1:length(result_mediate2)){
    print(i)
    features1<-result_mediate2[[i]][1]
    features2<-result_mediate2[[i]][2]
    #d.avg<-result_mediate2[[i]][3]
    #z.avg<-result_mediate2[[i]][7]
    d.avg.p<-result_mediate2[[i]][4]
    #z.avg.p<-result_mediate2[[i]][8]
    #coefficientsy1<-result_mediate2[[i]][9]
    #coefficientsy2<-result_mediate2[[i]][10]
    #coefficientsy3<-result_mediate2[[i]][13]
    #n.avg<-result_mediate2[[i]][12]
    temp<-c(allfeatures[match(features1,colnames(data))],
            allfeatures[match(features2,colnames(data))],
            #coefficientsy1,
            #coefficientsy2,
            #coefficientsy3,
            d.avg.p
            #z.avg.p,
            #d.avg,
            #z.avg,
            #n.avg)
            )%>%t()%>%as.data.frame()
    colnames(temp)<-c('genus','metabolites','d.avg.p')
    res_m3<-rbind(res_m3,temp)
}





write.csv(res_m2,'D:/microbiome/micom/202404/paper/result/1-1/res_m2.csv')

#########
library(igraph)
library(ggraph)
metabolite<-unique(res_m1$metabolites)
nodetype<-data.frame(nodes=c(unique(res_m1$genus),
                             unique(res_m1$metabolites),
                             "Depression"),
                     type=c(rep(c("Independent","Mediator","Dependent"),
                                c(length(unique(res_m1$genus)),length(unique(res_m1$metabolites)),1))))
res_m1<-rbind(res_m1,data.frame(genus=metabolite,metabolites="Depression",p.val=0))

graph<-graph_from_data_frame(res_m1,vertices = nodetype)
pos<-nodetype
pos$x<-c(seq(0.1,0.1*length(unique(res_m1$genus)),0.1),seq(0.2,0.2*length(unique(res_m1$metabolites)),0.2),1.2)
pos$y<-c(rep(0.8,23),rep(0.5,11),0.2)

p26<-ggraph(graph,
          layout = "manual",
          x = pos$x,
          y = pos$y) +
  geom_edge_link0(width=1,arrow = arrow(length = unit(3, 'mm'), type = "closed"),alpha=0.5) +
  geom_node_point(aes(fill = type), shape = 21, size = 6,alpha=0.5) +
  #scale_color_manual(values = c("Independent"="orangered4","Mediator"="orchid4","Dependent"="red")) +
  scale_fill_manual(values = c("Independent"="orangered4","Mediator"="orchid4","Dependent"="red")) +
  #scale_edge_colour_gradientn(
  #  colors = c("aquamarine3","grey90"),
  #  #limits = c(0, 0.28),
  #  space = "Lab",
  #  na.value = "grey50")+
  guides(fill = guide_legend(ncol = 1))+
  theme_graph(base_family = "sans")+#theme(legend.box = 'horizontal',
  #      legend.box.just = 'top')+
  theme(plot.title = element_text(hjust = -0.065,vjust = 0.25),title = element_text(size = 12))+
  geom_node_text(aes(label = pos$nodes), size = 5, vjust = ifelse(1:35 <= 23, -0.5, 1),
                 angle = ifelse(1:35 <= 23, 45, 0),hjust = ifelse(1:35 <= 23, 0, 0.5))+
  expand_limits(y=c(0.2,1.2),x=c(0,2.6))
#p26<-p26+ggtitle('E')
p26
pdf(file = 'D:/microbiome/micom/202404/paper/figures/fig4_9.pdf',width = 15,height = 25)
grid.arrange(arrangeGrob(p18,arrangeGrob(p22,p24,nrow = 2),ncol=2),p25,nrow=2)
dev.off()
################
meta_H<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/meta_H_used.csv',row.names = 1)
meta_D<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/meta_D_used.csv',row.names = 1)
colnames(meta_H)[5]<-'disease'
colnames(meta_D)[5]<-'disease'
meta_H$disease<-'Healthy'
meta_D$disease<-'Depression'
meta<-rbind(meta_D,meta_H)
abundance_D<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/depression_genus.csv',row.names = 1)
abundance_D<-abundance_D[match(rownames(meta_D),colnames(abundance_D))]
temp<-abundance_D
temp[temp!=0]<-1
abundance_D<-abundance_D[which(rowSums(temp)>ncol(abundance_D)*0.1)%>%as.numeric(),]
abundance_H<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/healthy_genus.csv',row.names = 1)
abundance_H<-abundance_H[match(rownames(meta_H),colnames(abundance_H))]
temp<-abundance_H
temp[temp!=0]<-1
#abundance_H<-abundance_H[which(rowSums(temp)>ncol(abundance_H)*0.1)%>%as.numeric(),]
abundance<-merge(abundance_D,abundance_H,by = "row.names",all = TRUE)
abundance[is.na(abundance)]<-0
rownames(abundance)<-abundance$Row.names
abundance$Row.names<-NULL
abundance_D<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/depression_genus.csv',row.names = 1)
abundance_D<-abundance_D[match(rownames(meta_D),colnames(abundance_D))]
temp<-abundance_D
temp[temp!=0]<-1
abundance_D<-abundance_D[which(rowSums(temp)>ncol(abundance_D)*0.1)%>%as.numeric(),]
abundance_H<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/healthy_genus.csv',row.names = 1)
abundance_H<-abundance_H[match(rownames(meta_H),colnames(abundance_H))]
temp<-abundance_H
temp[temp!=0]<-1
#abundance_H<-abundance_H[which(rowSums(temp)>ncol(abundance_H)*0.1)%>%as.numeric(),]
abundance<-merge(abundance_D,abundance_H,by = "row.names",all = TRUE)
abundance[is.na(abundance)]<-0
rownames(abundance)<-abundance$Row.names
abundance$Row.names<-NULL
files <- list.files('D:/microbiome/micom/202404/paper/result/1-1/FBA_result/depression/', full.names = TRUE)
df_list <- vector("list", length(files))
for(i in seq_along(files)){
  print(i)
  sampleid <- tools::file_path_sans_ext(basename(files[i]))
  temp_m <- read.csv(files[i], row.names = 1)
  colnames(temp_m) <- 'flux'
  temp_m$sampleid <- sampleid
  temp_m$genus <- NA
  temp_m$metabolite <- NA
  temp_a <- abundance[which(abundance[, sampleid, drop = FALSE] > 0), , drop = FALSE]
  genus <- rownames(temp_a)
  for(gen in genus){
    indices <- grep(gen, rownames(temp_m))
    if(length(indices) > 0){
      temp_m$genus[indices] <- gen
      temp_m$metabolite[indices] <- substr(rownames(temp_m)[indices], nchar(gen) + 1, nchar(rownames(temp_m)[indices]))
    }
  }
  df_list[[i]] <- temp_m
}
df_depression <- do.call(rbind, df_list)
##########
taxon<-read.csv('D:/microbiome/micom/202404/paper/data/taxa.csv')
tsnedata<-read.csv('D:/microbiome/micom/202404/paper/result/1-1/tsne_data.csv',row.names = 1)
tsnedata$group<-meta$disease[match(tsnedata$sampleid,rownames(meta))]
#tsnedata<-tsnedata[tsnedata$genus%in%diff,]
tsnedata$phyla<-taxon$p__[match(tsnedata$genus,taxon$g__)]
ggplot(tsnedata, aes(x = TSNE.1, y = TSNE.2,color=genus, shape = group)) +
  geom_point(size = 1,alpha=0.6) +
  theme_minimal() +
  theme(legend.position = 'none')
##########
mes<-read.csv('D:/microbiome/micom/202404/paper/result/1-1/mes.csv',row.names = 1)
mes$group<-meta$disease[match(mes$sampleid,rownames(meta))]
mes_diff<-mes[mes$metabolite%in%diff_metabolite,]
mes_diff<-mes_diff[mes_diff$MES>0,]
mes<-mes[mes$MES>0,]
ggplot(mes_diff,aes(x=MES,fill=group))+
  geom_density(alpha=0.6)+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line())+
  scale_fill_manual(values=c('Depression'='#E7B800', 'Healthy'='#00AFBB'))+
  xlim(c(0,6))
  

metabolite_D_diff<-metabolite_D[diff_metabolite,]
metabolite_H_diff<-metabolite_H[diff_metabolite,]
rowMeans(metabolite_D_diff)/rowMeans(metabolite_H_diff)

###########
all_genus<-unique(tsnedata$genus)
res_tsne<-c()
for(i in 1:length(all_genus)){
  print(i)
  tsne1_d<-tsnedata_d$TSNE.1[which(tsnedata_d$genus==all_genus[i])]
  tsne1_h<-tsnedata_h$TSNE.1[which(tsnedata_h$genus==all_genus[i])]
  freq_d<-length(tsne1_d)
  freq_h<-length(tsne1_h)
  if(freq_d>0&freq_h>0){
    print(shapiro.test(tsne1_d)$p.val)
    print(shapiro.test(tsne1_h)$p.val)
    temp<-c(t.test(tsne1_d,tsne1_h)$p.val,freq_d,freq_h)%>%as.data.frame()
  }else{
    temp<-c(NA,freq_d,freq_h)%>%as.data.frame()
  }
  temp<-t(temp)%>%as.data.frame()
  colnames(temp)<-c('p.val','freq_d','freq_h')
  rownames(temp)<-all_genus[i]
  res_tsne<-rbind(res_tsne,temp)
}
res_tsne$p.adj<-NA
res_tsne$p.adj<-p.adjust(res_tsne$p.val,method = 'BH')
######
taxon<-read.csv('D:/microbiome/micom/202404/paper/data/taxa.csv')
meta_H<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/meta_H_used.csv',row.names = 1)
meta_D<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/meta_D_used.csv',row.names = 1)
colnames(meta_H)[5]<-'disease'
colnames(meta_D)[5]<-'disease'
meta_H$disease<-'Healthy'
meta_D$disease<-'Depression'
meta<-rbind(meta_D,meta_H)
abundance_D<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/depression_genus.csv',row.names = 1)
abundance_D<-abundance_D[match(rownames(meta_D),colnames(abundance_D))]
temp<-abundance_D
temp[temp!=0]<-1
#abundance_D<-abundance_D[which(rowSums(temp)>ncol(abundance_D)*0.1)%>%as.numeric(),]
abundance_H<-read.csv('D:/microbiome/micom/202404/paper/data/1-1/healthy_genus.csv',row.names = 1)
abundance_H<-abundance_H[match(rownames(meta_H),colnames(abundance_H))]
temp<-abundance_H
temp[temp!=0]<-1
#abundance_H<-abundance_H[which(rowSums(temp)>ncol(abundance_H)*0.1)%>%as.numeric(),]
abundance<-merge(abundance_D,abundance_H,by = "row.names",all = TRUE)
abundance[is.na(abundance)]<-0
rownames(abundance)<-abundance$Row.names
abundance$Row.names<-NULL
abundance$phyla<-taxon$p__[match(rownames(abundance),taxon$g__)]
abundance<-na.omit(abundance)
Firm<-colSums(abundance[abundance$phyla=='Firmicutes',1:ncol(abundance)-1])
Bact<-colSums(abundance[abundance$phyla=='Bacteroidetes',1:ncol(abundance)-1])
F_B<-(Firm/Bact)%>%as.data.frame()
colnames(F_B)<-'Firmicutes/Bacteroidetes'
F_B$group<-meta$disease[match(rownames(F_B),rownames(meta))]
ggplot(F_B,aes(x=group,y=`Firmicutes/Bacteroidetes`,))+
  geom_boxplot()+
  ylim(c(0,10))
abundance$phyla<-NULL
shannon<-diversity(t(abundance))%>%as.data.frame()
simpson<-diversity(t(abundance),index = 'simpson')%>%as.data.frame()
abundance[abundance<0.01]<-0
temp<-abundance
temp$phyla<-NULL
temp[temp!=0]<-1
temp1<-rowSums(temp)
abundance<-abundance[rownames(temp)[temp1>53],]
abundance1<-rbind(abundance,shannon,simpson)
a<-colSums(abundance)

b<-colSums(abundance[,1:ncol(abundance)-1])
(a-b)%>%sort()%>%plot()
(a-b)%>%sort()%>%hist()

#t.test(F_B$`Firmicutes/Bacteroidetes`[F_B$group=='Depression'],
#       F_B$`Firmicutes/Bacteroidetes`[F_B$group=='Healthy'],paired = T)
#####
abundance1<-t(abundance1)%>%as.data.frame()
abundance1$`F/B`<-F_B$`Firmicutes/Bacteroidetes`
exchange_mat<-read.csv('D:/xy/20250108/exchanges_mat.csv',row.names = 1)
topo<-read.csv('D:/microbiome/Depression/topolog.csv',row.names = 1)
id<-intersect(rownames(abundance1),rownames(topo))
abundance1<-abundance1[id,]
exchange_mat<-exchange_mat[id,]
topo<-topo[id,]
shannon<-shannon[id,]
simpson<-simpson[id,]
data<-cbind(#abundance1,
            exchange_mat,
            topo
            )
#data<-cbind(shannon,simpson,exchange_mat,topo)

data$group<-meta$disease[match(rownames(data),rownames(meta))]
data$Disease.name<-NULL
data_split <- createDataPartition(data$group, p = .70, list = FALSE)
training_data <- data[ data_split,]
testing_data  <- data[-data_split,]
train_index <- createFolds(training_data$group, k = 10)

randomForestFit <- training_data %>% train(group ~ .,
                                           method = "rf",
                                           data = .,
                                           tuneLength = 5,
                                           trControl = trainControl(method = "cv", indexOut = train_index))
randomForestFit
pr_rf<-predict(randomForestFit,testing_data)
confusionMatrix(pr_rf,reference = testing_data$group%>%as.factor())
pr_rf<-predict(randomForestFit,testing_data,type = 'prob')
library(pROC)
roc_curve<-roc(testing_data$group,pr_rf[,1])
plot(roc_curve)
sp.obj <- ci.sp(roc_curve, sensitivities=seq(0, 1, .01), boot.n=100) 
plot(sp.obj, type="shape", col="#5bd1d7")
text(0,0.5,paste0("AUC:",round(roc_curve[["auc"]],4)))
varImp(randomForestFit)
