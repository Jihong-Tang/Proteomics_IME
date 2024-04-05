#### Visualization codes for Main Figure 1; 2024-04-04
#### Proteomics subgroups of IDHmut-noncodel gliomas with
#### divergent survival outcomes

# Fig 1A - NMF and genomics heatmap ----
## 1A-1 heatmap of clinical & pns ----
pth_cl <- "data/nmf_subtypes_clinical_updated_23Dec03.xlsx"
pth_ht <- "data/Heatmap_Input_MutA_Unbiased_4Groups_Aug18.rds"
ht_cl <- readxl::read_xlsx(pth_cl)
ht_input <- readRDS(pth_ht)
ht_input <- ht_input[, ht_cl$Cohort_ID]

ht_cl$Age <- as.numeric(ht_cl$Age)
age <- ht_cl$Age
ht_cl$Age_range <- ifelse(age<25, "[10, 25)",
                          ifelse(age <40, "[25, 40)",
                                 ifelse(age <55, "[40, 55)", "[55, 70]")))

sample_order_MutA <- ht_cl$Cohort_ID[order(ht_cl$Cohorts, decreasing = F)]

library(ComplexHeatmap)
library(dendextend)
library(circlize)
col_fun = colorRamp2(seq(-1.5, 1.5, 0.1), colorRampPalette(c('#29419E','#698DC9',"grey95",'#F67172','#DC2B18'))(31))
nmf_color <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')

top_anno <- HeatmapAnnotation(NMF_Cluster = ht_cl$NMF_Cluster,
                              Cohorts = ht_cl$Cohorts,
                              Grade_2016 = ht_cl$Grade_2016,
                              Grade_2021 = ht_cl$Grade_2021,
                              Race = ht_cl$Race,
                              Gender = ht_cl$Gender,
                              Age = ht_cl$Age_range,
                              
                              col = list(
                                NMF_Cluster = c("ProteoNMF1" = "#33a02c", "ProteoNMF2" = '#e31a1c', "ProteoNMF3" = '#ff7f00', "ProteoNMF4" = '#1f78b4'),
                                Cohorts = c("CGPA" = "#bc80bd", "CPTAC" = "#386cb0", "CRM2022_Jakob" = "#8dd3c7"), # nolint
                                Grade_2016 = c("WHO II" = "#a1d99b", "WHO III" = "#41ab5d", "WHO IV" = '#006d2c', "NA" = "grey"),
                                Grade_2021 = c("G2" = "#a1d99b", "G3" = "#41ab5d", "G4" = '#006d2c', "NA" = "grey"),
                                Race = c("Asian" = "#b15928", "Caucasian" = "#ffff99", "NA" = "grey"),
                                Gender = c("Male" = "#377eb8", "Female" = "#f781bf", "NA" = "grey"), 
                                Age = c("[10, 25)" = '#cbc9e2', "[25, 40)" = '#9e9ac8', "[40, 55)" = '#756bb1', "[55, 70]" = '#54278f')
                              ),
                              annotation_label = c("Proteomics subgroup", "Data cohorts", "Grade WHO2016", "Grade WHO2021", 
                                                   "Race", "Gender", "Age"), 
                              annotation_name_gp= gpar(fontsize = 8),
                              simple_anno_size = unit(.25, "cm"),
                              annotation_legend_param = list(direction = "horizontal"))

ht_pn <- Heatmap(as.matrix(ht_input), 
                 name = "Protein Relative Abundance", 
                 height = unit(8, "cm"), width = unit(13, "cm"),
                 column_split = ht_cl$NMF_Cluster, column_order = sample_order_MutA,
                 show_row_names = F, show_row_dend = F, top_annotation = top_anno,
                 show_column_names = F,
                 row_names_gp = gpar(fontsize = 10),
                 column_names_gp = gpar(fontsize = 10),
                 col = col_fun, 
                 cluster_rows = F, cluster_columns = F, column_title =' '
                 )                             
ht_pn                             
pdf("./figures/1203_Ht_cl&pn.pdf", width = 10, height = 8)
draw(ht_pn, heatmap_legend_side = "right", 
     annotation_legend_side = "bottom")
dev.off()

## 1A-2 heatmap of genomics information ----
pth_genomics <- "data/nmf_subtypes_genomics_updated_23Dec03.xlsx"
ht_genomics <- readxl::read_xlsx() %>% as.data.frame()
mut_gene <- c("IDH1", "TP53", "ATRX",
              "NF1","PIK3CA", "PIK3R1")

row.names(ht_genomics) <- ht_genomics$Cohort_ID
df_mut_MutA <- ht_genomics[sample_order_MutA, ]
df_mut_MutA[is.na(df_mut_MutA)] <- "NA"
mat_mut_MutA <- as.matrix(t(df_mut_MutA[, -seq(1,3)]))
mat_mut_MutA <- mat_mut_MutA[, ht_cl$Cohort_ID]

col_mut <- c("NA" = "#DDDDDD", No = "white",
             "Missense"= "#b2df8a", "Splice_Site" = "#a6cee3", 
             "Frameshift"= "#fb9a99", "Nonsense" = "#cab2d6", "Inframe_Indel" = "#fdbf6f"
)

col_mut <- c("NA" = "#DDDDDD", No = "white",
             "Missense"= "#74c476", "Splice_Site" = "#fd8d3c", 
             "Frameshift"= "#fb6a4a", "Nonsense" = "#9e9ac8", "Inframe_Indel" = "#6baed6"
)

top_simple <- HeatmapAnnotation(NMF_Cluster = ht_cl$NMF_Cluster,
                                
                                col = list(NMF_Cluster = c("ProteoNMF1" = "#33a02c", "ProteoNMF2" = '#e31a1c', "ProteoNMF3" = '#ff7f00', "ProteoNMF4" = '#1f78b4')
                                           ),
                                annotation_label = c("Proteomics subgroup"), 
                                annotation_name_gp= gpar(fontsize = 8),
                                simple_anno_size = unit(.25, "cm"),
                                annotation_legend_param = list(direction = "horizontal"))

ht_mut<- oncoPrint(mat_mut_MutA[mut_gene, ],
                              alter_fun = function(x, y, w, h, v) {
                                n = sum(v)  
                                h = h*0.9
                                if(n) grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.8, 1/n*h, 
                                                gp = gpar(fill = col_mut[names(which(v))], col = NA), just = "top")
                              }
                              , col = col_mut,
                              heatmap_legend_param = list(title = "Mutations",
                                                          at = c("Missense", "Splice_Site", "Frameshift", 
                                                                 "Nonsense", "Inframe_Indel", "No", "NA"),
                                                          labels = c("Missense", "Splice Site", "Frameshift", 
                                                                     "Nonsense", "Inframe Indel", "No", "Not Available")),
                              show_pct = F, row_names_side = "right", row_names_gp = gpar(fontface = "italic", fontsize = 10),
                              show_column_names = T, column_split = ht_cl$NMF_Cluster,
                              column_title = NULL,
                              row_order = mut_gene,
                              column_order = sample_order_MutA,
                              height = unit(0.4*length(mut_gene), "cm"), width = unit(13, "cm"),
                              right_annotation = NULL, top_annotation = top_simple,
                              border = T
)
ht_mut

cnv_gene <- c("CDKN2A", "PDGFRA", "CDK4", "CDK6",
              "CCND2", "MET", "MYC")

df_CNV_MutA <- ht_genomics[sample_order_MutA, cnv_gene]
df_CNV_MutA[is.na(df_CNV_MutA)] <- "NA"
mat_CNV_MutA <- as.matrix(t(df_CNV_MutA))
mat_CNV_MutA <- mat_CNV_MutA[, ht_cl$Cohort_ID]

col_CNV <- c("NA" = "#DDDDDD", "Neural" = "white", "Gain" = "#F67172", 
             "Amp" = "#DC2B18", "Loss" = "#698DC9", "Del" = "#29419E")

ht_CNV <- oncoPrint(mat_CNV_MutA[cnv_gene, ],
                         alter_fun = list(
                           #background = function(...) NULL,
                           background = function(x, y, w, h) grid.rect(x, y, w, h, 
                                                                       gp = gpar(fill = "white", col = "#f0f0f0")),
                           Neural = function(x, y, w, h) grid.rect(x, y, w*.6, h*.7, 
                                                                   gp = gpar(fill = col_CNV["Neural"], col = NA)),
                           Gain = function(x, y, w, h) grid.rect(x, y, w*.8, h*.9, 
                                                                 gp = gpar(fill = col_CNV["Gain"], col = NA)),
                           Amp = function(x, y, w, h) grid.rect(x, y, w*.8, h*.9, 
                                                                gp = gpar(fill = col_CNV["Amp"], col = NA)),
                           Loss = function(x, y, w, h) grid.rect(x, y, w*.8, h*.9, 
                                                                 gp = gpar(fill = col_CNV["Loss"], col = NA)),
                           Del = function(x, y, w, h) grid.rect(x, y, w*.8, h*.9, 
                                                                gp = gpar(fill = col_CNV["Del"], col = NA)),
                           "NA" = function(x, y, w, h) grid.rect(x, y, w*.8, h*.9, 
                                                                 gp = gpar(fill = col_CNV["NA"], col = NA))
                         ), col = col_CNV,
                         heatmap_legend_param = list(title = "CNV Classification",
                                                     at = c("Amp", "Gain", "Neural", "Loss", "Del", "NA"),
                                                     labels = c("Amp", "Gain", "Neural", "Loss", "Del","Not available")),
                         show_pct = F, row_names_side = "right", row_names_gp = gpar(fontface = "italic", fontsize = 10),
                         show_column_names = F, column_split = ht_cl$NMF_Cluster,
                         column_title = " ", 
                         row_order = cnv_gene,
                         column_order = sample_order_MutA,
                         right_annotation = NULL, top_annotation = NULL,
                         height = unit(0.4*length(cnv_gene), "cm"), width = unit(13, "cm"),
                         border = T 
)
ht_CNV
ht_Genomics_MutA_long <- ht_mut %v%  ht_CNV
ht_Genomics_MutA_long

pdf("./figures/1203_Ht_genomics.pdf", width = 10, height = 8)
draw(ht_Genomics_MutA_long, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
dev.off()

# Fig 1B - Clinical and genomics association ---- 
## 1B-1 prepare the correlation data ---- 
library(tidyverse)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

ht_genomics <- ht_genomics[ht_cl$Cohort_ID, ]
sum(ht_genomics$Cohort_ID == ht_cl$Cohort_ID)

new1 <- ht_cl
ProteoNMF1<-rep(0,length(new1$Cohort_ID))
ProteoNMF1[which(new1$NMF_Cluster=="ProteoNMF1")]=1
ProteoNMF2<-rep(0,length(new1$Cohort_ID))
ProteoNMF2[which(new1$NMF_Cluster=="ProteoNMF2")]=1
ProteoNMF3<-rep(0,length(new1$Cohort_ID))
ProteoNMF3[which(new1$NMF_Cluster=="ProteoNMF3")]=1
ProteoNMF4<-rep(0,length(new1$Cohort_ID))
ProteoNMF4[which(new1$NMF_Cluster=="ProteoNMF4")]=1

df_nmf <- data.frame(ProteoNMF1, ProteoNMF2, ProteoNMF3,
                          ProteoNMF4)

new2 <- ht_cl %>% as.data.frame()
Cohort.CGPA <- rep(NA, length(new2[, 1]))
Cohort.CGPA[which(new2$Cohorts == "CGPA")] = 1
Cohort.CGPA[which(new2$Cohorts != "CGPA")] = 0

Cohort.CPTAC <- rep(NA, length(new2[, 1]))
Cohort.CPTAC[which(new2$Cohorts == "CPTAC")] = 1
Cohort.CPTAC[which(new2$Cohorts != "CPTAC")] = 0

Cohort.CRM2022 <- rep(NA, length(new2[, 1]))
Cohort.CRM2022[which(new2$Cohorts == "CRM2022_Jakob")] = 1
Cohort.CRM2022[which(new2$Cohorts != "CRM2022_Jakob")] = 0

Grade.II <- rep(NA, length(new2[, 1]))
Grade.II[which(new2$Grade_2021 == "G2")] = 1
Grade.II[which(new2$Grade_2021 != "G2")] = 0

Grade.III <- rep(NA, length(new2[, 1]))
Grade.III[which(new2$Grade_2021 == "G3")] = 1
Grade.III[which(new2$Grade_2021 != "G3")] = 0

Grade.IV <- rep(NA, length(new2[, 1]))
Grade.IV[which(new2$Grade_2021 == "G4")] = 1
Grade.IV[which(new2$Grade_2021 != "G4")] = 0

Male <- rep(NA, length(new2[, 1]))
Male[which(new2$Gender == "Male")] = 1
Male[which(new2$Gender != "Male")] = 0

Female <- rep(NA, length(new2[, 1]))
Female[which(new2$Gender == "Female")] = 1
Female[which(new2$Gender != "Female")] = 0

Asian <- rep(NA, length(new2[, 1]))
Asian[which(new2$Race == "Asian")] = 1
Asian[which(new2$Race != "Asian")] = 0

Caucasian<- rep(NA, length(new2[, 1]))
Caucasian[which(new2$Race == "Caucasian")] = 1
Caucasian[which(new2$Race != "Caucasian")] = 0

df_cli <- data.frame(
  Cohort.CGPA, Cohort.CPTAC, Cohort.CRM2022, Grade.II, Grade.III, Grade.IV, 
  Male, Female, Asian, Caucasian)


new3 <- ht_genomics %>% as.data.frame()

#df_nmf_mut <- ht_genomics[, mut_gene] %>% as.data.frame()
IDH1 <- rep(NA,length(new3[,1]))
IDH1[which(new3$IDH1 == "No")] = 0
IDH1[which(new3$IDH1 != "No" & new3$IDH1 != "NA")] = 1

TP53 <- rep(NA,length(new3[,1]))
TP53[which(new3$TP53 == "No")] = 0
TP53[which(new3$TP53 != "No" & new3$TP53 != "NA")] = 1

ATRX <- rep(NA,length(new3[,1]))
ATRX[which(new3$ATRX == "No")] = 0
ATRX[which(new3$ATRX != "No" & new3$ATRX != "NA")] = 1

NF1 <- rep(NA,length(new3[,1]))
NF1[which(new3$NF1 == "No")] = 0
NF1[which(new3$NF1 != "No" & new3$NF1 != "NA")] = 1

PIK3CA <- rep(NA,length(new3[,1]))
PIK3CA[which(new3$PIK3CA == "No")] = 0
PIK3CA[which(new3$PIK3CA != "No" & new3$PIK3CA != "NA")] = 1

PIK3R1 <- rep(NA,length(new3[,1]))
PIK3R1[which(new3$PIK3R1 == "No")] = 0
PIK3R1[which(new3$PIK3R1 != "No" & new3$PIK3R1 != "NA")] = 1


CDKN2A.Loss <- rep(NA,length(new3[,1]))
CDKN2A.Loss[which(new3$CDKN2A == "Loss" | new3$CDKN2A == "Del")] = 1
CDKN2A.Loss[which(new3$CDKN2A == "Neural" | new3$CDKN2A == "Gain" | new3$CDKN2A == "Amp")] = 0

PDGFRA.Gain <- rep(NA,length(new3[,1]))
PDGFRA.Gain[which(new3$PDGFRA == "Gain" | new3$PDGFRA == "Amp")] = 1
PDGFRA.Gain[which(new3$PDGFRA == "Neural" | new3$PDGFRA == "Loss" | new3$PDGFRA == "Del")] = 0

CDK4.Gain <- rep(NA,length(new3[,1]))
CDK4.Gain[which(new3$CDK4 == "Gain" | new3$CDK4 == "Amp")] = 1
CDK4.Gain[which(new3$CDK4 == "Neural" | new3$CDK4 == "Loss" | new3$CDK4 == "Del")] = 0

CDK6.Gain <- rep(NA,length(new3[,1]))
CDK6.Gain[which(new3$CDK6 == "Gain" | new3$CDK6 == "Amp")] = 1
CDK6.Gain[which(new3$CDK6 == "Neural" | new3$CDK6 == "Loss" | new3$CDK6 == "Del")] = 0

CCND2.Gain <- rep(NA,length(new3[,1]))
CCND2.Gain[which(new3$CCND2 == "Gain" | new3$CCND2 == "Amp")] = 1
CCND2.Gain[which(new3$CCND2 == "Neural" | new3$CCND2 == "Loss" | new3$CCND2 == "Del")] = 0

MET.Gain <- rep(NA,length(new3[,1]))
MET.Gain[which(new3$MET == "Gain" | new3$MET == "Amp")] = 1
MET.Gain[which(new3$MET == "Neural" | new3$MET == "Loss" | new3$MET == "Del")] = 0

MYC.Gain <- rep(NA,length(new3[,1]))
MYC.Gain[which(new3$MYC == "Gain" | new3$MYC == "Amp")] = 1
MYC.Gain[which(new3$MYC == "Neural" | new3$MYC == "Loss" | new3$MYC == "Del")] = 0

df_genomics <- data.frame(IDH1, TP53, ATRX, NF1, PIK3CA, PIK3R1, 
                          CDKN2A.Loss,PDGFRA.Gain, CDK4.Gain, 
                          CDK6.Gain, CCND2.Gain, MET.Gain, MYC.Gain)

correlation_table_MutA <- cbind(df_nmf, df_cli, df_genomics)
muGene.P.MutA <- matrix(nrow=nrow(correlation_table_MutA),ncol=ncol(correlation_table_MutA))
muGene.P.MutA[correlation_table_MutA ==1 | correlation_table_MutA == "Yes"] <- 'Y'
muGene.P.MutA[correlation_table_MutA == 0 | correlation_table_MutA == 'No'] <- 'N'
colnames(muGene.P.MutA) <- colnames(correlation_table_MutA)

## 1B-2 calculate correlation between each feature ----
plot_data_MutA <- data.frame( rep(NA,length(correlation_table_MutA[1,])^2),rep(NA,length(correlation_table_MutA[1,])^2),rep(NA,length(correlation_table_MutA[1,])^2),rep(NA,length(correlation_table_MutA[1,])^2),rep(NA,length(correlation_table_MutA[1,])^2))
n<-1
#nfeaA <- 17##length(correlation_table_MutA[1,]) -1
nfeaA <- length(correlation_table_MutA[1,])-23

nfeaB <- length(correlation_table_MutA[1,])
for(i in 1:nfeaA){
  x<-i+1
  for(j in (nfeaA+1):nfeaB){
    plot_data_MutA[n,1] <- colnames(correlation_table_MutA)[i]
    plot_data_MutA[n,2] <- colnames(correlation_table_MutA)[j]
    
    Mu.FEtest <- cbind(c(0,0),c(0,0))
    Mu.FEtest[1,1] <- length(which(muGene.P.MutA[,i] == 'Y' & muGene.P.MutA[,j] == 'Y'))
    Mu.FEtest[1,2] <- length(which(muGene.P.MutA[,i] == 'Y' & muGene.P.MutA[,j] == 'N'))
    Mu.FEtest[2,1] <- length(which(muGene.P.MutA[,i] == 'N' & muGene.P.MutA[,j] == 'Y'))
    Mu.FEtest[2,2] <- length(which(muGene.P.MutA[,i] == 'N' & muGene.P.MutA[,j] == 'N'))
    coMutation <- (((Mu.FEtest[1,1])*(Mu.FEtest[2,2]))+1) / (((Mu.FEtest[1,2])*(Mu.FEtest[2,1]))+1) ## Odds ratio
    
    #fisher exact test
    pValue <- fisher.test(Mu.FEtest,alternative ="two.sided")$p.value
    
    if(pValue < 0.05){
      plot_data_MutA[n,3] <- -log10(pValue) #for dot size, but if pValue > 0.1 , filling with an instinct number
      if(coMutation > 1){
        plot_data_MutA[n,4] <- 'B_red'
      }
      else{
        plot_data_MutA[n,4] <- 'A_green'
      }
    }
    else{
      plot_data_MutA[n,3] <- 1#for dot size, but if pValue > 0.1 , filling with an instinct number
      plot_data_MutA[n,4] <- 'C_grey'
    }
    plot_data_MutA[n,5] <- log10(coMutation)
    n<-n+1
  }
}
plot_data_MutA<-plot_data_MutA[which(!(is.na(plot_data_MutA[,1]))),]
colnames(plot_data_MutA) <- c('FeatureA','FeatureB','dotSize','dotColor',"OR")
plot_data_MutA$orderID<-c(1:nrow(plot_data_MutA))

## 1B-3 plot co-mutation ----
library(tidyverse)
plot_2.data <- plot_data_MutA[which(plot_data_MutA$dotSize > 1),]
plot_3.data <- plot_data_MutA[which(plot_data_MutA$dotSize == 1),]
coMu.plot<-ggplot()+theme_classic()
coMu.plot<-coMu.plot+geom_vline(xintercept = unique(plot_data_MutA$FeatureA), linetype=2,color="grey",size=0.3)
coMu.plot<-coMu.plot+geom_hline(yintercept = unique(plot_data_MutA$FeatureB), linetype=2,color="grey",size=0.3)

coMu.plot<-coMu.plot+geom_point(data = plot_data_MutA,aes(y=reorder(FeatureB,-orderID),x=reorder(FeatureA,orderID),size=dotSize,fill=OR), shape = 21, color = "white")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
coMu.plot
coMu.plot<-coMu.plot+geom_point(data = plot_2.data,aes(y=reorder(FeatureB,-orderID),x=reorder(FeatureA,orderID),size=dotSize,fill=OR),shape=21)+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
coMu.plot
coMu.plot<-coMu.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(1,4,1,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                           text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(1,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="bottom",
                           legend.text=element_text(size=12,face='plain'),axis.text.y=element_text(size=16,vjust=0.5,hjust=1,face='plain',color='black'),legend.title=element_text(size=12,vjust=0.5,hjust=0,face='plain'),
                           axis.text.x=element_text(size=16,angle=45,vjust=1,hjust=1,face='plain',color='black'),axis.title.x=element_text(size=20,vjust=0,hjust=0.5,face='plain',color='black'),
                           axis.title.y=element_text(size=20,hjust=0.5,vjust=2,face='plain',color='black'),
                           strip.text = element_text(size=18,face='bold',vjust=0.5,hjust=0.5),strip.background = element_rect(colour="black", fill=gg_color_hue(3)))
coMu.plot
coMu.plot<-coMu.plot+scale_shape_manual(name=NULL,values=c(A_green=15,B_red=17,C_grey=16),
                                        labels=c(A_green='Mutual exclusion',B_red='Co-occurrence',C_grey='Not significant'),
                                        guide = guide_legend(override.aes=list(size=4),nrow=4),na.translate = F)
coMu.plot<-coMu.plot+scale_size(name   = "Fisher's test\n   p-value",  breaks = c(1, 1.3,2 , 4, 6,8),labels = expression(1,0.05,10^-2, 10^-4,10^-6,10^-8) ,guide = guide_legend(nrow=2))
coMu.plot<-coMu.plot+scale_fill_gradientn(name="        Odds ratio\n(log10 transformed)",colours=c('#276419','#4d9221','#7fbc41','white','#de77ae','#c51b7d','#8e0152'),limits=c(-2.5,2.5),breaks=seq(-2,2,1))
coMu.plot<-coMu.plot+scale_x_discrete(position = "bottom") 
coMu.plot
plotAll<-rbind(ggplotGrob(coMu.plot),size="first")
plotAll
ggsave(file="./figures/1203_correlation_cli_genomics_V2.pdf", plot=plotAll,bg = 'white', width = 15, height = 20, units = 'cm', dpi = 600)

coMu.plot.NoL <- coMu.plot+theme(axis.text.x  = element_blank(), 
                                 axis.text.y  = element_blank(),
                                 legend.position = "none")
coMu.plot.NoL
ggsave(file="./figures/1203_correlation_cli_genomics_noL.pdf", plot=coMu.plot.NoL,bg = 'white', width = 9, height = 10, units = 'cm', dpi = 600)


# Fig 1C - Survival analysis ---- 
ht_cl <- readxl::read_xlsx("data/nmf_subtypes_clinical_updated_23Dec03.xlsx")
library(survival)
library(survminer)
#cl_TME_4G_IDHmut_r4$Mixed <- ifelse(cl_TME_4G_IDHmut_r4$membership < 0.5, "Mixed", "Non-Mixed")
#cl_5k_clV2_noncodel_rank4
OS_cl <- ht_cl[, c("OS_day", "Censor_OS", "NMF_Cluster")]
colnames(OS_cl) <- c("OS", "Censor", "Expre")
OS_cl$OS <- as.numeric(OS_cl$OS)
OS_cl$Censor <- as.numeric(OS_cl$Censor)
Pth_color <- c('#33a02c','#e31a1c','#ff7f00','#1f78b4')

custom_theme <- function() {
  theme_classic() %+replace%
    theme(
      text = 	element_text(size = 14, face = 'bold'),
      plot.title=element_text(hjust=0.5)
    )
}
p <- ggsurvplot(
  fit = survfit(Surv(OS, Censor) ~ Expre, data = OS_cl), 
  size = .75, risk.table = F, pval = F, ggtheme = custom_theme(), 
  #title = paste0("IDH-Mut-Noncodel Overall Survival"),
  xlab = "Years from diagnosis", ylab = "Survival Probability",
  
  legend = "none", legend.title = "", font.legend = c(10, "bold", "black"), 
  
  break.x.by = 365.25, xscale = "d_y",break.y.by = 0.2, axes.offset = T,
  font.x = c(14, "plain", "black"), font.y = c(14, "plain", "black"),
  font.tickslab = c(12, "bold", "black"), 
  palette = Pth_color
)
p
ggsave(paste0("./figures/1211_OS_4groups.pdf"), width = 7.5, height =6, units = 'cm', dpi = 600)


PFS_cl <- ht_cl[, c("PFS_day", "Censor_PFS", "NMF_Cluster")]
colnames(PFS_cl) <- c("PFS", "Censor", "Expre")
PFS_cl$PFS <- as.numeric(PFS_cl$PFS)
PFS_cl$Censor <- as.numeric(PFS_cl$Censor)
Pth_color <- c('#33a02c','#e31a1c','#ff7f00','#1f78b4')

p <- ggsurvplot(
  fit = survfit(Surv(PFS, Censor) ~ Expre, data = PFS_cl), 
  size = .75, risk.table = F, pval = F, ggtheme = custom_theme(), 
  #title = paste0("IDH-Mut-Noncodel Overall Survival"),
  xlab = "Years from diagnosis", ylab = "Survival Probability",
  
  legend = "none", legend.title = "", font.legend = c(10, "bold", "black"), 
  
  break.x.by = 365.25/1, xscale = "d_y",break.y.by = 0.2, axes.offset = T,
  font.x = c(14, "plain", "black"), font.y = c(14, "plain", "black"),
  font.tickslab = c(12, "bold", "black"), 
  palette = Pth_color
)
p
ggsave(paste0("./figures/1211_PFS_4groups.pdf"), width = 7.5, height =6, units = 'cm', dpi = 600)
