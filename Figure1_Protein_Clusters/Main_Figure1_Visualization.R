#### Visualization codes for Main Figure 1; 2024-04-04
#### Protein clusters of IDH-mutant astrocytoma with divergent survival outcomes
#### Author: Jihong TANG; Jiguang WANG

# Fig 1A - Schematic workflow for data cohorts ----
## All the figures were manually created using PowerPoint software.

# Fig 1B - NMF and genomics heatmap ----
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
pdf("./figures/Fig1B_Ht_cl&pn.pdf", width = 10, height = 8)
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

pdf("./figures/Fig1B_Ht_genomics.pdf", width = 10, height = 8)
draw(ht_Genomics_MutA_long, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
dev.off()

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
