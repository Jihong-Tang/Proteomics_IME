#### Visualization codes for Figure 8
#### AI-aided multi-omics diagnostic utility for the IDHm-IME gliomas
#### Author: Jihong TANG; Jiguang WANG

# Fig8A DNA MNP projection ----
library(tidyverse)
methyl.meta <- read_delim("./data/metadata_new_version.txt")

table(methyl.meta$tsne_Finallabel_DNAmethylation_class)
table(methyl.meta$WHO2021_group)

methyl.cgpa <- read.csv("./data/cgga.MutA.methyl.1029_axis.csv")
colnames(methyl.cgpa)[1] <- "Cohort_ID"
nmf.cgpa <- read_delim("../Figure1_Protein_Clusters/data/nmf_subtypes_clinical_updated_23Dec03.xlsx")

axis.cgpa <- merge(methyl.cgpa, nmf.cgpa, by = "Cohort_ID")
Pth_color <- c()
methyl.meta$fill <- ifelse(methyl.meta$WHO2021_group == "Adult-type diffuse gliomas", methyl.meta$tsne_Finallabel_DNAmethylation_class, "Other")
axis.cgpa <- axis.cgpa[order(axis.cgpa$NMF_Cluster, decreasing = T), ] %>% filter(Cohort_ID != "CGGA_P84")


methyl.cgpa.WT <- read.csv("./data/cgga.WT.methyl.1211.axis.csv") %>% filter(X %in% c("CGGA_P712", "CGGA_P924", "CGGA_P836", "CGGA_P1152"))

ggplot() +
  geom_point(data = methyl.meta, aes(x = TsneP1, y = TsenP2, color = fill), shape = 16, size = 1, alpha = .7) +
  geom_point(data = axis.cgpa, aes(x = axis1, y = axis2, fill = NMF_Cluster), shape = 24, size = 2) +
  scale_fill_manual(
    values = c(
      "ProteoNMF1" = "#33a02c", "ProteoNMF2" = "#e31a1c",
      "ProteoNMF3" = "#ff7f00", "ProteoNMF4" = "#1f78b4"
    ),
    name = "Proteomics subgroup"
  ) +
  scale_color_manual(
    values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#ffd92f", "#b3b3b3"),
    name = "CNS reference"
  ) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
# ggsave("figures/Fig3/1217_Fig3A_cgga_methyl_reference_mapping.pdf", width = 8, height = 6)
ggsave("figures/Fig8A_cgga_methyl_reference_mapping.pdf", width = 8, height = 6)

# Fig8B DNA methylation heatmap ---- 
library(tidyverse)

library(tidyverse)
cl_3C_IDHmut <- read_delim("..//data/nmf_subtypes_clinical_updated_23Dec03.xlsx")

df.methyl.sigs <- read_delim("../../241009_MultiOmics_Classifier_forIDHmes/results/1008_ProteoNMF_M450K_signatures_W.txt")

df_CGGA_methyl <- readRDS("../../241009_MultiOmics_Classifier_forIDHmes/data/CGGA_Methylation_Data_103samples_Filtered_1225.rds")
df_MutA_methyl <- df_CGGA_methyl[, colnames(df_CGGA_methyl) %in% cl_3C_IDHmut$Cohort_ID]
load("../../241009_MultiOmics_Classifier_forIDHmes/results/1009_res_DE_ProteoNMF_methylation_signatures.Rdata")

top <- 750

methyl_Proteo1_W <- res_DE_ProteoNMF1_methyl$pid[order(res_DE_ProteoNMF1_methyl$fc, decreasing = T)][1:top]
methyl_Proteo2_W <- res_DE_ProteoNMF2_methyl$pid[order(res_DE_ProteoNMF2_methyl$fc, decreasing = T)][1:top]
methyl_Proteo3_W <- res_DE_ProteoNMF3_methyl$pid[order(res_DE_ProteoNMF3_methyl$fc, decreasing = T)][1:top]
methyl_Proteo4_W <- res_DE_ProteoNMF4_methyl$pid[order(res_DE_ProteoNMF4_methyl$fc, decreasing = T)][1:top]

methyl_Proteo1_W <- res_DE_ProteoNMF1_methyl$pid[order(res_DE_ProteoNMF1_methyl$tstat, decreasing = T)][1:top]
methyl_Proteo2_W <- res_DE_ProteoNMF2_methyl$pid[order(res_DE_ProteoNMF2_methyl$tstat, decreasing = T)][1:top]
methyl_Proteo3_W <- res_DE_ProteoNMF3_methyl$pid[order(res_DE_ProteoNMF3_methyl$tstat, decreasing = T)][1:top]
methyl_Proteo4_W <- res_DE_ProteoNMF4_methyl$pid[order(res_DE_ProteoNMF4_methyl$tstat, decreasing = T)][1:top]

#list_cpgs <- c(df.methyl.sigs$ProteoNMF1, df.methyl.sigs$ProteoNMF2, df.methyl.sigs$ProteoNMF3, df.methyl.sigs$ProteoNMF4)
list_cpgs <- c(methyl_Proteo1_W, methyl_Proteo2_W, methyl_Proteo3_W, methyl_Proteo4_W)

plt_mutA <- df_MutA_methyl[list_cpgs, ]
plt_cl <- cl_3C_IDHmut[cl_3C_IDHmut$Cohort_ID %in% colnames(plt_mutA), c("Cohort_ID", "NMF_Cluster")]

plt_mutA <- plt_mutA[, plt_cl$Cohort_ID]
plt_mutA_scale <- apply(plt_mutA, 1, scale) %>% t()
colnames(plt_mutA_scale) <- colnames(plt_mutA)
sample_order_MutA <- plt_cl$Cohort_ID[order(plt_cl$NMF_Cluster, decreasing = F)]

library(ComplexHeatmap)
library(dendextend)
library(circlize)
col_fun = colorRamp2(seq(-1, 1, 0.1), colorRampPalette(c('#29419E','#698DC9',"grey95",'#F67172','#DC2B18'))(21))
top_anno <- HeatmapAnnotation(NMF_Cluster = plt_cl$NMF_Cluster,
                              
                              col = list(
                                NMF_Cluster = c("ProteoNMF1" = "#33a02c", "ProteoNMF2" = '#e31a1c', "ProteoNMF3" = '#ff7f00', "ProteoNMF4" = '#1f78b4')
                              ),
                              annotation_label = c("Proteomics subgroup"), 
                              annotation_name_gp= gpar(fontsize = 8),
                              simple_anno_size = unit(.25, "cm"),
                              annotation_legend_param = list(direction = "horizontal"))

ht_methyl <- Heatmap(as.matrix(plt_mutA_scale), 
                 name = "DNA methylation", 
                 height = unit(4, "cm"), width = unit(6, "cm"),
                 column_split = plt_cl$NMF_Cluster, column_order = sample_order_MutA,
                 show_row_names = F, show_row_dend = F, top_annotation = top_anno,
                 show_column_names = F,
                 row_names_gp = gpar(fontsize = 10),
                 column_names_gp = gpar(fontsize = 10),
                 col = col_fun, 
                 cluster_rows = F, cluster_columns = F, column_title =' ',
                 use_raster = F
)                             
ht_methyl                             
pdf("./figures/Fig8B_methyl_heatmap.pdf", width = 5, height = 4)
draw(ht_methyl, heatmap_legend_side = "right", 
     annotation_legend_side = "bottom")
dev.off()

### Fig8D Balanced Accuracy ----
library(tidyverse)
per_all <- read_delim("data/1021_performance_all.txt", delim = "\t")

plt_acc <- per_all[, c("Balanced Accuracy", "omics")] %>% as.data.frame()

ggplot(plt_acc, aes(x = omics, y = `Balanced Accuracy`)) +
  geom_bar(stat = "identity", fill = "#a6cee3", width = 0.7) +
  # geom_text(aes(label = round(`Balanced Accuracy`, 3)), vjust = -0.3) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete(limits = c("I", "T", "IT", "M", "IM", "TM", "ITM", "P")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Balanced Accuracy", x = "Omics", y = "Balanced Accuracy") +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "transparent", color = "transparent"), plot.margin = unit(c(2, 2, 2, 2), "lines"),
    plot.title = element_text(size = 34, vjust = 0.5, hjust = 0.5, face = "bold.italic", color = "transparent"), text = element_text(size = 14, face = "bold"),
    legend.key.width = unit(0.6, "cm"), legend.key.height = unit(0.6, "cm"), legend.position = "none", legend.text = element_text(size = 14, hjust = 0, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"), axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_text(size = 14, face = "plain", color = "black"), axis.title.y = element_text(size = 14, face = "plain", color = "black")
  )
ggsave("figures/Fig8D_Balanced_Accuracy.pdf", width = 5, height = 4, dpi = 300)

### Fig8E F1 score ----
plt_F1 <- per_all[, c("F1", "omics")] %>% as.data.frame()

ggplot(plt_F1, aes(x = omics, y = F1)) +
  geom_bar(stat = "identity", fill = "#fb9a99", width = 0.7) +
  # geom_text(aes(label = round(F1, 3)), vjust = -0.3) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete(limits = c("I", "T", "IT", "M", "IM", "TM", "ITM", "P")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Fa", x = "Omics", y = "F1") +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "transparent", color = "transparent"), plot.margin = unit(c(2, 2, 2, 2), "lines"),
    plot.title = element_text(size = 34, vjust = 0.5, hjust = 0.5, face = "bold.italic", color = "transparent"), text = element_text(size = 14, face = "bold"),
    legend.key.width = unit(0.6, "cm"), legend.key.height = unit(0.6, "cm"), legend.position = "none", legend.text = element_text(size = 14, hjust = 0, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"), axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_text(size = 14, face = "plain", color = "black"), axis.title.y = element_text(size = 14, face = "plain", color = "black")
  )
ggsave("figures/Fig8E_F1score.pdf", width = 5, height = 4, dpi = 300)
