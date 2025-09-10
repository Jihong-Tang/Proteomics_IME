#### Visualization codes for Figure 3
#### Single-cell analysis reveals distinct lymphocyte infiltration and enhanced immune activity in IDHm-IME gliomas
#### Author: Jihong TANG; Jiguang WANG
library(tidyverse)
library(Seurat)

# Fig3A DimPlot---- 
obj.use <- readRDS("./data/scRNA_revise/0316_Merged.Major.18MutA.QC.rds")
table(obj.use$orig.ident, obj.use$RNA_snn_res.2)

color_scheme2 <- c('Tumor' = '#ea3323', 'Oligo' = '#8ab0d0', 'Myeloid' = '#b0d667',
                   'Lymph' = '#ed926b')
order <- c("Tumor", "Oligo", "Myeloid", "Lymph")

order <- rev(order)
DimPlot(object = obj.use, 
        group.by = "cell_type", 
        cols = color_scheme2,
        label = T, 
        pt.size = .1, #order = order,
        label.size = 0) & NoAxes()
ggsave("./figures/Fig3A_sc_major.DimPlot.pt0.1.pdf", width = 7, height = 6)

genes_to_plot = c("CD3D", "CD79A", "NKG7", "MAG", "CSF1R", "PTPRZ1", "SOX2", "BCAN", "TENM1", "TENM4")
genes_to_plot = rev(genes_to_plot)

tmp = DotPlot(obj.use, 
              features = genes_to_plot,
              cols = "RdBu",
              group.by = "cell_type")

tmp = tmp + 
  scale_y_discrete(limits = c("Lymph", "Oligo", "Myeloid", "Tumor")) + 
  coord_flip() + 
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1),
        axis.text.y = element_text(size=14))

tmp
ggsave("figures/DotPlot.Major.Cluster.pdf", width = 4.5, height = 4)

# Fig 3B - downsampling four samples ---- 
library(Seurat)
library(tidyverse)
obj.use <- readRDS("./data/scRNA_revise/0316_Merged.Major.18MutA.QC.rds")

xxx <- read_delim("data/scRNA_cluster.txt")
#obj.use$bulkMES <- ifelse(obj.use$orig.ident %in% c("P698", "P714", "P737"), "MES", "NotMES")
obj.use$bulkProt <- ifelse(obj.use$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "MIX"], "MIX", 
                           ifelse(obj.use$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "PPR"], "PPR",
                                  ifelse(obj.use$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "IME"], "IME", "NEU")))
DimPlot(obj.use, group.by = "bulkProt", order = c( "IME","NEU", "PPR", "MIX"), 
        pt.size = .5, cols = c("IME" = "#ff7f00", "NEU" = "#1f78b4", "PPR" =  "#e31a1c", "MIX" = "#33a02c")
) + NoAxes()

set.seed(123)
ds <- rep("No", length(obj.use$orig.ident))
idx <- sample(1:length(obj.use$orig.ident), length(obj.use$orig.ident)/4)
ds[idx] <- "Yes"

obj.use$downsamplig <- ds

obj.ds <- subset(obj.use, downsamplig == "Yes")
DimPlot(obj.ds, group.by = "bulkProt", order = c( "IME", "MIX","PPR", "NEU" ), 
        pt.size = .1, cols = c("IME" = "#ff7f00", "NEU" = "#1f78b4", "PPR" =  "#e31a1c", "MIX" = "#33a02c")
) + NoAxes()

ggsave("figures/Fig3B_downsampling.pt0.1.pdf", width = 7.5, height = 6)

# Fig 3C - major fraction comparison ---- 
count.all <- as.data.frame(table(obj.use$cell_type, obj.use$orig.ident))
count.wide <- reshape2::dcast(count.all, Var2 ~ Var1)
rownames(count.wide) <- count.wide$Var2
frac.all <- apply(count.wide[, -1], 1, function(x)x/sum(x))
frac.long <- reshape2::melt(frac.all)
colnames(frac.long) <- c("Var2", "Cohort_ID", "frac")

xxx <- read_delim("data/scRNA_cluster.txt")
colnames(frac.all)[1] <- "Cohort_ID"
frac.long <- merge(frac.long, xxx, by = "Cohort_ID")


library(ggbeeswarm)
NAN_plot <- ggplot(data=frac.long,aes(x=NMF_Cluster, y=frac)) + theme_classic() 
#NAN_plot
NAN_plot <- NAN_plot + 
  geom_boxplot(data=frac.long,aes(x = factor(NMF_Cluster, levels = c("MIX", "PPR", "IME", "NEU")), y = frac),width = 0.6,size=0.5,fill="transparent", outlier.color = "white")+
  geom_quasirandom(data=frac.long,aes(x = factor(NMF_Cluster, levels = c("MIX", "PPR", "IME", "NEU")), y = frac, color = NMF_Cluster, fill = NMF_Cluster),width = 0.25,size=0.75,alpha=.6,stroke=0.8, varwidth = T) + 
  facet_grid(.~factor(Var2)) 
NAN_plot
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks = seq(-10,20,0.2)) 
NAN_plot<- NAN_plot +ylab("Fraction") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(c("IME" = "#ff7f00", "NEU" = "#1f78b4", "PPR" =  "#e31a1c", "MIX" = "#33a02c")))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(c("IME" = "#ff7f00", "NEU" = "#1f78b4", "PPR" =  "#e31a1c", "MIX" = "#33a02c")))

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
NAN_plot
ggsave(file = "./figures/Fig3C_boxplot_compare_major.pdf", plot = figure_2, bg = "white", width = 13.5, height = 9, units = "cm", dpi = 600)

library(tidyverse)
library(ggpubr)
frac.long$compare <- ifelse(frac.long$NMF_Cluster == "IME", "Yes", "No")
compare_means(method = "wilcox.test", frac ~ compare,   
              data = frac.long %>% filter(Var2 == "Lymph") ,paired = F )
compare_means(method = "kruskal.test", frac ~ NMF_Cluster,   
              data = frac.long %>% filter(Var2 == "Lymph"),paired = F )
compare_means(method = "kruskal.test", frac ~ NMF_Cluster,   
              data = frac.long %>% filter(Var2 == "Myeloid"),paired = F )
compare_means(method = "kruskal.test", frac ~ NMF_Cluster,  
              data = frac.long %>% filter(Var2 == "Oligo"),paired = F )
compare_means(method = "kruskal.test", frac ~ NMF_Cluster,  
              data = frac.long %>% filter(Var2 == "Tumor"),paired = F )

compare_means(method = "anova", frac ~ NMF_Cluster,   
              data = frac.long %>% filter(Var2 == "Lymph") ,paired = F )
compare_means(method = "anova", frac ~ NMF_Cluster,   
              data = frac.long %>% filter(Var2 == "Myeloid"),paired = F )
compare_means(method = "anova", frac ~ NMF_Cluster,  
              data = frac.long %>% filter(Var2 == "Oligo"),paired = F )
compare_means(method = "anova", frac ~ NMF_Cluster,  
              data = frac.long %>% filter(Var2 == "Tumor"),paired = F )

rm(obj.use)
# Fig 3D tumor Dimplot ----
obj.tumor <- readRDS("data/scRNA_revise/0316_Merged.18MutA.tumor.rds")
DimPlot(obj.tumor, group.by = "cluster")

FeaturePlot(obj.tumor, features = c('ADGRB1', 'ADGRB2', "ADGRB3"), order = T)
color_scheme2 <- c('OC-like' = '#e41a1c', 'NPC-like' = '#377eb8', 'MES-like' = '#4daf4a',
                   'G2M' = '#ffff33',
                   'G1S' = '#a65628', 'AC-like' = '#f781bf')
order <- c(
           'G1S',
           'G2M','MES-like','OC-like','NPC-like', 'AC-like')
order <- rev(order)
DimPlot(object = obj.tumor, 
        group.by = "cluster", 
        cols = color_scheme2,
        label = T, 
        pt.size = .5, order = order,
        label.size = 0) & NoAxes()
ggsave("./figures/Fig3D_sc_tmr.DimPlot.pdf", width = 7, height = 6)

genes_to_plot = c("RND3", "CCND2", "STMN2", 
                  "APOE", "CRYAB", 'CLU',
                  "VIM", "S100A10", "IGFBP2", "CD74",
                  "OLIG1", "POLR2F", "DLL3", 
                  "TYMS", "PCNA", "TOP2A", "UBE2C")
genes_to_plot = rev(genes_to_plot)

tmp = DotPlot(obj.tumor, 
              features = genes_to_plot,
              cols = "RdBu",
              group.by = "cluster")
tmp + 
  scale_y_discrete(limits = c('NPC-like','AC-like', 'MES-like', 'OC-like', 'G1S',
                              'G2M')) + 
  coord_flip() + 
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(size=20, angle=45, hjust=1),
        axis.text.y = element_text(size=20))

ggsave("figures/DotPlot_marker_TumorCellStates.6CS.pdf", width = 6, height = 6.5)

# Fig 3E tumor states comparison ---- 
library(ggbeeswarm)
count.all <- as.data.frame(table(obj.tumor$cluster, obj.tumor$orig.ident))
count.wide <- reshape2::dcast(count.all, Var2 ~ Var1)
rownames(count.wide) <- count.wide$Var2
frac.all <- apply(count.wide[, -1], 1, function(x)x/sum(x))
frac.long <- reshape2::melt(frac.all)
colnames(frac.long) <- c("Var2", "Cohort_ID", "frac")

xxx <- read_delim("data/scRNA_cluster.txt")
colnames(frac.all)[1] <- "Cohort_ID"
frac.tmr <- merge(frac.long, xxx, by = "Cohort_ID")

NAN_plot <- ggplot(data=frac.tmr,aes(x=NMF_Cluster, y=frac)) + theme_classic() 
#NAN_plot
NAN_plot <- NAN_plot + 
  geom_boxplot(data=frac.tmr,aes(x = factor(NMF_Cluster, levels = c("MIX", "PPR", "IME", "NEU")), y = frac),width = 0.6,size=0.5,fill="transparent", outlier.color = "white")+
  geom_quasirandom(data=frac.tmr,aes(x = factor(NMF_Cluster, levels = c("MIX", "PPR", "IME", "NEU")), y = frac, color = NMF_Cluster, fill = NMF_Cluster),width = 0.25,size=0.75,alpha=.6,stroke=0.8, varwidth = T) + 
  facet_grid(.~ factor(Var2, levels = c("MES-like", "AC-like", "NPC-like", "OC-like", "OC_prey-like", "Phagocyte-like", 
                                        "G1S", "G2M"))) 
NAN_plot
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0,0.91),breaks = seq(-10,20,0.2)) 
NAN_plot<- NAN_plot +ylab("Fraction") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(c("IME" = "#ff7f00", "NEU" = "#1f78b4", "PPR" =  "#e31a1c", "MIX" = "#33a02c")))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(c("IME" = "#ff7f00", "NEU" = "#1f78b4", "PPR" =  "#e31a1c", "MIX" = "#33a02c")))

figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
NAN_plot
ggsave(file = "./figures/Fig3E_boxplot_compare_tumor.frac.pdf", plot = figure_2, bg = "white", width = 23, height = 9, units = "cm", dpi = 600)


library(ggpubr)
compare_means(method = "anova", frac ~ NMF_Cluster,   
              data = frac.tmr %>% filter(Var2 == "MES-like") ,paired = F )
compare_means(method = "anova", frac ~ NMF_Cluster,   
              data = frac.tmr %>% filter(Var2 == "AC-like"),paired = F )
compare_means(method = "anova", frac ~ NMF_Cluster,  
              data = frac.tmr %>% filter(Var2 == "NPC-like"),paired = F )
compare_means(method = "anova", frac ~ NMF_Cluster,  
              data = frac.tmr %>% filter(Var2 == "OC-like"),paired = F )

compare_means(method = "anova", frac ~ NMF_Cluster,   
              data = frac.tmr %>% filter(Var2 == "OC_prey-like") ,paired = F )
compare_means(method = "anova", frac ~ NMF_Cluster,   
              data = frac.tmr %>% filter(Var2 == "Phagocyte-like"),paired = F )
compare_means(method = "anova", frac ~ NMF_Cluster,  
              data = frac.tmr %>% filter(Var2 == "G1/S"),paired = F )
compare_means(method = "anova", frac ~ NMF_Cluster,  
              data = frac.tmr %>% filter(Var2 == "G2/M"),paired = F )


compare_means(method = "kruskal.test", frac ~ NMF_Cluster,   
              data = frac.tmr %>% filter(Var2 == "MES-like") ,paired = F )
compare_means(method = "kruskal.test", frac ~ NMF_Cluster,   
              data = frac.tmr %>% filter(Var2 == "AC-like"),paired = F )
compare_means(method = "kruskal.test", frac ~ NMF_Cluster,  
              data = frac.tmr %>% filter(Var2 == "NPC-like"),paired = F )
compare_means(method = "kruskal.test", frac ~ NMF_Cluster,  
              data = frac.tmr %>% filter(Var2 == "OC-like"),paired = F )

compare_means(method = "kruskal.test", frac ~ NMF_Cluster,   
              data = frac.tmr %>% filter(Var2 == "OC_prey-like") ,paired = F )
compare_means(method = "kruskal.test", frac ~ NMF_Cluster,   
              data = frac.tmr %>% filter(Var2 == "Phagocyte-like"),paired = F )
compare_means(method = "kruskal.test", frac ~ NMF_Cluster,  
              data = frac.tmr %>% filter(Var2 == "G1S"),paired = F )
compare_means(method = "kruskal.test", frac ~ NMF_Cluster,  
              data = frac.tmr %>% filter(Var2 == "G2M"),paired = F )

# Fig 3F Expression comparison between all cell types ---- 
library(Seurat)
library(scCustomize)
library(tidyverse)
obj.all <- readRDS("data/scRNA_revise/0323_Merged.Major.18MutA.QC.rds")
cgpa.lymphocyte <- readRDS("data/scRNA_revise/0422_merged_lymphocyte_annotated.rds")

cgpa.tumor <- readRDS("data/scRNA_revise/0316_Merged.18MutA.tumor.rds")
xxx <- read_delim("data/scRNA_cluster.txt")
tumor.meta <- cgpa.tumor@meta.data
tumor.meta$cid <- rownames(tumor.meta)
aaa <- tumor.meta[, c("cid", "cluster")]

lym.meta <- cgpa.lymphocyte@meta.data
lym.meta <- lym.meta[!(lym.meta$cluster %in% c("tumor_like", "Unknown")), ]

lym.meta$cid <- rownames(lym.meta)
aaa <- rbind(aaa, lym.meta[, c("cid", "cluster")])
table(aaa$cluster)

all.meta <- obj.all@meta.data
alll <- all.meta
alll$cid <- rownames(alll)

# alll$cluster <- alll$cell_type
blll <- merge(alll, aaa, by = "cid", all.x = TRUE)
blll$cluster <- ifelse(is.na(blll$cluster), blll$cell_type, blll$cluster)
table(blll$cluster)
rownames(blll) <- blll$cid
blll <- blll[rownames(all.meta), ]
sum(rownames(blll) == rownames(all.meta))

obj.all@meta.data <- blll
obj.all$bulkProt <- ifelse(obj.all$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "MIX"], "AFM",
  ifelse(obj.all$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "PPR"], "PPR",
    ifelse(obj.all$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "IME"], "IME", "NEU")
  )
)

obj.all$use <- ifelse(obj.all$cluster %in% c("Tumor", "Lyphm", "myeloid_like"), "No", "Yes")
obj.use <- subset(obj.all, subset = use == "Yes")

obj.use$plotcluster <- obj.use$cluster
obj.use$plotcluster <- ifelse(obj.use$cell_type == "Lymph", "Lymph", obj.use$plotcluster)
obj.use$ttt <- paste0(obj.use$plotcluster, "_", obj.use$bulkProt)

plt_genes <- c(
  "S100A10", "CHI3L2", # MES-like
  "APOE", "CLU", # AC-like
  "RND3", "CCND2", # NPC-like
  "OLIG1", "POLR2F", # OC-like
  "TYMS", "PCNA", "TOP2A", "UBE2C", # Cycling
  "CSF1R", "P2RY12", "CD68", "S100A9", # Myeloid
  "CD3E", "CD8A", "CD4", "FOXP3", "NKG7", "CD79A", "IGKC", # Lymphocyte
  "MAG", # Oligo
  "GBP2", "GBP1", "SERPINA3", "ISG15", "IFI44L", "IFI6", "SERPING1", "TYMP", "CD44", "CD74", "C1R",
  "HLA-DRA", "HLA-DRB1", "HLA-A", "B2M", "CD163", "TGFBI", "GZMK", "IGHG1", "PDCD1", "CD274"
)

df_gns <- data.frame(plt_genes, cls = c(
  rep("1Celltype", 24),
  rep("2IMEbulk", 11),
  rep("3ImmuneHot", 10)
))

xx <- rev(c("NEU", "IME", "PPR", "AFM"))
Idents(obj.use) <- obj.use$plotcluster
DotPlot(obj.use,
  group.by = "ttt", features = split(plt_genes, df_gns$cls),
  cols = "RdBu"
) + # , col.min = -4, col.max = 4) +
  scale_y_discrete(limits = rev(c(
    paste0("MES-like_", xx),
    paste0("AC-like_", xx),
    paste0("NPC-like_", xx),
    paste0("OC-like_", xx),
    paste0("G1S_", xx),
    paste0("G2M_", xx),
    paste0("Myeloid_", xx),
    paste0("Lymph_", xx),
    paste0("Oligo_", xx)
  ))) +
  theme(
    panel.border = element_rect(color = "black", linewidth = 1),
    panel.spacing = unit(1, "mm"),
    strip.text = element_text(margin = margin(b = 3, unit = "mm")),
    strip.placement = "outlet",
    legend.position = "bottom",
    axis.line = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain", color = "black"),
  ) + labs(x = "", y = "")
ggsave("figures/Fig3F_sc.compare_dotplot_majorLegend.pdf", width = 28, height = 25, units = "cm")