#### Visualization codes for Main Figure 3; 2024-12-09
#### Single-cell analysis revealed distinct lymphocyte infiltration and enhanced immune activity in the IDHm-IME gliomas
#### Author: Jihong TANG; Jiguang WANG

# Fig 3A - Major clustering ---- 
library(Seurat)
library(tidyverse)
obj.use <- readRDS("../../../scRNA_MutA/Windows_Data/Sep18_Version_Results/Merged.Major.18MutA.seurat.object.QC.rds")

# visulize the clusters 
color_scheme2 <- c('Tumor' = '#ea3323', 'Oligo' = '#8ab0d0', 'Myeloid' = '#b0d667',
                   'Lymph' = '#ed926b')
order <- c("Tumor", "Oligo", "Myeloid", "Lymph")

order <- rev(order)
DimPlot(object = obj.use, 
        group.by = "cell_type", 
        cols = color_scheme2,
        label = T, 
        pt.size = .5, #order = order,
        label.size = 0) & NoAxes()
ggsave("./figures/Fig3A_sc_major.DimPlot.pdf", width = 7, height = 6)

# Fig 3B - downsampling four samples ---- 
library(Seurat)
library(tidyverse)
obj.use <- readRDS("../../../scRNA_MutA/Windows_Data/Sep18_Version_Results/Merged.Major.18MutA.seurat.object.QC.rds")
DimPlot(obj.use)

xxx <- read_delim("data/Fig3/scRNA_cluster.txt")
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
        pt.size = .5, cols = c("IME" = "#ff7f00", "NEU" = "#1f78b4", "PPR" =  "#e31a1c", "MIX" = "#33a02c")
) + NoAxes()

ggsave("figures/Fig3B_downsampling.pdf", width = 8.5, height = 6.5)

# Fig 3C major cell fraction comparison ----
# now in windows workstation
load("../Fig_codes/data/Fig3/0318_Figure3F_H_rawdata.Rdata")
xxx <- read_delim("data/Fig3/scRNA_cluster.txt")
colnames(frac.all)[1] <- "Cohort_ID"
frac.all <- merge(frac.all, xxx, by = "Cohort_ID")
library(ggbeeswarm)
NAN_plot <- ggplot(data = frac.all, aes(x = NMF_Cluster, y = frac)) +
  theme_classic()
# NAN_plot
NAN_plot <- NAN_plot +
  geom_boxplot(data = frac.all, aes(x = factor(NMF_Cluster, levels = c("MIX", "PPR", "IME", "NEU")), y = frac), width = 0.6, size = 0.5, fill = "transparent", outlier.color = "white") +
  geom_quasirandom(data = frac.all, aes(x = factor(NMF_Cluster, levels = c("MIX", "PPR", "IME", "NEU")), y = frac, color = NMF_Cluster, fill = NMF_Cluster), width = 0.25, size = 0.75, alpha = .6, stroke = 0.8, varwidth = T) +
  facet_grid(. ~ factor(Var2))
NAN_plot
NAN_plot <- NAN_plot + theme(
  panel.background = element_rect(fill = "transparent", color = "transparent"), plot.margin = unit(c(2, 2, 2, 2), "lines"),
  plot.title = element_text(size = 34, vjust = 0.5, hjust = 0.5, face = "bold.italic", color = "transparent"), text = element_text(size = 14, face = "bold"),
  legend.key.width = unit(0.6, "cm"), legend.key.height = unit(0.6, "cm"), legend.position = "none", legend.text = element_text(size = 14, hjust = 0, face = "bold"),
  axis.text.x = element_text(size = 12, face = "bold", color = "black"), axis.text.y = element_text(size = 12, face = "bold", color = "black"),
  axis.title.x = element_text(size = 14, face = "plain", color = "black"), axis.title.y = element_text(size = 14, face = "plain", color = "black")
)

NAN_plot <- NAN_plot + scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks = seq(-10, 20, 0.2))
# NAN_plot <- NAN_plot + scale_x_discrete(limits = c("MES", "NotMES"))
NAN_plot <- NAN_plot + ylab("Fraction") + xlab(NULL)
# NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(c("IME" = "#fdbf6f", "NEU" = "#a6cee3", "PPR" =  "#fb9a99", "MIX" = "#b2df8a")))
# NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(c("IME" = "#fdbf6f", "NEU" = "#a6cee3", "PPR" =  "#fb9a99", "MIX" = "#b2df8a")))
NAN_plot <- NAN_plot + scale_fill_manual(name = NULL, values = c(c("IME" = "#ff7f00", "NEU" = "#1f78b4", "PPR" = "#e31a1c", "MIX" = "#33a02c")))
NAN_plot <- NAN_plot + scale_color_manual(name = NULL, values = c(c("IME" = "#ff7f00", "NEU" = "#1f78b4", "PPR" = "#e31a1c", "MIX" = "#33a02c")))

# NAN_plot <- NAN_plot +coord_flip()
figure_2 <- rbind(ggplotGrob(NAN_plot), size = "last")
NAN_plot
ggsave(file = "./figures/Fig3C_boxplot_compare_all.frac_4subtype.pdf", plot = figure_2, bg = "white", width = 13.5, height = 9, units = "cm", dpi = 600)

library(tidyverse)
library(ggpubr)
compare_means(method = "kruskal.test", frac ~ NMF_Cluster,   
              data = frac.all %>% filter(Var2 == "Lymph") ,paired = F )
compare_means(method = "kruskal.test", frac ~ NMF_Cluster,   
              data = frac.all %>% filter(Var2 == "Myeloid"),paired = F )
compare_means(method = "kruskal.test", frac ~ NMF_Cluster,  
              data = frac.all %>% filter(Var2 == "Oligo"),paired = F )
compare_means(method = "kruskal.test", frac ~ NMF_Cluster,  
              data = frac.all %>% filter(Var2 == "Tumor"),paired = F )

compare_means(method = "anova", frac ~ NMF_Cluster,   
              data = frac.all %>% filter(Var2 == "Lymph") ,paired = F )
compare_means(method = "anova", frac ~ NMF_Cluster,   
              data = frac.all %>% filter(Var2 == "Myeloid"),paired = F )
compare_means(method = "anova", frac ~ NMF_Cluster,  
              data = frac.all %>% filter(Var2 == "Oligo"),paired = F )
compare_means(method = "anova", frac ~ NMF_Cluster,  
              data = frac.all %>% filter(Var2 == "Tumor"),paired = F )

# Fig 3D tumor cell states ---- 
load("../../../scRNA_MutA/1216_MES_sc_comapre/Fig3_V4/0321_sc.combine.CGPA.MESNonMES.compare.RData")

tt <- cgpa.tumor@reductions[["umap"]]@cell.embeddings %>% as.data.frame()
remove <- rownames(tt)[tt$UMAP_2 > 6]
cgpa.tumor$subset <- ifelse(colnames(cgpa.tumor) %in% remove, "No", "Yes")
cgpa.tumor.new <- subset(cgpa.tumor, subset = subset == "Yes")
Idents(cgpa.tumor.new) <- "cluster"

DimPlot(cgpa.tumor.new)

ddd <- cgpa.tumor.new 
ddd <- RunUMAP(ddd, dims = 1:12, reduction = "harmony")
DimPlot(ddd, group.by = "cluster")

color_scheme2 <- c('OC-like' = '#e41a1c', 'NPC-like' = '#377eb8', 'MES-like' = '#4daf4a',
                   'G2/M' = '#ffff33',
                   'G1/S' = '#a65628', 'AC-like' = '#f781bf')
order <- c('OC-like','NPC-like', 'MES-like','AC-like', 
             'G1/S',
           'G2/M')
order <- rev(order)
DimPlot(object = ddd, 
        group.by = "cluster", 
        cols = color_scheme2,
        label = T, 
        pt.size = .5, order = order,
        label.size = 0) & NoAxes()
ggsave("./figures/Fig3D_sc_tmr.DimPlot_reduction12.pdf", width = 7, height = 6)

# Fig 3E tumor cell fraction compare ----
frac.tmr <- read_delim("../Fig_codes/data/Fig3/Fig3_cgga_scRNA_tmrstates_fractions.txt")
library(ggbeeswarm)
colnames(frac.tmr)[1] <- "Cohort_ID"
frac.tmr <- merge(frac.tmr, xxx, by = "Cohort_ID")

frac.tmr <- frac.tmr[!(frac.tmr$Var2 %in% c("Phagocyte-like", "OC_prey-like")), ]
NAN_plot <- ggplot(data=frac.tmr,aes(x=NMF_Cluster, y=frac)) + theme_classic() 
#NAN_plot
NAN_plot <- NAN_plot + 
  geom_boxplot(data=frac.tmr,aes(x = factor(NMF_Cluster, levels = c("MIX", "PPR", "IME", "NEU")), y = frac),width = 0.6,size=0.5,fill="transparent", outlier.color = "white")+
  geom_quasirandom(data=frac.tmr,aes(x = factor(NMF_Cluster, levels = c("MIX", "PPR", "IME", "NEU")), y = frac, color = NMF_Cluster, fill = NMF_Cluster),width = 0.25,size=0.75,alpha=.6,stroke=0.8, varwidth = T) + 
  facet_grid(.~ factor(Var2, levels = c("MES-like", "AC-like", "NPC-like", "OC-like", "OC_prey-like", "Phagocyte-like", 
                                        "G1/S", "G2/M"))) 
NAN_plot
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0,0.81),breaks = seq(-10,20,0.2)) 
#NAN_plot <- NAN_plot + scale_x_discrete(limits = c("MES", "NotMES"))
NAN_plot<- NAN_plot +ylab("Fraction") +xlab(NULL)
#NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(c("IME" = "#fdbf6f", "NEU" = "#a6cee3", "PPR" =  "#fb9a99", "MIX" = "#b2df8a")))
#NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(c("IME" = "#fdbf6f", "NEU" = "#a6cee3", "PPR" =  "#fb9a99", "MIX" = "#b2df8a")))
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(c("IME" = "#ff7f00", "NEU" = "#1f78b4", "PPR" =  "#e31a1c", "MIX" = "#33a02c")))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(c("IME" = "#ff7f00", "NEU" = "#1f78b4", "PPR" =  "#e31a1c", "MIX" = "#33a02c")))

#NAN_plot <- NAN_plot +coord_flip()
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
NAN_plot
ggsave(file="./figures/Fig3E_boxplot_compare_tumor.frac_ProMES.pdf", plot=figure_2,bg = 'white', width =23, height = 9, units = 'cm', dpi = 600)

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
              data = frac.tmr %>% filter(Var2 == "G1/S"),paired = F )
compare_means(method = "kruskal.test", frac ~ NMF_Cluster,  
              data = frac.tmr %>% filter(Var2 == "G2/M"),paired = F )

# Fig 3F expression comparison ---- 
library(Seurat)
library(scCustomize)
library(tidyverse)
load("../../../scRNA_MutA/1216_MES_sc_comapre/Fig3_V4/0321_sc.combine.CGPA.MESNonMES.compare.RData")

xxx <- read_delim("data/Fig3/scRNA_cluster.txt")
#obj.use$bulkMES <- ifelse(obj.use$orig.ident %in% c("P698", "P714", "P737"), "MES", "NotMES")
cgpa.myeloid$bulkProt <- ifelse(cgpa.myeloid$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "MIX"], "MIX", 
                           ifelse(cgpa.myeloid$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "PPR"], "PPR",
                                  ifelse(cgpa.myeloid$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "IME"], "IME", "NEU")))
cgpa.lymphocyte$bulkProt <- ifelse(cgpa.lymphocyte$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "MIX"], "MIX", 
                           ifelse(cgpa.lymphocyte$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "PPR"], "PPR",
                                  ifelse(cgpa.lymphocyte$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "IME"], "IME", "NEU")))
cgpa.oligo$bulkProt <- ifelse(cgpa.oligo$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "MIX"], "MIX", 
                           ifelse(cgpa.oligo$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "PPR"], "PPR",
                                  ifelse(cgpa.oligo$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "IME"], "IME", "NEU")))
cgpa.tumor$bulkProt <- ifelse(cgpa.tumor$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "MIX"], "MIX", 
                           ifelse(cgpa.tumor$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "PPR"], "PPR",
                                  ifelse(cgpa.tumor$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "IME"], "IME", "NEU")))

plt_genes <- c("CD74","HLA-DRA", "HLA-DPA1", "HLA-DRB1", "HLA-A","HLA-B","HLA-C","B2M", "C1QA", "C1QB", "C1QC", 
               "SEPP1", "IFI44L","IFI6", "ISG15","SERPING1", "SPP1", "XIST", "VIM", "CD44",# tot
               "GBP2", "GBP1","CHI3L2", "C1R", "C1S", "IGFBP2", #tmr
               "F13A1","CD163","STAB1","MS4A6A","TMEM176B","DAB2","HLA-DRB5","CCL3L1",# mye,
               "PDCD1", "CD274", "GZMK", "CD8A", "CD4", "FOXP3"
)

df_gns <- data.frame(plt_genes, cls = c(rep("atot", 20), rep("btmr", 6),
                                        rep("cmye", 8),rep("dlym", 6)))

DotPlot(cgpa.myeloid, group.by = "bulkProt", features = split(plt_genes, df_gns$cls), 
        cols = "RdBu") + 
  scale_y_discrete(limits = c("NEU", "IME", "PPR", "MIX")) + 
  theme(
    panel.border = element_rect(color="black", linewidth =1),
    panel.spacing = unit(1, "mm"), 
    
    strip.text = element_text(margin=margin(b=3, unit="mm")),
    strip.placement = 'outlet', 
    
    legend.position = "none", 
    
    axis.line = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain", color = "black"),
    
  )+labs(x="", y="")
ggsave("figures/Fig3F_sc.compare_dotplot_myetmpV2.pdf", width = 24, height = 5.5, units = "cm")

DotPlot(cgpa.lymphocyte, group.by = "bulkProt", features = split(plt_genes, df_gns$cls), 
        cols = "RdBu") + 
  scale_y_discrete(limits = c("NEU", "IME", "PPR", "MIX")) + 
  theme(
    panel.border = element_rect(color="black", linewidth =1),
    panel.spacing = unit(1, "mm"), 
    
    strip.text = element_text(margin=margin(b=3, unit="mm")),
    strip.placement = 'outlet', 
    
    legend.position = "null", 
    
    axis.line = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain", color = "black"),
    
  )+labs(x="", y="")
ggsave("figures/Fig3F_sc.compare_dotplot_lymV2.pdf", width = 24, height = 5.5, units = "cm")

DotPlot(cgpa.oligo, group.by = "bulkProt", features = split(plt_genes, df_gns$cls), 
        cols = "RdBu") + 
  scale_y_discrete(limits = c("NEU", "IME", "PPR", "MIX")) + 
  theme(
    panel.border = element_rect(color="black", linewidth =1),
    panel.spacing = unit(1, "mm"), 
    
    strip.text = element_text(margin=margin(b=3, unit="mm")),
    strip.placement = 'outlet', 
    
    legend.position = "null", 
    
    axis.line = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain", color = "black"),
    
  )+labs(x="", y="")
ggsave("figures/Fig3F_sc.compare_dotplot_olgV2.pdf", width = 24, height = 5.5, units = "cm")

DotPlot(cgpa.tumor, group.by = "bulkProt", features = split(plt_genes, df_gns$cls), 
        cols = "RdBu") + 
  scale_y_discrete(limits = c("NEU", "IME", "PPR", "MIX")) + 
  theme(
    panel.border = element_rect(color="black", linewidth =1),
    panel.spacing = unit(1, "mm"), 
    
    strip.text = element_text(margin=margin(b=3, unit="mm")),
    strip.placement = 'outlet', 
    
    legend.position = "null", 
    
    axis.line = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain", color = "black"),
    
  )+labs(x="", y="")
ggsave("figures/Fig3F_sc.compare_dotplot_tmrV2.pdf", width = 24, height = 5.5, units = "cm")

unique(cgpa.tumor$cluster)
tmr.MES <- subset(cgpa.tumor, subset = cluster == "MES-like")
DotPlot(tmr.MES, group.by = "bulkProt", features = split(plt_genes, df_gns$cls), 
        cols = "RdBu") + 
  scale_y_discrete(limits = c("NEU", "IME", "PPR", "MIX")) + 
  theme(
    panel.border = element_rect(color="black", linewidth =1),
    panel.spacing = unit(1, "mm"), 
    
    strip.text = element_text(margin=margin(b=3, unit="mm")),
    strip.placement = 'outlet', 
    
    legend.position = "null", 
    
    axis.line = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain", color = "black"),
    
  )+labs(x="", y="")
ggsave("figures/Fig3F_sc.compare_dotplot_tmr.MES.pdf", width = 24, height = 5.5, units = "cm")

tmr.AC <- subset(cgpa.tumor, subset = cluster == "AC-like")
DotPlot(tmr.AC, group.by = "bulkProt", features = split(plt_genes, df_gns$cls), 
        cols = "RdBu") + 
  scale_y_discrete(limits = c("NEU", "IME", "PPR", "MIX")) + 
  theme(
    panel.border = element_rect(color="black", linewidth =1),
    panel.spacing = unit(1, "mm"), 
    
    strip.text = element_text(margin=margin(b=3, unit="mm")),
    strip.placement = 'outlet', 
    
    legend.position = "null", 
    
    axis.line = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain", color = "black"),
    
  )+labs(x="", y="")
ggsave("figures/Fig3F_sc.compare_dotplot_tmr.AC.pdf", width = 24, height = 5.5, units = "cm")

tmr.Phagocyte <- subset(cgpa.tumor, subset = cluster == "Phagocyte-like")
DotPlot(tmr.Phagocyte, group.by = "bulkProt", features = split(plt_genes, df_gns$cls), 
        cols = "RdBu") + 
  scale_y_discrete(limits = c("NEU", "IME", "PPR", "MIX")) + 
  theme(
    panel.border = element_rect(color="black", linewidth =1),
    panel.spacing = unit(1, "mm"), 
    
    strip.text = element_text(margin=margin(b=3, unit="mm")),
    strip.placement = 'outlet', 
    
    legend.position = "null", 
    
    axis.line = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain", color = "black"),
    
  )+labs(x="", y="")
ggsave("figures/Fig3F_sc.compare_dotplot_tmr.Phagocyte.pdf", width = 24, height = 5.5, units = "cm")

tmr.NPC <- subset(cgpa.tumor, subset = cluster == "NPC-like")
DotPlot(tmr.NPC, group.by = "bulkProt", features = split(plt_genes, df_gns$cls), 
        cols = "RdBu") + 
  scale_y_discrete(limits = c("NEU", "IME", "PPR", "MIX")) + 
  theme(
    panel.border = element_rect(color="black", linewidth =1),
    panel.spacing = unit(1, "mm"), 
    
    strip.text = element_text(margin=margin(b=3, unit="mm")),
    strip.placement = 'outlet', 
    
    legend.position = "null", 
    
    axis.line = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain", color = "black"),
    
  )+labs(x="", y="")
ggsave("figures/Fig3F_sc.compare_dotplot_tmr.NPC.pdf", width = 24, height = 5.5, units = "cm")

tmr.OC <- subset(cgpa.tumor, subset = cluster == "OC-like")
DotPlot(tmr.OC, group.by = "bulkProt", features = split(plt_genes, df_gns$cls), 
        cols = "RdBu") + 
  scale_y_discrete(limits = c("NEU", "IME", "PPR", "MIX")) + 
  theme(
    panel.border = element_rect(color="black", linewidth =1),
    panel.spacing = unit(1, "mm"), 
    
    strip.text = element_text(margin=margin(b=3, unit="mm")),
    strip.placement = 'outlet', 
    
    legend.position = "null", 
    
    axis.line = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain", color = "black"),
    
  )+labs(x="", y="")
ggsave("figures/Fig3F_sc.compare_dotplot_tmr.OC.pdf", width = 24, height = 5.5, units = "cm")

tmr.OC_prey <- subset(cgpa.tumor, subset = cluster == "OC_prey-like")
DotPlot(tmr.OC_prey, group.by = "bulkProt", features = split(plt_genes, df_gns$cls), 
        cols = "RdBu") + 
  scale_y_discrete(limits = c("NEU", "IME", "PPR", "MIX")) + 
  theme(
    panel.border = element_rect(color="black", linewidth =1),
    panel.spacing = unit(1, "mm"), 
    
    strip.text = element_text(margin=margin(b=3, unit="mm")),
    strip.placement = 'outlet', 
    
    legend.position = "null", 
    
    axis.line = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain", color = "black"),
    
  )+labs(x="", y="")
ggsave("figures/Fig3F_sc.compare_dotplot_tmr.OC_prey.pdf", width = 24, height = 5.5, units = "cm")

tmr.G1S <- subset(cgpa.tumor, subset = cluster %in% c("G1/S"))
DotPlot(tmr.G1S, group.by = "bulkProt", features = split(plt_genes, df_gns$cls), 
        cols = "RdBu") + 
  scale_y_discrete(limits = c("NEU", "IME", "PPR", "MIX")) + 
  theme(
    panel.border = element_rect(color="black", linewidth =1),
    panel.spacing = unit(1, "mm"), 
    
    strip.text = element_text(margin=margin(b=3, unit="mm")),
    strip.placement = 'outlet', 
    
    legend.position = "null", 
    
    axis.line = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain", color = "black"),
    
  )+labs(x="", y="")
ggsave("figures/Fig3F_sc.compare_dotplot_tmr.G1S.pdf", width = 24, height = 5.5, units = "cm")

tmr.G2M <- subset(cgpa.tumor, subset = cluster %in% c("G2/M"))
DotPlot(tmr.G2M, group.by = "bulkProt", features = split(plt_genes, df_gns$cls), 
        cols = "RdBu") + 
  scale_y_discrete(limits = c("NEU", "IME", "PPR", "MIX")) + 
  theme(
    panel.border = element_rect(color="black", linewidth =1),
    panel.spacing = unit(1, "mm"), 
    
    strip.text = element_text(margin=margin(b=3, unit="mm")),
    strip.placement = 'outlet', 
    
    legend.position = "null", 
    
    axis.line = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, face = "plain", color = "black"),
    
  )+labs(x="", y="")
ggsave("figures/Fig3F_sc.compare_dotplot_tmr.G2M.pdf", width = 24, height = 5.5, units = "cm")