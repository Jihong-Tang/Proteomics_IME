# 3A - Major clustering ---- 
library(Seurat)
library(tidyverse)
obj.use <- readRDS("../../../scRNA_MutA/Windows_Data/Sep18_Version_Results/Merged.Major.18MutA.seurat.object.QC.rds")

### 1.1 - Dimplot UMAP ---- 
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
ggsave("./figures/Fig3/1023_sc_major.DimPlot.pdf", width = 7, height = 6)

# 3B - downsampling four samples ---- 
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

ggsave("figures/Fig3/Fig3B_downsampling.pdf", width = 8.5, height = 6.5)

# Fig 3D 1127 version cutoff OC and OC_prey_like ---- 
load("../../../scRNA_MutA/1216_MES_sc_comapre/Fig3_V4/0321_sc.combine.CGPA.MESNonMES.compare.RData")

tt <- cgpa.tumor@reductions[["umap"]]@cell.embeddings %>% as.data.frame()
remove <- rownames(tt)[tt$UMAP_2 > 6]
cgpa.tumor$subset <- ifelse(colnames(cgpa.tumor) %in% remove, "No", "Yes")
cgpa.tumor.new <- subset(cgpa.tumor, subset = subset == "Yes")
Idents(cgpa.tumor.new) <- "cluster"

DimPlot(cgpa.tumor.new)

ddd <- cgpa.tumor.new 
# ddd <- NormalizeData(object = ddd, normalization.method = "LogNormalize", scale.factor = 10000)
# ddd <- FindVariableFeatures(object = ddd, selection.method = "vst", nfeatures = 2000)
# ddd <- ScaleData(object = ddd)
# ddd <- RunPCA(object = ddd, features = VariableFeatures(object = ddd), npcs = 50)
# DimPlot(ddd, reduction = "pca")
# ElbowPlot(ddd)
# 
# library(harmony)
# ddd <- RunHarmony(object = ddd, group.by.vars = "orig.ident", max.iter.harmony = 25)
# ddd <- RunTSNE(object = ddd, dims = 1:15, reduction = "harmony")
# ddd <- RunUMAP(object = ddd, dims = 1:15, reduction = "harmony")
ddd <- RunUMAP(ddd, dims = 1:12, reduction = "harmony")
DimPlot(ddd, group.by = "cluster")
FeaturePlot(ddd, features = c("CD74", "CD44"), order = T)

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
ggsave("./figures/Fig3/1127_sc_tmr.DimPlot_reduction12.pdf", width = 7, height = 6)

genes_to_plot = c("RND3", "CCND2", "STMN2", 
                  "APOE", "CRYAB", 'CLU',
                  "VIM", "S100A10", "IGFBP2", "CD74",
                  "OLIG1", "POLR2F", "DLL3", 
                  "TYMS", "PCNA", "TOP2A", "UBE2C")
genes_to_plot = rev(genes_to_plot)

tmp = DotPlot(ddd, 
              features = genes_to_plot,
              cols = "RdBu",
              group.by = "cluster")
tmp + 
  scale_y_discrete(limits = c('NPC-like','AC-like', 'MES-like', 'OC-like', 'G1/S',
                              'G2/M')) + 
  coord_flip() + 
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(size=20, angle=45, hjust=1),
        axis.text.y = element_text(size=20))

ggsave("figures/Fig3/1128_DotPlot_marker_TumorCellStates.6CS.pdf", width = 6, height = 6.5)

# Fig 3E single cell fraction comparison; 0316 ----
# now in windows workstation 
load("../Fig_codes/data/Fig3/0318_Figure3F_H_rawdata.Rdata")
xxx <- read_delim("data/Fig3/scRNA_cluster.txt")
colnames(frac.all)[1] <- "Cohort_ID"
frac.all <- merge(frac.all, xxx, by = "Cohort_ID")
library(ggbeeswarm)
NAN_plot <- ggplot(data=frac.all,aes(x=NMF_Cluster, y=frac)) + theme_classic() 
#NAN_plot
NAN_plot <- NAN_plot + 
  geom_boxplot(data=frac.all,aes(x = factor(NMF_Cluster, levels = c("MIX", "PPR", "IME", "NEU")), y = frac),width = 0.6,size=0.5,fill="transparent", outlier.color = "white")+
  geom_quasirandom(data=frac.all,aes(x = factor(NMF_Cluster, levels = c("MIX", "PPR", "IME", "NEU")), y = frac, color = NMF_Cluster, fill = NMF_Cluster),width = 0.25,size=0.75,alpha=.6,stroke=0.8, varwidth = T) + 
  facet_grid(.~factor(Var2)) 
NAN_plot
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks = seq(-10,20,0.2)) 
#NAN_plot <- NAN_plot + scale_x_discrete(limits = c("MES", "NotMES"))
NAN_plot<- NAN_plot +ylab("Fraction") +xlab(NULL)
#NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(c("IME" = "#fdbf6f", "NEU" = "#a6cee3", "PPR" =  "#fb9a99", "MIX" = "#b2df8a")))
#NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(c("IME" = "#fdbf6f", "NEU" = "#a6cee3", "PPR" =  "#fb9a99", "MIX" = "#b2df8a")))
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(c("IME" = "#ff7f00", "NEU" = "#1f78b4", "PPR" =  "#e31a1c", "MIX" = "#33a02c")))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(c("IME" = "#ff7f00", "NEU" = "#1f78b4", "PPR" =  "#e31a1c", "MIX" = "#33a02c")))

#NAN_plot <- NAN_plot +coord_flip()
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
NAN_plot
ggsave(file="./figures/Fig3/1223_boxplot_compare_all.frac_4subtype.pdf", plot=figure_2,bg = 'white', width =13.5, height = 9, units = 'cm', dpi = 600)

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

frac.tmr <- read_delim("../Fig_codes/data/Fig3/Fig3_cgga_scRNA_tmrstates_fractions.txt")
library(ggbeeswarm)
colnames(frac.tmr)[1] <- "Cohort_ID"
frac.tmr <- merge(frac.tmr, xxx, by = "Cohort_ID")

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
ggsave(file="./figures/Fig3/1023_boxplot_compare_tumor.frac_ProMES.pdf", plot=figure_2,bg = 'white', width =24, height = 9, units = 'cm', dpi = 600)

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
ggsave(file="./figures/Fig3/1127_boxplot_compare_tumor.frac_ProMES.pdf", plot=figure_2,bg = 'white', width =23, height = 9, units = 'cm', dpi = 600)




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
              data = frac.tmr %>% filter(Var2 == "G1/S"),paired = F )
compare_means(method = "kruskal.test", frac ~ NMF_Cluster,  
              data = frac.tmr %>% filter(Var2 == "G2/M"),paired = F )

# Fig 3F expression comparison ---- 
## 1023 updates ---- 
# compare only between IDHmut MES and IDHmut non-MES 
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
               "GBP2", "GBP1","CHI3L2", "C1R", "C1S", "IGFBP2", "CXCL9", "CXCL10", #tmr
               "F13A1","CD163","STAB1","MS4A6A","TMEM176B","DAB2","HLA-DRB5","CCL3L1",# mye,
               "PDCD1", "CD274", "GZMK", "CD8A", "CD4", "FOXP3", "CXCR3", "CXCR5", "CTLA4", "TIGIT", "LAG3"
)

df_gns <- data.frame(plt_genes, cls = c(rep("atot", 20), rep("btmr", 8),
                                        rep("cmye", 8),rep("dlym", 11)))

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
ggsave("figures/Fig3/1023_sc.compare_dotplot_myetmpV2.pdf", width = 24, height = 5.5, units = "cm")

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
ggsave("figures/Fig3/1023_sc.compare_dotplot_lymV2.pdf", width = 24, height = 5.5, units = "cm")

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
ggsave("figures/Fig3/1023_sc.compare_dotplot_olgV2.pdf", width = 24, height = 5.5, units = "cm")

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
ggsave("figures/Fig3/1023_sc.compare_dotplot_tmrV2.pdf", width = 24, height = 5.5, units = "cm")

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
ggsave("figures/Fig3/1023_sc.compare_dotplot_tmr.MES.pdf", width = 24, height = 5.5, units = "cm")

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
ggsave("figures/Fig3/1023_sc.compare_dotplot_tmr.AC.pdf", width = 24, height = 5.5, units = "cm")

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
ggsave("figures/Fig3/1023_sc.compare_dotplot_tmr.Phagocyte.pdf", width = 24, height = 5.5, units = "cm")

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
ggsave("figures/Fig3/1023_sc.compare_dotplot_tmr.NPC.pdf", width = 24, height = 5.5, units = "cm")

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
ggsave("figures/Fig3/1023_sc.compare_dotplot_tmr.OC.pdf", width = 24, height = 5.5, units = "cm")

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
ggsave("figures/Fig3/1023_sc.compare_dotplot_tmr.OC_prey.pdf", width = 24, height = 5.5, units = "cm")

tmr.Cycle <- subset(cgpa.tumor, subset = cluster %in% c("G1/S", "G2/M"))
DotPlot(tmr.Cycle, group.by = "bulkProt", features = split(plt_genes, df_gns$cls), 
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
ggsave("figures/Fig3/1023_sc.compare_dotplot_tmr.Cycle.pdf", width = 24, height = 5.5, units = "cm")

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
ggsave("figures/Fig3/1127_sc.compare_dotplot_tmr.G1S.pdf", width = 24, height = 5.5, units = "cm")

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
ggsave("figures/Fig3/1127_sc.compare_dotplot_tmr.G2M.pdf", width = 24, height = 5.5, units = "cm")


# Supple Tmr TME correlation ---- 
library(tidyverse)

# 1- Load data ---- 
tumor.meta <- readRDS("../../../scRNA_MutA//Windows_Data/Sep18_Version_Results/18MutA_tumor.meta.Oct13V2.rds")
tb.tmr <- table(tumor.meta$orig.ident, tumor.meta$cluster) %>% as.data.frame()
df.tmr <- reshape(tb.tmr, idvar = "Var1", timevar = "Var2", direction = "wide") %>% as.data.frame()
row.names(df.tmr) <- df.tmr$Var1
df.tmr <- df.tmr[, -1]
colnames(df.tmr) <- gsub("Freq.", "", colnames(df.tmr))

df.tmr.count <- apply(df.tmr, 1, sum) %>% as.data.frame()
colnames(df.tmr.count) <- "Tumor"

myeloid.meta <- readRDS("../../../scRNA_MutA/Windows_Data/Sep18_Version_Results/18MutA_myeloid.meta.Sep18V1.rds")

tb.mye <- table(myeloid.meta$orig.ident, myeloid.meta$cluster) %>% as.data.frame()
df.mye <- reshape(tb.mye, idvar = "Var1", timevar = "Var2", direction = "wide") %>% as.data.frame()
row.names(df.mye) <- df.mye$Var1
df.mye <- df.mye[, -1]
df.mye <- df.mye[, -5]
colnames(df.mye) <- c("MDM", "MG", "Mono", "Neu")

lymph.meta <- readRDS("../../../scRNA_MutA/Windows_Data/Sep18_Version_Results/18MutA_lymph.metaSep18.rds")

tb.lym <- table(lymph.meta$orig.ident, lymph.meta$cluster) %>% as.data.frame()
df.lym <- reshape(tb.lym, idvar = "Var1", timevar = "Var2", direction = "wide") %>% as.data.frame()
row.names(df.lym) <- df.lym$Var1
df.lym <- df.lym[, -1]
df.lym2 <- df.lym %>% mutate(Tcell = Freq.Treg + Freq.CD4T + Freq.CD8T ) %>% 
  mutate(NK = Freq.NK, NKT = Freq.NKT, Bcell = Freq.Bcell) %>% 
  select(Tcell, NK, NKT, Bcell)


oligo.meta <- readRDS("../../../scRNA_MutA/Windows_Data/Sep18_Version_Results/Merged.meta.18MutA.rds")
oligo.meta <- oligo.meta[oligo.meta$cell_type == "Oligo", ]
tb.oligo <- table(oligo.meta$orig.ident) %>% as.data.frame()
df.oligo <- tb.oligo %>% as.data.frame()
row.names(df.oligo) <- df.oligo$Var1
df.oligo <- df.oligo[, -1] %>% as.data.frame()
row.names(df.oligo) <- tb.oligo$Var1; colnames(df.oligo) <- "Oligo"
tmp.df <- data.frame(matrix(rep(0,1), ncol=1))
colnames(tmp.df) <- colnames(df.oligo)
row.names(tmp.df) <- row.names(df.tmr)[!(row.names(df.tmr) %in% row.names(df.oligo))]
df.oligo <- rbind(df.oligo, tmp.df)

# 2 - Combine and calculate the fractions ----
df.mye <- df.mye[row.names(df.tmr.count),]; df.lym2 <- df.lym2[row.names(df.tmr.count),]; 
df.oligo <- df.oligo[row.names(df.tmr.count),] %>% as.data.frame(); 
colnames(df.oligo) <- "Oligo"; rownames(df.oligo) <- rownames(df.tmr.count)

df.cell.all <- cbind(df.tmr.count, df.mye)
df.cell.all <- cbind(df.mye)
df.cell.all <- cbind(df.cell.all, df.lym2)
df.cell.all <- cbind(df.cell.all, df.oligo)
frac.cell.all <- apply(df.cell.all, 1, function(x)x/sum(x)) %>% t() %>% as.data.frame()

df.tmr2 <- df.tmr[, c("AC-like", "MES-like", "OC-like", "NPC-like", "G1/S", "G2/M")]
frac.tmr.states <- apply(df.tmr2, 1, function(x)x/sum(x)) %>% t() %>% as.data.frame()
#frac.tmr.states$MES <- frac.tmr.states$Tmr.State3 + frac.tmr.states$Tmr.State4
row.names(frac.tmr.states) == row.names(frac.cell.all)

# 3- calculate the correlations ---- 
tmrstates <- colnames(frac.tmr.states)
cells <- colnames(frac.cell.all)

get_cor_tmr <- function(tmr, cell){
  scc.test <- cor.test(frac.tmr.states[[tmr]], frac.cell.all[[cell]], method = "spearman")
  scc <- scc.test$estimate; pvalue.scc <- scc.test$p.value
  pcc.test <- cor.test(frac.tmr.states[[tmr]], frac.cell.all[[cell]], method = "pearson")
  pcc <- pcc.test$estimate; pvalue.pcc <- pcc.test$p.value
  
  out <- data.frame(tmr, cell, scc, pcc, pvalue.scc, pvalue.pcc)
  return(out)
}

get_cor_all <- function(cell){
  res_tmr <- do.call(rbind, lapply(tmrstates, function(x)get_cor_tmr(x, cell)))
}

res_cor_all <- do.call(rbind, lapply(cells, function(x)get_cor_all(x)))

# 4- Visualization ---- 

## 4.1 stacked barplot ---- 
pt_stack <- data.frame()

frac.tmr.filter <- frac.tmr.states
for (i in 1:ncol(frac.tmr.filter)){
  tmp <- data.frame(row.names(frac.tmr.filter), frac.tmr.filter[, i])
  colnames(tmp) <- c("Sample.ID", "Frac")
  #tmp <- merge(tmp, , by = "Sample.ID")
  tmp$Group <- colnames(frac.tmr.filter)[i]
  pt_stack <- rbind(pt_stack, tmp)
}

col_order <- frac.tmr.filter[order(frac.tmr.filter$`MES-like` + frac.tmr.filter$`AC-like`, decreasing = F), ] %>% row.names()
color_scheme2 <- c('OC-like' = '#e41a1c', 'NPC-like' = '#377eb8', 'MES-like' = '#4daf4a',
                   'Phagocyte-like' = '#984ea3', 'OC_prey-like' = '#ff7f00', 'G2/M' = '#ffff33',
                   'G1/S' = '#a65628', 'AC-like' = '#f781bf')
pt_stack %>% 
  #filter(Region == "A") %>% 
  mutate(Group = factor(Group, levels = c('AC-like', 'MES-like', 'NPC-like','OC-like',
                                          'OC_prey-like','Phagocyte-like',  'G1/S',
                                          'G2/M'))) %>% 
  ggplot(aes(x = Sample.ID, y = Frac, fill = Group)) + 
  geom_bar(position="fill", stat="identity", show.legend = T, width = .6, color = "black") + 
  #geom_text(aes(x = Sample.ID, y = 1, label = Pred), hjust =1.5, angle = 90, size = 2) + 
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.2), expand = c(0, 0)) + 
  scale_x_discrete(limits = col_order) + 
  scale_fill_manual(values = color_scheme2) + 
  labs(x = "", y = "% Tumor States") + 
  theme_classic() + 
  #coord_flip() +
  #facet_grid(.~Region) + 
  # theme(
  #   axis.text.x = element_text(size = 12, angle = 45, color = "black", hjust = 1),
  #   #axis.text.x = element_blank(),
  #   axis.ticks.x = element_blank(),
  #   axis.text.y = element_text(size = 12, color = "black"),
  #   axis.title = element_text(size = 14)
  # ) 
  theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(1,4,1,1),'lines'),
        #text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),
        legend.key.width=unit(.5,'cm'),legend.key.height=unit(0.5,'cm'),legend.position="bottom",
        legend.text=element_text(size=12,face='plain'),
        axis.text.y=element_text(size=16,vjust=0.5,hjust=1,face='plain',color='black'),
        legend.title=element_text(size=12,vjust=0.5,hjust=0,face='plain'),
        axis.text.x=element_text(size=16,angle=45,vjust=,hjust=1,face='plain',color='black'),
        axis.title.x=element_text(size=20,vjust=0,hjust=0.5,face='plain',color='black'),
        axis.title.y=element_text(size=20,hjust=0.5,vjust=2,face='plain',color='black'))

ggsave("figures/stack_MutA_tmr_CellState_1016.pdf", height = 6, width = 10)


## 4.2 correlation heatmap ---- 
res_cor_plt <- res_cor_all %>% 
  filter(tmr != "MES") 

tmp.pcc <- res_cor_plt[, c("tmr", "cell", "pcc")]
pcc.mt <- reshape(tmp.pcc, idvar = "tmr", timevar = "cell", direction = "wide") %>% as.data.frame()
row.names(pcc.mt) <- pcc.mt$tmr; pcc.mt <- pcc.mt[,-1]
colnames(pcc.mt) <- c( "Marcophage", "Microglia", "Monocyte", "Neutrophil", "T cells", "NK cells", "NKT cells", "B cells", "Oligodendrocyte" )

tmp.pvalue.pcc <- res_cor_plt[, c("tmr", "cell", "pvalue.pcc")]
pvalue.pcc.mt <- reshape(tmp.pvalue.pcc, idvar = "tmr", timevar = "cell", direction = "wide") %>% as.data.frame()
row.names(pvalue.pcc.mt) <- pvalue.pcc.mt$tmr; pvalue.pcc.mt <- pvalue.pcc.mt[,-1]
colnames(pvalue.pcc.mt) <- c("Marcophage", "Microglia", "Monocyte", "Neutrophil", "T cells", "NK cells", "NKT cells", "B cells", "Oligodendrocyte" )

library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(seq(0.6, -0.6, -0.05), colorRampPalette(c('#e66101','#fdb863','#f7f7f7','#b2abd2','#5e3c99'))(25))

ht_inter_pcc <- Heatmap(as.matrix(pcc.mt), 
                        name = "PCC", cluster_rows = T, cluster_columns = T, show_column_names = T, show_row_names = T,
                        col = col_fun, #rect_gp = gpar(),
                        #column_names_side = "top", column_names_gp = gpar(fontsize = 12), column_names_rot = 90,
                        width = ncol(pcc.mt) * unit(8, "mm"), height = nrow(pcc.mt) * unit(8, "mm"),
                        #row_split = df_IDHmut_path$NMF_Cluster, 
                        row_title = " ", 
                        #left_annotation = hl_IDHmut, right_annotation = path_anno, top_annotation = ht,
                        show_heatmap_legend = T, 
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.rect(x = x, y = y, width = width, height = height, 
                                    gp = gpar(col = "grey", fill = NA, lwd = 1.5))
                          #make the * middle
                          gb = textGrob("*")
                          gb_w = convertWidth(grobWidth(gb), "mm")
                          gb_h = convertHeight(grobHeight(gb), "mm")
                          if(pvalue.pcc.mt[i, j] < 0.05) {
                            grid.text("**", x, y - gb_h*0.5 + gb_w*0.4)
                          } else if(pvalue.pcc.mt[i, j] < 0.1) {
                            grid.text("*", x, y - gb_h*0.5 + gb_w*0.4)
                          }
                        },
                        border = T)
ht_inter_pcc
pdf("./figures/Fig3/1128_18MutA_tmr_TME_Interactions_PCC.pdf", width = 6, height = 7)
draw(ht_inter_pcc, heatmap_legend_side = "right", 
     annotation_legend_side = "right")
dev.off()



tmp.scc <- res_cor_plt[, c("tmr", "cell", "scc")]
scc.mt <- reshape(tmp.scc, idvar = "tmr", timevar = "cell", direction = "wide") %>% as.data.frame()
row.names(scc.mt) <- scc.mt$tmr; scc.mt <- scc.mt[,-1]
colnames(scc.mt) <- c( "Marcophage", "Microglia", "Monocyte", "Neutrophil", "T cells", "NK cells", "NKT cells", "B cells", "Oligodendrocyte" )

tmp.pvalue.scc <- res_cor_plt[, c("tmr", "cell", "pvalue.scc")]
pvalue.scc.mt <- reshape(tmp.pvalue.scc, idvar = "tmr", timevar = "cell", direction = "wide") %>% as.data.frame()
row.names(pvalue.scc.mt) <- pvalue.scc.mt$tmr; pvalue.scc.mt <- pvalue.scc.mt[,-1]
colnames(pvalue.scc.mt) <- c( "Marcophage", "Microglia", "Monocyte", "Neutrophil", "T cells", "NK cells", "NKT cells", "B cells", "Oligodendrocyte" )

library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(seq(0.6, -0.6, -0.05), colorRampPalette(c('#e66101','#fdb863','#f7f7f7','#b2abd2','#5e3c99'))(25))
col_fun = colorRamp2(seq(0.6, -0.6, -0.05), colorRampPalette(c('#ca0020','#f4a582','#f7f7f7','#92c5de','#0571b0'))(25))
col_fun = colorRamp2(seq(-0.6, 0.6, 0.05), colorRampPalette(c('#29419E','#698DC9',"grey95",'#F67172','#DC2B18'))(25))


ht_inter_scc <- Heatmap(as.matrix(scc.mt), 
                        name = "SCC", cluster_rows = T, cluster_columns = T, show_column_names = T, show_row_names = T,
                        col = col_fun, #rect_gp = gpar(),
                        #column_names_side = "top", column_names_gp = gpar(fontsize = 12), column_names_rot = 90,
                        width = ncol(scc.mt) * unit(8, "mm"), height = nrow(scc.mt) * unit(8, "mm"),
                        #row_split = df_IDHmut_path$NMF_Cluster, 
                        row_title = " ", 
                        #left_annotation = hl_IDHmut, right_annotation = path_anno, top_annotation = ht,
                        show_heatmap_legend = T, 
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.rect(x = x, y = y, width = width, height = height, 
                                    gp = gpar(col = "grey", fill = NA, lwd = 1.5))
                          #make the * middle
                          gb = textGrob("*")
                          gb_w = convertWidth(grobWidth(gb), "mm")
                          gb_h = convertHeight(grobHeight(gb), "mm")
                          if(pvalue.scc.mt[i, j] < 0.05) {
                            grid.text("**", x, y - gb_h*0.5 + gb_w*0.4)
                          } else if(pvalue.scc.mt[i, j] < 0.1) {
                            grid.text("*", x, y - gb_h*0.5 + gb_w*0.4)
                          }
                        },
                        border = T)
ht_inter_scc
pdf("./figures/Fig3/1028_18MutA_tmr_TME_Interactions_scc.pdf", width = 6, height = 6)
draw(ht_inter_scc, heatmap_legend_side = "right", 
     annotation_legend_side = "right")
dev.off()

plot(frac.tmr.filter$`MES-like`, frac.cell.all$MDM)

