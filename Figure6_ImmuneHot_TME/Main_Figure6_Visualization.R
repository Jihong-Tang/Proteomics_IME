#### Visualization codes for Figure 6
#### Spatial lymphocytes infiltration in the IDHm-IME gliomas
#### Author: Jihong TANG; Jiguang WANG

# Fig 6E IHC ----
library(tidyverse)
library(ggbeeswarm)
## CD3+ ---- 
dd <- readxl::read_xlsx("data/0421_IHC_RawData.xlsx", sheet = "CD3Raw")
tt <- wilcox.test(dd$Number[dd$Group == "GD"], dd$Number[dd$Group == "non-GD"], alternative = "two.sided", paired = FALSE, var.equal = FALSE)
tt$p.value

plt_pn <- dd
plt_pn$group <- ifelse(plt_pn$Group == "GD", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=Number)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_quasirandom(data=plt_pn,aes(x=group,y=Number, color=group, fill=group),width = 0.25,size= 0.1,alpha=0.7,stroke=0.8, varwidth = T)+
  geom_boxplot(data=plt_pn,aes(x=group,y=Number),width = 0.4,size=0.5,fill="transparent", outlier.colour = NA)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-20, 700),breaks = seq(0,600,150)) 
NAN_plot <- NAN_plot + scale_x_discrete(limits = c("WT", "Mut"), labels = c("Mut", "WT") )
NAN_plot<- NAN_plot +ylab(paste0("CD3+ Number")) +xlab(NULL) + 
  coord_flip()
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./CD3_IHC.pdf"), plot=figure_2,bg = 'white', width =11, height = 5.5, units = 'cm', dpi = 600)

## CD20+ ---- 
dd <- readxl::read_xlsx("data/0421_IHC_RawData.xlsx", sheet = "CD20Raw")
tt <- wilcox.test(dd$Number[dd$Group == "GD"], dd$Number[dd$Group == "non-GD"], alternative = "two.sided", paired = FALSE, var.equal = FALSE)
tt$p.value
plt_pn <- dd
plt_pn$group <- ifelse(plt_pn$Group == "GD", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=Number)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_quasirandom(data=plt_pn,aes(x=group,y=Number, color=group, fill=group),width = 0.25,size= 0.1,alpha=0.7,stroke=0.8, varwidth = T)+
  geom_boxplot(data=plt_pn,aes(x=group,y=Number),width = 0.4,size=0.5,fill="transparent", outlier.colour = NA)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-20, 320),breaks = seq(0,600,100)) 
NAN_plot <- NAN_plot + scale_x_discrete(limits = c("WT", "Mut"), labels = c("Mut", "WT") )
NAN_plot<- NAN_plot +ylab(paste0("CD20+ Number")) +xlab(NULL) + 
  coord_flip()
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file= paste0("./figures/CD20_IHC.pdf"), plot=figure_2,bg = 'white', width =11, height = 5.5, units = 'cm', dpi = 600)

## CD8+ ---- 
dd <- readxl::read_xlsx("data/0421_IHC_RawData.xlsx", sheet = "CD8Raw")
tt <- wilcox.test(dd$Number[dd$Group == "GD"], dd$Number[dd$Group == "non-GD"], alternative = "two.sided", paired = FALSE, var.equal = FALSE)
tt$p.value
plt_pn <- dd
plt_pn$group <- ifelse(plt_pn$Group == "GD", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=Number)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_quasirandom(data=plt_pn,aes(x=group,y=Number, color=group, fill=group),width = 0.25,size= 0.1,alpha=0.7,stroke=0.8, varwidth = T)+
  geom_boxplot(data=plt_pn,aes(x=group,y=Number),width = 0.4,size=0.5,fill="transparent", outlier.colour = NA)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-20, 300),breaks = seq(0,600,100)) 
NAN_plot <- NAN_plot + scale_x_discrete(limits = c("WT", "Mut"), labels = c("Mut", "WT") )
NAN_plot<- NAN_plot +ylab(paste0("CD8+ Number")) +xlab(NULL) + 
  coord_flip()
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file= paste0("./figures/CD8_IHC.pdf"), plot=figure_2,bg = 'white', width =10, height = 5.5, units = 'cm', dpi = 600)


## CD4+ ---- 
dd <- readxl::read_xlsx("data/0421_IHC_RawData.xlsx", sheet = "CD4Raw")
tt <- wilcox.test(dd$Number[dd$Group == "GD"], dd$Number[dd$Group == "non-GD"], alternative = "two.sided", paired = FALSE, var.equal = FALSE)
tt$p.value
plt_pn <- dd
plt_pn$group <- ifelse(plt_pn$Group == "GD", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=Number)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_quasirandom(data=plt_pn,aes(x=group,y=Number, color=group, fill=group),width = 0.25,size= 0.1,alpha=0.7,stroke=0.8, varwidth = T)+
  geom_boxplot(data=plt_pn,aes(x=group,y=Number),width = 0.4,size=0.5,fill="transparent", outlier.colour = NA)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-20, 240),breaks = seq(0,600,50)) 
NAN_plot <- NAN_plot + scale_x_discrete(limits = c("Mut", "WT"), labels = c("Mut", "WT") )
NAN_plot<- NAN_plot +ylab(paste0("CD4+ Number")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3")) + 
  coord_flip()
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/CD4_IHC.pdf"), plot=figure_2,bg = 'white', width =11, height = 5.5, units = 'cm', dpi = 600)

## CD38+ ---- 
library(tidyverse)
library(ggbeeswarm)
dd <- readxl::read_xlsx("data/0421_IHC_RawData.xlsx", sheet = "CD38Raw")
tt <- wilcox.test(dd$Number[dd$Group == "GD"], dd$Number[dd$Group == "non-GD"], alternative = "two.sided", paired = FALSE, var.equal = FALSE)
tt$p.value
plt_pn <- dd
plt_pn$group <- ifelse(plt_pn$Group == "GD", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=Number)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_quasirandom(data=plt_pn,aes(x=group,y=Number, color=group, fill=group),width = 0.25,size= 0.2,alpha=0.7,stroke=0.8, varwidth = T)+
  geom_boxplot(data=plt_pn,aes(x=group,y=Number),width = 0.4,size=0.5,fill="transparent", outlier.colour = NA)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-5, 100),breaks = seq(0,600,20)) 
NAN_plot <- NAN_plot + scale_x_discrete(limits = c("Mut", "WT"), labels = c("Mut", "WT") )
NAN_plot<- NAN_plot +ylab(paste0("CD38+ Number")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3")) + 
  coord_flip()
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file= paste0("./figures/CD38_IHC.pdf"), plot=figure_2,bg = 'white', width =11, height = 6.5, units = 'cm', dpi = 600)

## CD138+ ---- 
library(tidyverse)
library(ggbeeswarm)
dd <- readxl::read_xlsx("data/0421_IHC_RawData.xlsx", sheet = "CD138Raw")
tt <- wilcox.test(dd$Number[dd$Group == "GD"], dd$Number[dd$Group == "non-GD"], alternative = "two.sided", paired = FALSE, var.equal = FALSE)
tt$p.value
plt_pn <- dd
plt_pn$group <- ifelse(plt_pn$Group == "GD", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=Number)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_quasirandom(data=plt_pn,aes(x=group,y=Number, color=group, fill=group),width = 0.25,size= 0.2,alpha=0.7,stroke=0.8, varwidth = T)+
  geom_boxplot(data=plt_pn,aes(x=group,y=Number),width = 0.4,size=0.5,fill="transparent", outlier.colour = NA)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2, 28),breaks = seq(0,600,5)) 
NAN_plot <- NAN_plot + scale_x_discrete(limits = c("Mut", "WT"), labels = c("Mut", "WT") )
NAN_plot<- NAN_plot +ylab(paste0("CD138+ Number")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3")) + 
  coord_flip()
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file= paste0("./figures/CD138_IHC.pdf"), plot=figure_2,bg = 'white', width =11, height = 5.5, units = 'cm', dpi = 600)

# Figure 6F spatial lymphocyte statistics ----
library(tidyverse)
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}

specie <- c(rep("Near", 2), rep("Middle", 2), rep("Far", 2))
condition <- rep(c("CD4", "CD8"), 3)
# value <- c(21, 520-21, 62, 283-62)

value <- c(97, 79, 26, 37, 47, 97)
data <- data.frame(specie, condition, value)

# Stacked + percent
ggplot(data, aes(fill = condition, y = value, x = specie)) +
  geom_bar(position = "fill", stat = "identity", show.legend = F, width = .7) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.2), labels = scales::percent) +
  scale_x_discrete(limits = c("Near", "Middle", "Far")) +
  scale_fill_manual(values = c(
    "CD4" = "#e41a1c",
    "CD8" = "#fffe54"
  )) +
  labs(x = "", y = "Percentage") +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "transparent", color = "transparent"), plot.margin = unit(c(2, 2, 2, 2), "lines"),
    plot.title = element_text(size = 34, vjust = 0.5, hjust = 0.5, face = "bold.italic", color = "transparent"), text = element_text(size = 14, face = "bold"),
    legend.key.width = unit(0.6, "cm"), legend.key.height = unit(0.6, "cm"), legend.position = "none", legend.text = element_text(size = 14, hjust = 0, face = "bold"),
    axis.text.y = element_blank(), axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_text(size = 14, face = "plain", color = "black"), axis.title.y = element_text(size = 14, face = "plain", color = "black")
  ) +
  coord_flip()

ggsave("figures/Fig6F_PLC_region_CD4CD8_comparison.pdf", width = 4, height = 2.5)
M <- as.table(rbind(c(97, 26, 47), c(176, 63, 144)))
chisq.test(M)

# Fig6G Dimplot ----
obj.T.no <- readRDS("./data/0704_obj_T_noUnknown.rds")
color_scheme2 <- c(
  "CD4.Treg"   = "#e41a1c",
  "CD4.ISG+T"  = "#ffff33",
  "CD4.Tn"     = "#ff7f00",
  "CD4.Tcm"    = "#377eb8",
  "CD4.Tem"    = "#4daf4a",
  "CD8.Trm"    = "#984ea3",
  "CD8.NKT"    = "#999999",
  "CD8.Tex"    = "#a65628",
  "CD8.Tem"    = "#f781bf",
  "Unknown"    = "#ffffff"
)
DimPlot(
  object = obj.T.no, reduction = "umap",
  group.by = "T_cluster",
  cols = color_scheme2,
  label = T,
  pt.size = 1, # order = order,
  label.size = 0
) & NoAxes() & NoLegend()
ggsave("./figures/Fig6G_sc.Tdetails.DimPlot.pdf", width = 5.5, height = 5)

# Fig6H Barplots ----
xxx <- read_delim("data/scRNA_cluster.txt")
obj.T.no$pc <- ifelse(obj.T.no$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "IME"], "IME", "others")
table(obj.T.no$pc, obj.T.no$T_cluster)
table(obj.T.no$pc)
# prop.test for CD8.Tex in two clusters
prop.test(x = c(329, 73), n = c(1194, 528), alternative = "two.sided")
prop.test(x = c(280, 51), n = c(1194, 528), alternative = "two.sided")
p.test <- prop.test(x = c(74, 126), n = c(1194, 528), alternative = "two.sided")

count.data <- data.frame(
  Var1 = c(
    "2CD4.ISG+T", "4CD4.Tcm", "5CD4.Tem", "3CD4.Tn",
    "1CD4.Treg", "9CD8.NKT", "7CD8.Tem", "8CD8.Tex", "6CD8.Trm"
  ),
  Freq = c(21, 92, 88, 132, 41, 74, 280, 329, 137)
)
count.data <- count.data[order(-count.data$Freq), ]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.25 * Freq) %>%
  mutate(percent = Freq / sum(count.data$Freq))
count.data

# order <- c("CD4T","CD56bright_NK","CD56dim_NK", "CD8T", "B", "NKT",  "PC","Treg" )
ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size = 0.5, stat = "identity", color = "white", alpha = 0.85) +
  scale_fill_manual(values = c(
    "1CD4.Treg" = "#e41a1c",
    "2CD4.ISG+T" = "#ffff33",
    "3CD4.Tn" = "#ff7f00",
    "4CD4.Tcm" = "#377eb8",
    "5CD4.Tem" = "#4daf4a",
    "6CD8.Trm" = "#984ea3",
    "9CD8.NKT" = "#999999",
    "8CD8.Tex" = "#a65628",
    "7CD8.Tem" = "#f781bf"
  )) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_blank(),
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.key.width = unit(0.6, "cm"), legend.key.height = unit(0.6, "cm"), legend.key.size = unit(5, "lines"), legend.key = element_rect(size = 0.1, color = NA),
    legend.position = "right", legend.text = element_text(size = 14, face = "bold.italic"), legend.margin = margin(t = 0.1, r = 0.1, b = 0, l = 0.1, unit = "cm")
  )

ggsave("./figures/Fig6H_sc.Tcell.IME.bar.pdf", width = 3, height = 4)

count.data <- data.frame(
  Var1 = c(
    "2CD4.ISG+T", "4CD4.Tcm", "5CD4.Tem", "3CD4.Tn",
    "1CD4.Treg", "9CD8.NKT", "7CD8.Tem", "8CD8.Tex", "6CD8.Trm"
  ),
  Freq = c(3, 64, 58, 87, 17, 126, 51, 73, 49)
)
count.data <- count.data[order(-count.data$Freq), ]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.25 * Freq) %>%
  mutate(percent = Freq / sum(count.data$Freq))
count.data

ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size = 0.5, stat = "identity", color = "white", alpha = 0.85) +
  scale_fill_manual(values = c(
    "1CD4.Treg" = "#e41a1c",
    "2CD4.ISG+T" = "#ffff33",
    "3CD4.Tn" = "#ff7f00",
    "4CD4.Tcm" = "#377eb8",
    "5CD4.Tem" = "#4daf4a",
    "6CD8.Trm" = "#984ea3",
    "9CD8.NKT" = "#999999",
    "8CD8.Tex" = "#a65628",
    "7CD8.Tem" = "#f781bf"
  )) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_blank(),
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.key.width = unit(0.6, "cm"), legend.key.height = unit(0.6, "cm"), legend.key.size = unit(5, "lines"), legend.key = element_rect(size = 0.1, color = NA),
    legend.position = "right", legend.text = element_text(size = 14, face = "bold.italic"), legend.margin = margin(t = 0.1, r = 0.1, b = 0, l = 0.1, unit = "cm")
  )
ggsave("./figures/Fig6H_sc.Tcell.other.bar.pdf", width = 3, height = 4)

# FIg6I T cell state plots ----
gns_naive <- c("CCR7", "SELL", "LEF1", "TCF7")
gns_effector <- c("IFNG", "PRF1", "GNLY", "NKG7", "GZMA", "GZMB", "GZMK", "GZMH", "CX3CR1", "CST7", "CTSW")
gns_exhausted <- c("PDCD1", "LAG3", "TIGIT", "CTLA4", "HAVCR2", "TOX", "CXCL13", "ENTPD1")
list_Tcell_states <- list(gns_naive, gns_effector, gns_exhausted)
obj.T.no <- AddModuleScore(obj.T.no, features = list_Tcell_states, name = c("Naive", "Effector", "Exhausted"))
df_Tstate <- obj.T.no@meta.data[, c("T_cluster", "pc", "Naive1", "Effector2", "Exhausted3")]

library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)

# Plot
ggplot(df_Tstate, aes(x = Exhausted3, y = pc, fill = pc)) +
  # geom_density_ridges(alpha=0.6, stat="binline", bins=30) +
  geom_density_ridges(alpha = 1) +
  scale_y_discrete(limits = c("others", "IME")) +
  scale_fill_manual(values = c("#f9da56", "#8b351f")) +
  theme_ridges() +
  theme(
    legend.position = "none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )
ggsave("figures/FIg6I_Tcell_Exhausted3_ridgeline.pdf", width = 6, height = 4)

library(ggpubr)
ggecdf(df_Tstate,
  x = "Exhausted3", color = "pc",
  size = 0.8, alpha = 0.3
) +
  geom_hline(yintercept = c(0, 1), linetype = "dashed", color = "#dddddd") +
  scale_color_manual(values = c("#f9da56", "#8b351f")) +
  theme_classic() +
  theme( # panel.background=element_rect(fill='transparent',color='white'),
    # text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),
    legend.position = "none",
    legend.key.width = unit(.5, "cm"), legend.key.height = unit(0.5, "cm"),
    legend.text = element_text(size = 12, face = "plain"),
    axis.text.y = element_text(size = 16, vjust = 0.5, hjust = 1, face = "plain", color = "black"),
    legend.title = element_text(size = 12, vjust = 0.5, hjust = 0, face = "plain"),
    axis.text.x = element_text(size = 16, angle = 0, face = "plain", color = "black"),
    axis.title.x = element_text(size = 20, vjust = 0, hjust = 0.5, face = "plain", color = "black"),
    axis.title.y = element_text(size = 20, hjust = 0.5, vjust = 2, face = "plain", color = "black")
  )
ggsave("figures/Fig6I_Tcell_Exhausted3_cdf.pdf", width = 4, height = 3)