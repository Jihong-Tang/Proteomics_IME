library(tidyverse)
library(ggpubr)
library(preprocessCore)
library(ggbeeswarm)

# 1- Load all the data ---- 
## 1.1 TCGA data ----
tcga_glioma <- read_delim("../../0427_TMEsubtypes_Characterize_RNA/data/TCGAGBMLGG.genexpression.RPKM.txt") %>% as.data.frame()
tcga_glioma <- tcga_glioma[!duplicated(tcga_glioma$gene_name),]
row.names(tcga_glioma) <- tcga_glioma$gene_name
tcga_glioma <- tcga_glioma[, -c(1,2)]
tcga_P <- tcga_glioma[,substr(names(tcga_glioma),14,15)=='01']
tcga_P <- tcga_glioma

tcga_cl <- readxl::read_xlsx("../../0427_TMEsubtypes_Characterize_RNA/data/TCGA_GBMLGG_Clinical_Cell2016.xlsx")

sum(tcga_cl$RNAseq == "Yes")

tcga_idhmutA <- tcga_cl$Case[which(tcga_cl$`IDH/codel subtype`=="IDHmut-non-codel")]
dfall = tcga_P#backup the whole data set,we only need the IDHmut for analysis

df.mutA.tcga = dfall[,substr(names(dfall),1,12) %in% tcga_idhmutA]

df.mutA.tcga_vst <- log2(df.mutA.tcga + 1)
mt_lgg_qn <- normalize.quantiles(as.matrix(df.mutA.tcga_vst), copy = T, keep.names = T)
df_tcga_lgg <- data.frame((mt_lgg_qn), check.names = F)

tmp_scale <- apply(mt_lgg_qn, 1, scale)
df.mutA.tcga_scale <- data.frame(t(tmp_scale), check.names = F)
colnames(df.mutA.tcga_scale) <- colnames(df.mutA.tcga_vst)

nmf.tcga <- readxl::read_xlsx("data/IME_Signatures/nmf_lgg_res_1123.xlsx")

## 1.2 CGGA325 CGGA 693 data ---- 
cgga325 <- read.table("../../0427_TMEsubtypes_Characterize_RNA//data/CGGA/CGGA.mRNAseq_325.RSEM-genes.20200506.txt", header = T)
row.names(cgga325) <- cgga325$Gene_Name
cl.cgga325 <- read.delim("../../0427_TMEsubtypes_Characterize_RNA//data/CGGA/CGGA.mRNAseq_325_clinical.20200506.txt", check.names = F)
colnames(cl.cgga325)[8] <- "Censor"

cl.cgga325.update <- read.delim("../../!IMPORTANT_FILES/clinical_Info_update_20231120/CGGA.mRNAseq_325_clinical.20231120.txt", check.names = F)
colnames(cl.cgga325.update)[8] <- "Censor"

id_mutA <- cl.cgga325$CGGA_ID[cl.cgga325$IDH_mutation_status == "Mutant" & cl.cgga325[12] == "Non-codel"]

cgga325_mutA <- cgga325[, id_mutA[!is.na(id_mutA)]]

cgga693 <- read.table("../../0427_TMEsubtypes_Characterize_RNA//data/CGGA/CGGA.mRNAseq_693.RSEM-genes.20200506.txt", header = T)
row.names(cgga693) <- cgga693$Gene_Name
cl.cgga693 <- read.delim("../../0427_TMEsubtypes_Characterize_RNA//data/CGGA/CGGA.mRNAseq_693_clinical.20200506.txt", check.names = F)
colnames(cl.cgga693)[8] <- "Censor"

cl.cgga693.update <- read.delim("../../!IMPORTANT_FILES/clinical_Info_update_20231120/CGGA.mRNAseq_693_clinical.20231120.txt", check.names = F)
colnames(cl.cgga693.update)[8] <- "Censor"

id_mutA <- cl.cgga693$CGGA_ID[cl.cgga693$IDH_mutation_status == "Mutant" & cl.cgga693[12] == "Non-codel"]
cgga693_mutA <- cgga693[, id_mutA[!is.na(id_mutA)]]

df_cgga325_vst <- log2(cgga325_mutA + 1)
mt_cgga325_qn <- normalize.quantiles(as.matrix(df_cgga325_vst), copy = T, keep.names = T)
df_cgga325 <- data.frame((mt_cgga325_qn), check.names = F)

tmp_scale <- apply(mt_cgga325_qn, 1, scale)
df_cgga325_scale <- data.frame(t(tmp_scale), check.names = F)
colnames(df_cgga325_scale) <- colnames(df_cgga325_vst)

df_cgga693_vst <- log2(cgga693_mutA + 1)
mt_cgga693_qn <- normalize.quantiles(as.matrix(df_cgga693_vst), copy = T, keep.names = T)
df_cgga693 <- data.frame((mt_cgga693_qn), check.names = F)

tmp_scale <- apply(mt_cgga693_qn, 1, scale)
df_cgga693_scale <- data.frame(t(tmp_scale), check.names = F)
colnames(df_cgga693_scale) <- colnames(df_cgga693_vst)

nmf.cgga325 <- read_delim("data/IME_Signatures/CGGA325_proteosubgroup_Vali_1121.txt")
nmf.cgga693 <- read_delim("data/IME_Signatures/CGGA693_proteosubgroup_Vali_1121.txt", delim = "\t")

# 1.3 CGPA protein data ---- 
pn_cgpa <- read.table("data/IME_Signatures/CGGA_7kProteins_Imputeall_Feb23.txt",
                      check.names = F)

nmf.cgpa <- readxl::read_xlsx("data/IME_Signatures/nmf_subtypes_clinical_updated_23Dec03.xlsx")

# 2- Protein Expression Data ---- 
## IGHG1 ----
pn <- "IGHG1"
dd.pn <- data.frame(Cohort_ID = rownames(pn_cgpa), pn = pn, pn_val = pn_cgpa[, pn])

plt.pn <- merge(dd.pn, nmf.cgpa, by = "Cohort_ID")
plt.pn$group <- ifelse(plt.pn$NMF_Cluster == "ProteoNMF3", "IME", "Non-IME")
plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "greater")

plt_pn <- plt.pn

plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2.5,2),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Abundance")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_CGPA_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)

## IGHG2 ----
pn <- "IGHG2"
dd.pn <- data.frame(Cohort_ID = rownames(pn_cgpa), pn = pn, pn_val = pn_cgpa[, pn])

plt.pn <- merge(dd.pn, nmf.cgpa, by = "Cohort_ID")
plt.pn$group <- ifelse(plt.pn$NMF_Cluster == "ProteoNMF3", "IME", "Non-IME")
plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "greater")

plt_pn <- plt.pn

plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2.5,2),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Abundance")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_CGPA_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## IGHG3 ----
pn <- "IGHG3"
dd.pn <- data.frame(Cohort_ID = rownames(pn_cgpa), pn = pn, pn_val = pn_cgpa[, pn])

plt.pn <- merge(dd.pn, nmf.cgpa, by = "Cohort_ID")
plt.pn$group <- ifelse(plt.pn$NMF_Cluster == "ProteoNMF3", "IME", "Non-IME")
plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "greater")

plt_pn <- plt.pn

plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2.2,2.5),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Abundance")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_CGPA_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## FCGR2A ----
pn <- "FCGR2A"
dd.pn <- data.frame(Cohort_ID = rownames(pn_cgpa), pn = pn, pn_val = pn_cgpa[, pn])

plt.pn <- merge(dd.pn, nmf.cgpa, by = "Cohort_ID")
plt.pn$group <- ifelse(plt.pn$NMF_Cluster == "ProteoNMF3", "IME", "Non-IME")
plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "greater")

plt_pn <- plt.pn
plt_pn <- plt_pn[plt_pn$Cohort_ID != "CGGA_P778",]
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2, 2),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Abundance")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_CGPA_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## SERPINA3 ----
pn <- "SERPINA3"
dd.pn <- data.frame(Cohort_ID = rownames(pn_cgpa), pn = pn, pn_val = pn_cgpa[, pn])

plt.pn <- merge(dd.pn, nmf.cgpa, by = "Cohort_ID")
plt.pn$group <- ifelse(plt.pn$NMF_Cluster == "ProteoNMF3", "IME", "Non-IME")
plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "greater")

plt_pn <- plt.pn
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2.2, 2.2),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Abundance")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_CGPA_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## IGKC ----
pn <- "IGKC"
dd.pn <- data.frame(Cohort_ID = rownames(pn_cgpa), pn = pn, pn_val = pn_cgpa[, pn])

plt.pn <- merge(dd.pn, nmf.cgpa, by = "Cohort_ID")
plt.pn$group <- ifelse(plt.pn$NMF_Cluster == "ProteoNMF3", "IME", "Non-IME")
plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "greater")

plt_pn <- plt.pn
#plt_pn <- plt_pn[plt_pn$Cohort_ID != "CGGA_P778",]
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2.2, 2.2),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Abundance")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_CGPA_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## GBP1 ----
pn <- "GBP1"
dd.pn <- data.frame(Cohort_ID = rownames(pn_cgpa), pn = pn, pn_val = pn_cgpa[, pn])

plt.pn <- merge(dd.pn, nmf.cgpa, by = "Cohort_ID")
plt.pn$group <- ifelse(plt.pn$NMF_Cluster == "ProteoNMF3", "IME", "Non-IME")
plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "greater")

plt_pn <- plt.pn
#plt_pn <- plt_pn[plt_pn$Cohort_ID != "CGGA_P778",]
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2.2, 3),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Abundance")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_CGPA_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## GBP2 ----
pn <- "GBP2"
dd.pn <- data.frame(Cohort_ID = rownames(pn_cgpa), pn = pn, pn_val = pn_cgpa[, pn])

plt.pn <- merge(dd.pn, nmf.cgpa, by = "Cohort_ID")
plt.pn$group <- ifelse(plt.pn$NMF_Cluster == "ProteoNMF3", "IME", "Non-IME")
plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "greater")

plt_pn <- plt.pn
#plt_pn <- plt_pn[plt_pn$Cohort_ID != "CGGA_P778",]
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2.2, 3.3),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Abundance")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_CGPA_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)

# 3- Gene Expression Data ---- 
# CGGA693 ----
## IGHG1 ----
pn <- "IGHG1"
dd.pn <- data.frame(CGGA_ID = colnames(df_cgga693_scale), pn_val = df_cgga693_scale[pn,] %>% as.numeric())

plt.pn <- merge(dd.pn, nmf.cgga693, by = "CGGA_ID")
plt.pn <- plt.pn[plt.pn$nmf_subtype != "Mixed",]
plt.pn$group <- ifelse(plt.pn$nmf_subtype == "ProteoNMF3", "IME", "Non-IME")
#plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "greater")

plt_pn <- plt.pn
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-1.2, 3),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Expression")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_CGGA693_gene_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)

## IGHG2 ----
pn <- "IGHG2"
dd.pn <- data.frame(CGGA_ID = colnames(df_cgga693_scale), pn_val = df_cgga693_scale[pn,] %>% as.numeric())

plt.pn <- merge(dd.pn, nmf.cgga693, by = "CGGA_ID")
plt.pn <- plt.pn[plt.pn$nmf_subtype != "Mixed",]
plt.pn$group <- ifelse(plt.pn$nmf_subtype == "ProteoNMF3", "IME", "Non-IME")
#plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "greater")

plt_pn <- plt.pn
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-1, 3.6),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Expression")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_CGGA693_gene_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## IGHG3 ----
pn <- "IGHG3"
dd.pn <- data.frame(CGGA_ID = colnames(df_cgga693_scale), pn_val = df_cgga693_scale[pn,] %>% as.numeric())

plt.pn <- merge(dd.pn, nmf.cgga693, by = "CGGA_ID")
plt.pn <- plt.pn[plt.pn$nmf_subtype != "Mixed",]
plt.pn$group <- ifelse(plt.pn$nmf_subtype == "ProteoNMF3", "IME", "Non-IME")
#plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "greater")

plt_pn <- plt.pn
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-1, 3.7),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Expression")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_CGGA693_gene_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## FCGR2A ----
pn <- "FCGR2A"
dd.pn <- data.frame(CGGA_ID = colnames(df_cgga693_scale), pn_val = df_cgga693_scale[pn,] %>% as.numeric())

plt.pn <- merge(dd.pn, nmf.cgga693, by = "CGGA_ID")
plt.pn <- plt.pn[plt.pn$nmf_subtype != "Mixed",]
plt.pn$group <- ifelse(plt.pn$nmf_subtype == "ProteoNMF3", "IME", "Non-IME")
#plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "greater")

plt_pn <- plt.pn
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2.1, 3),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Expression")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_CGGA693_gene_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## SERPINA3 ----
pn <- "SERPINA3"
dd.pn <- data.frame(CGGA_ID = colnames(df_cgga693_scale), pn_val = df_cgga693_scale[pn,] %>% as.numeric())

plt.pn <- merge(dd.pn, nmf.cgga693, by = "CGGA_ID")
plt.pn <- plt.pn[plt.pn$nmf_subtype != "Mixed",]
plt.pn$group <- ifelse(plt.pn$nmf_subtype == "ProteoNMF3", "IME", "Non-IME")
#plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "greater")

plt_pn <- plt.pn
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2.2, 3),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Expression")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_CGGA693_gene_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## CCL2 ----
pn <- "CCL2"
dd.pn <- data.frame(CGGA_ID = colnames(df_cgga693_scale), pn_val = df_cgga693_scale[pn,] %>% as.numeric())

plt.pn <- merge(dd.pn, nmf.cgga693, by = "CGGA_ID")
plt.pn <- plt.pn[plt.pn$nmf_subtype != "Mixed",]
plt.pn$group <- ifelse(plt.pn$nmf_subtype == "ProteoNMF3", "IME", "Non-IME")
#plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "two.sided")

plt_pn <- plt.pn
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2.2, 3.5),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Expression")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_CGGA693_gene_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## CCR2 ----
pn <- "CCR2"
dd.pn <- data.frame(CGGA_ID = colnames(df_cgga693_scale), pn_val = df_cgga693_scale[pn,] %>% as.numeric())

plt.pn <- merge(dd.pn, nmf.cgga693, by = "CGGA_ID")
plt.pn <- plt.pn[plt.pn$nmf_subtype != "Mixed",]
plt.pn$group <- ifelse(plt.pn$nmf_subtype == "ProteoNMF3", "IME", "Non-IME")
#plt.pn$pn_val <- scale(plt.pn$pn_val)
plt.pn <- plt.pn[plt.pn$pn_val <=3,]
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "two.sided")

plt_pn <- plt.pn
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-1, 3),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Expression")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_CGGA693_gene_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## IGKC ----
pn <- "IGKC"
dd.pn <- data.frame(CGGA_ID = colnames(df_cgga693_scale), pn_val = df_cgga693_scale[pn,] %>% as.numeric())

plt.pn <- merge(dd.pn, nmf.cgga693, by = "CGGA_ID")
plt.pn <- plt.pn[plt.pn$nmf_subtype != "Mixed",]
plt.pn$group <- ifelse(plt.pn$nmf_subtype == "ProteoNMF3", "IME", "Non-IME")
#plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "two.sided")

plt_pn <- plt.pn
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-1.5, 3),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Expression")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_CGGA693_gene_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## GBP1 ----
pn <- "GBP1"
dd.pn <- data.frame(CGGA_ID = colnames(df_cgga693_scale), pn_val = df_cgga693_scale[pn,] %>% as.numeric())

plt.pn <- merge(dd.pn, nmf.cgga693, by = "CGGA_ID")
plt.pn <- plt.pn[plt.pn$nmf_subtype != "Mixed",]
plt.pn$group <- ifelse(plt.pn$nmf_subtype == "ProteoNMF3", "IME", "Non-IME")
#plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "two.sided")

plt_pn <- plt.pn
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-1.6, 4.5),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Expression")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_CGGA693_gene_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## GBP2 ----
pn <- "GBP2"
dd.pn <- data.frame(CGGA_ID = colnames(df_cgga693_scale), pn_val = df_cgga693_scale[pn,] %>% as.numeric())

plt.pn <- merge(dd.pn, nmf.cgga693, by = "CGGA_ID")
plt.pn <- plt.pn[plt.pn$nmf_subtype != "Mixed",]
plt.pn$group <- ifelse(plt.pn$nmf_subtype == "ProteoNMF3", "IME", "Non-IME")
#plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "two.sided")

plt_pn <- plt.pn
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2.2, 2.5),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Expression")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_CGGA693_gene_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


# tcga ----
colnames(nmf.tcga)[1] <- "CGGA_ID"
## FCGR2A ----
pn <- "FCGR2A"
dd.pn <- data.frame(CGGA_ID = substr(colnames(df.mutA.tcga_scale),1,12), pn_val = df.mutA.tcga_scale[pn,] %>% as.numeric())

plt.pn <- merge(dd.pn, nmf.tcga, by = "CGGA_ID")
plt.pn <- plt.pn[plt.pn$nmf_subtype != "Mixed",]
plt.pn$group <- ifelse(plt.pn$nmf_subtype == "ProteoNMF3", "IME", "Non-IME")
#plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "greater")

plt_pn <- plt.pn
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-1.8, 4),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Expression")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_tcga_gene_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## SERPINA3 ----
pn <- "SERPINA3"
dd.pn <- data.frame(CGGA_ID = substr(colnames(df.mutA.tcga_scale),1,12), pn_val = df.mutA.tcga_scale[pn,] %>% as.numeric())

plt.pn <- merge(dd.pn, nmf.tcga, by = "CGGA_ID")
plt.pn <- plt.pn[plt.pn$nmf_subtype != "Mixed",]
plt.pn$group <- ifelse(plt.pn$nmf_subtype == "ProteoNMF3", "IME", "Non-IME")
#plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "greater")

plt_pn <- plt.pn
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-1, 5),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Expression")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_tcga_gene_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## CCL2 ----
pn <- "CCL2"
dd.pn <- data.frame(CGGA_ID = substr(colnames(df.mutA.tcga_scale),1,12), pn_val = df.mutA.tcga_scale[pn,] %>% as.numeric())

plt.pn <- merge(dd.pn, nmf.tcga, by = "CGGA_ID")
plt.pn <- plt.pn[plt.pn$nmf_subtype != "Mixed",]
plt.pn$group <- ifelse(plt.pn$nmf_subtype == "ProteoNMF3", "IME", "Non-IME")
#plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "two.sided")

plt_pn <- plt.pn
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2.2, 4),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Expression")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_tcga_gene_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## CCR2 ----
pn <- "CCR2"
dd.pn <- data.frame(CGGA_ID = substr(colnames(df.mutA.tcga_scale),1,12), pn_val = df.mutA.tcga_scale[pn,] %>% as.numeric())

plt.pn <- merge(dd.pn, nmf.tcga, by = "CGGA_ID")
plt.pn <- plt.pn[plt.pn$nmf_subtype != "Mixed",]
plt.pn$group <- ifelse(plt.pn$nmf_subtype == "ProteoNMF3", "IME", "Non-IME")
#plt.pn$pn_val <- scale(plt.pn$pn_val)
plt.pn <- plt.pn[plt.pn$pn_val <=3,]
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "two.sided")

plt_pn <- plt.pn
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-1, 2),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Expression")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_tcga_gene_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## GBP1 ----
pn <- "GBP1"
dd.pn <- data.frame(CGGA_ID = substr(colnames(df.mutA.tcga_scale),1,12), pn_val = df.mutA.tcga_scale[pn,] %>% as.numeric())

plt.pn <- merge(dd.pn, nmf.tcga, by = "CGGA_ID")
plt.pn <- plt.pn[plt.pn$nmf_subtype != "Mixed",]
plt.pn$group <- ifelse(plt.pn$nmf_subtype == "ProteoNMF3", "IME", "Non-IME")
#plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "two.sided")

plt_pn <- plt.pn
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-1.6, 3.5),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Expression")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_tcga_gene_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)


## GBP2 ----
pn <- "GBP2"
dd.pn <- data.frame(CGGA_ID = substr(colnames(df.mutA.tcga_scale),1,12), pn_val = df.mutA.tcga_scale[pn,] %>% as.numeric())

plt.pn <- merge(dd.pn, nmf.tcga, by = "CGGA_ID")
plt.pn <- plt.pn[plt.pn$nmf_subtype != "Mixed",]
plt.pn$group <- ifelse(plt.pn$nmf_subtype == "ProteoNMF3", "IME", "Non-IME")
#plt.pn$pn_val <- scale(plt.pn$pn_val)
wilcox.test(plt.pn$pn_val[plt.pn$group == "IME"], plt.pn$pn_val[plt.pn$group == "Non-IME"], alternative = "two.sided")

plt_pn <- plt.pn
plt_pn$group <- ifelse(plt_pn$group == "IME", "Mut", "WT")
NAN_plot <- ggplot(data=plt_pn,aes(x=group,y=pn_val)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_pn,aes(x=group,y=pn_val),width = 0.4,size=0.5,fill="transparent", outlier.colour = "white")+
  geom_quasirandom(data=plt_pn,aes(x=group,y=pn_val, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-1.6, 3.5),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab(paste0(pn, " Relative Expression")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Plasma_mecha/0302_tcga_gene_IME_", pn, ".pdf"), plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)

# 4 HeatMap dotplots for visualize important genes DE results ---- 
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

de_pns <- readxl::read_xlsx("data/Fig6_mechanism_candidates.xlsx")
de_pns <- de_pns[de_pns$Cohort == "CGPA", seq(1,5)]
colnames(de_pns) <- c("Gene", "pvalue", "fc", "qvalue", "Cohort")

de_gns <- readxl::read_xlsx("data/Fig6_mechanism_candidates.xlsx", sheet = "RNA")
de_gns_cgga <- de_gns[de_gns$Cohort == "CGGA693",seq(1,5)]
colnames(de_gns_cgga) <- c("Gene", "pvalue", "fc", "qvalue", "Cohort")
de_gns_tcga <- de_gns[de_gns$Cohort == "TCGA", seq(1,5)]
colnames(de_gns_tcga) <- c("Gene", "pvalue", "fc", "qvalue", "Cohort")

gns <- c("IGHG1", "IGHG2", "IGHG3", "IGKC", "FCGR2A", "CCL2", "CCR2", "GBP1", "GBP2", "SERPINA3")

plt_pns <- de_pns[de_pns$Gene %in% gns,]
plt_gns_cgga <- de_gns_cgga[de_gns_cgga$Gene %in% gns,]
plt_gns_tcga <- de_gns_tcga[de_gns_tcga$Gene %in% gns,]

plt_all <- rbind(plt_pns, plt_gns_cgga, plt_gns_tcga)

## version 1 ----
# Transform pvalue to -log10 scale
df <- plt_all %>%
  mutate(log10_p = -log10(pvalue)) %>% 
  mutate(log10_p = pmin(-log10(pvalue), 6))  # cap at 8

# Prepare matrices for color and size
fc_matrix <- reshape2::acast(df, Cohort ~ Gene, value.var = "fc")
size_matrix <- reshape2::acast(df, Cohort ~ Gene, value.var = "log10_p")
q_matrix <- reshape2::acast(df, Cohort ~ Gene, value.var = "pvalue")

# Set color scale for fold change
col_fun <- colorRamp2(c(0, 2.5), c("white", "red"))
col_fun = colorRamp2(seq(0,2.5,0.1), colorRampPalette(c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026'))(26))
col_fun = colorRamp2(seq(0,2.5,0.1), colorRampPalette(c('#fff5f0','#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d'))(26))
# Draw heatmap using dots
ht <- Heatmap(matrix = fc_matrix,
        name = "log2FC",
        col = col_fun,
        height = unit(nrow(fc_matrix)*1.25, "cm"), width = unit(ncol(fc_matrix)*1, "cm"),
        rect_gp = gpar(col = "#969696", fill = NA, lwd = 1.5, lty = "dashed"),
        row_order = c("CGPA", "CGGA693", "TCGA"),
        column_order = c("IGHG1", "IGHG2", "IGHG3", "IGKC", "FCGR2A", "CCL2", "CCR2", "GBP1", "GBP2", "SERPINA3"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (is.na(fc_matrix[i, j]) || is.na(size_matrix[i, j])) {
            # NA case → draw small grey dot
            grid.circle(x = x, y = y,
                        r = unit(0.1, "cm"),  # small fixed radius
                        gp = gpar(fill = "#dddddd", col = NA))
          } else {
          #dot_size <- unit(size_matrix[i, j] / max(size_matrix, na.rm = TRUE) * 0.4, "cm")
          dot_radius <- scales::rescale(size_matrix[i, j],
                                        to = c(0.1, 0.4),
                                        from = range(size_matrix, na.rm = TRUE))
          dot_size <- unit(dot_radius, "cm")
          dot_color <- col_fun(fc_matrix[i, j])
          # Check if q-value is below threshold
          dot_border <- if (!is.na(q_matrix[i, j]) && q_matrix[i, j] < 0.05) "black" else NA
          
          grid.circle(x = x, y = y,
                      r = dot_size,
                      gp = gpar(fill = dot_color, col = dot_border, lwd = 1))
          }
        },
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_names_rot = 45,
        row_names_side = "left",
        column_names_gp = gpar(fontsize = 12, fontface = "bold"),
        row_names_gp = gpar(fontsize = 12, fontface = "bold"),
        border = T,
        border_gp = gpar(col = "black", lwd = 1.5)
        )
pdf("figures/Plasma_mecha/0416_genes_heatmap.pdf", width = 6, height = 4)
draw(ht)
dev.off()

## version 2 ----
# Transform values
df <- plt_all %>%
  mutate(fc_capped = pmin(fc, 3)) %>%                      # cap fold change
  mutate(log10_p = pmin(-log10(pvalue), 6))                # cap -log10(pvalue)

# Matrices
color_matrix <- reshape2::acast(df, Cohort ~ Gene, value.var = "log10_p")  # color by p-value
size_matrix  <- reshape2::acast(df, Cohort ~ Gene, value.var = "fc_capped")  # size by fc
p_matrix     <- reshape2::acast(df, Cohort ~ Gene, value.var = "pvalue")     # for borders

# Color function for p-value (white to red)
col_fun <- colorRamp2(c(0, 4), c("white", "red"))

# Heatmap
Heatmap(matrix = color_matrix,
        name = "-log10(p)",
        col = col_fun,
        height = unit(nrow(color_matrix)*1.5, "cm"),
        width = unit(ncol(color_matrix)*1, "cm"),
        rect_gp = gpar(col = "grey", fill = NA),  # keep border, remove fill
        cell_fun = function(j, i, x, y, width, height, fill) {
          # Size is based on capped FC
          dot_size <- unit(size_matrix[i, j] / max(size_matrix, na.rm = TRUE) * 0.4, "cm")
          # Color is based on -log10(p)
          dot_color <- col_fun(color_matrix[i, j])
          # Border if significant
          dot_border <- if (!is.na(p_matrix[i, j]) && p_matrix[i, j] < 0.05) "black" else NA
          
          grid.circle(x = x, y = y,
                      r = dot_size,
                      gp = gpar(fill = dot_color, col = dot_border, lwd = 1))
        },
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_names_rot = 45,
        row_names_side = "left",
        column_names_gp = gpar(fontsize = 12, fontface = "bold"),
        row_names_gp = gpar(fontsize = 12, fontface = "bold"))


# 5- scRNA analysis ---- 
library(Seurat)
obj.lymphocyte <- readRDS("data/scRNA_revise/0422_merged_lymphocyte_annotated.rds")

lymphocyte.noTmr <- subset(obj.lymphocyte, subset = cluster != "tumor_like")
lymphocyte.noTmr <- subset(lymphocyte.noTmr, subset = cluster != "myeloid_like")
lymphocyte.noTmr <- subset(lymphocyte.noTmr, subset = cluster != "Unknown")

FeaturePlot(lymphocyte.noTmr, features = "IGKC", reduction = "tsne")
ttt <- FetchData(lymphocyte.noTmr, vars = "rna_IGKC")
lymphocyte.noTmr$IGKC_zero <- ifelse(ttt$rna_IGKC == 0, "zero", "non-zero")
table(lymphocyte.noTmr$IGKC_zero, lymphocyte.noTmr$cluster)

library(viridisLite)
library(scCustomize)

feature = c("IGKC")
pp <- Plot_Density_Custom(seurat_object = lymphocyte.noTmr,pt.size =1,  features = feature,  reduction = "tsne", combine = F, aspect_ratio = 0.9) +
  NoLegend() 

pp + theme(
  axis.line = element_blank(), axis.text.x = element_blank(), 
  axis.text.y = element_blank(), axis.title = element_blank(), 
  axis.ticks = element_blank(),
  title = element_blank()) 
ggsave("./figures/Plasma_mecha/0416_sc.lymph.DimPlot.IGKC.pdf", width = 7*0.7, height = 6*0.7)

table(lymphocyte.noTmr$orig.ident, lymphocyte.noTmr$cluster)

xxx <- read_delim("data/scRNA_revise//scRNA_cluster.txt")
lymphocyte.noTmr$IME <- ifelse(lymphocyte.noTmr$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "IME"], 
                               "IME", "non-IME")

lymphocyte.noTmr$cluster2 <- lymphocyte.noTmr$cluster
lymphocyte.noTmr$cluster2[!(lymphocyte.noTmr$cluster2 %in% c("plasmaB", "matureB"))] <- "Other"

Idents(lymphocyte.noTmr) <- "IME"
DimPlot(lymphocyte.noTmr, reduction = "tsne", group.by = "IME", label = T, pt.size = 0.1) +
  theme(legend.position = "none") +
  ggtitle("IME vs non-IME") +
  theme(plot.title = element_text(hjust = 0.5))


color_scheme2 <- c(#"CD4T" = '#66c2a5', "CD8T" = '#fc8d62',
                   "matureB"  = '#8da0cb', "plasmaB" ='#e78ac3', "Other" = "grey" 
                   #"NK"= '#a6d854', "NKT" = '#ffd92f', "Treg"  = "#80b1d3"
                   )

order <- c("matureB", "plasmaB", "Other")
order <- rev(order)
DimPlot(object = lymphocyte.noTmr, reduction = "tsne",
        group.by = "cluster2", 
        cols = color_scheme2,
        label = T, 
        pt.size = 1, order = order,
        label.size = 0) & NoAxes()
ggsave("./figures/Plasma_mecha/0416_sc.lymph.DimPlot_cluster2.pdf", width = 7, height = 6)

DimPlot(object = lymphocyte.noTmr, reduction = "tsne",
        group.by = "IME", 
        cols = c("IME" = "red", "non-IME" = "black"),
        label = T, 
        pt.size = 1, order = rev(c("IME", "non-IME")),
        label.size = 0) & NoAxes()
ggsave("./figures/Plasma_mecha/0416_sc.lymph.DimPlot_IME.pdf", width = 7, height = 6)


# pie for lymphocyte ---- 
xxx <- read_delim("data/scRNA_revise//scRNA_cluster.txt")
lymphocyte.noTmr$pc <- ifelse(lymphocyte.noTmr$orig.ident %in% xxx$Cohort_ID[xxx$NMF_Cluster == "IME"], "IME", "others")
table(lymphocyte.noTmr$pc, lymphocyte.noTmr$cluster)

count.data <- data.frame(
  Var1 = c("3CD4T","8CD56bright_NK","7CD56dim_NK", "5CD8T", "1B", "6NKT",  "2PC","4Treg" ),
  Freq = c(274, 41, 23, 987, 293, 36, 84, 46)
)
count.data <- count.data[order(-count.data$Freq),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.25*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

#order <- c("CD4T","CD56bright_NK","CD56dim_NK", "CD8T", "B", "NKT",  "PC","Treg" )
ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=0.5,stat = "identity", color = "white",alpha=0.85) +
  #coord_polar("y", start = 0)+
  #geom_text(aes(x=1.1,y = lab.ypos, label = paste0(Freq," (",round(percent,digits = 2),")") ), color = "black")+
  #scale_fill_manual(name=NULL,values=gg_color_hue(7)) +
  scale_fill_manual(values = c("3CD4T" = '#66c2a5', "5CD8T" = '#fc8d62',
                               "1B"  = '#8da0cb', "2PC" ='#e78ac3', "8CD56bright_NK" = "#ccebc5",
                               "7CD56dim_NK"= '#a6d854', "6NKT" = '#ffd92f', "4Treg"  = "#80b1d3")) + 
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(), 
        legend.title = element_blank(),
        #panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
        legend.position='right',legend.text=element_text(size=14,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm')) 

ggsave("./figures/Plasma_mecha/0422_sc.lymph.IME.bar.pdf", width = 3, height = 4)
count.data <- data.frame(
  Var1 = c("3CD4T","8CD56bright_NK","7CD56dim_NK", "5CD8T", "1B", "6NKT",  "2PC","4Treg"),
  Freq = c(138, 69, 96, 286, 27, 108, 5, 25)
)
count.data <- count.data[order(-count.data$Freq),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.25*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=0.5,stat = "identity", color = "white",alpha=0.85) +
  #coord_polar("y", start = 0)+
  #geom_text(aes(x=1.1,y = lab.ypos, label = paste0(Freq," (",round(percent,digits = 2),")") ), color = "black")+
  #scale_fill_manual(name=NULL,values=gg_color_hue(7)) +
  scale_fill_manual(values = c("3CD4T" = '#66c2a5', "5CD8T" = '#fc8d62',
                               "1B"  = '#8da0cb', "2PC" ='#e78ac3', "8CD56bright_NK" = "#ccebc5",
                               "7CD56dim_NK"= '#a6d854', "6NKT" = '#ffd92f', "4Treg"  = "#80b1d3")) + 
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(), 
        legend.title = element_blank(),
        #panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
        legend.position='right',legend.text=element_text(size=14,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm')) 
ggsave("./figures/Plasma_mecha/0422_sc.lymph.other.bar.pdf", width = 3, height = 4)

