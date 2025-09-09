#### Visualization codes for Figure 6
#### Spatial lymphocytes infiltration in the IDHm-IME gliomas
#### Author: Jihong TANG; Jiguang WANG

# IHC ----
library(tidyverse)
library(ggbeeswarm)
## CD3+ ---- 
dd <- readxl::read_xlsx("data/Pathology_IHC/0421_IHC_RawData.xlsx", sheet = "CD3Raw")
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

ggsave(file= paste0("./figures/Pathology_IHC_All/0421_CD3_IHC.pdf"), plot=figure_2,bg = 'white', width =11, height = 5.5, units = 'cm', dpi = 600)


## CD20+ ---- 
dd <- readxl::read_xlsx("data/Pathology_IHC/0421_IHC_RawData.xlsx", sheet = "CD20Raw")
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


ggsave(file= paste0("./figures/Pathology_IHC_All/0421_CD20_IHC.pdf"), plot=figure_2,bg = 'white', width =11, height = 5.5, units = 'cm', dpi = 600)


## CD8+ ---- 
dd <- readxl::read_xlsx("data/Pathology_IHC/0421_IHC_RawData.xlsx", sheet = "CD8Raw")
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

ggsave(file= paste0("./figures/Pathology_IHC_All/0421_CD8_IHC_v2.pdf"), plot=figure_2,bg = 'white', width =10, height = 5.5, units = 'cm', dpi = 600)


## CD4+ ---- 
dd <- readxl::read_xlsx("data/Pathology_IHC/0421_IHC_RawData.xlsx", sheet = "CD4Raw")
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

ggsave(file= paste0("./figures/Pathology_IHC_All/0421_CD4_IHC.pdf"), plot=figure_2,bg = 'white', width =11, height = 5.5, units = 'cm', dpi = 600)


## CD38+ ---- 
library(tidyverse)
library(ggbeeswarm)
dd <- readxl::read_xlsx("data/Pathology_IHC/0421_IHC_RawData.xlsx", sheet = "CD38Raw")
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

ggsave(file= paste0("./figures/Pathology_IHC_All/0421_CD38_IHC.pdf"), plot=figure_2,bg = 'white', width =11, height = 6.5, units = 'cm', dpi = 600)



## CD138+ ---- 
library(tidyverse)
library(ggbeeswarm)
dd <- readxl::read_xlsx("data/Pathology_IHC/0421_IHC_RawData.xlsx", sheet = "CD138Raw")
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

ggsave(file= paste0("./figures/Pathology_IHC_All/0421_CD138_IHC.pdf"), plot=figure_2,bg = 'white', width =11, height = 5.5, units = 'cm', dpi = 600)



## PD1+ ---- 
library(tidyverse)
library(ggbeeswarm)
dd <- readxl::read_xlsx("data/Pathology_IHC/0421_IHC_RawData.xlsx", sheet = "PD1Raw")
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
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2, 115),breaks = seq(0,600,20)) 
NAN_plot <- NAN_plot + scale_x_discrete(limits = c("Mut", "WT"), labels = c("Mut", "WT") )
NAN_plot<- NAN_plot +ylab(paste0("PD1+ Number")) +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3")) + 
  coord_flip()
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")

ggsave(file= paste0("./figures/Pathology_IHC_All/0421_PD1_IHC.pdf"), plot=figure_2,bg = 'white', width =11, height = 5.5, units = 'cm', dpi = 600)
