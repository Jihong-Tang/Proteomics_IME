library(tidyverse)
# 2C - survival validation ---- 

custom_theme <- function() {
  theme_classic() %+replace%
    theme(
      text = 	element_text(size = 14, face = 'bold'),
      plot.title=element_text(hjust=0.5)
    )
}
# Fig 2F -> 2G - Survival survival validation ----
## 2F TCGA ---- 
library(tidyverse)
library(survival)
library(survminer)
nmf_lgg_1123 <- readxl::read_xlsx("../Fig_codes/data/Fig2/nmf_lgg_res_1123.xlsx")
pid_ProMES <- nmf_lgg_1123$Case[nmf_lgg_1123$nmf_subtype == "ProteoNMF3"]
pid_NonProMES <- nmf_lgg_1123$Case[nmf_lgg_1123$nmf_subtype %in% c("ProteoNMF1", "ProteoNMF2", "ProteoNMF4")]
cl.lgg.p <- read_delim("../Fig_codes/data/Fig2/tcga_cbioportal_firehoce/data_clinical_patient.txt", skip = 4)
cl.lgg.p <- cl.lgg.p %>% filter(PATIENT_ID %in% pid_ProMES | PATIENT_ID %in% pid_NonProMES)

pid_NMF1 <- pid_NonProMES[which(pid_NonProMES %in% nmf_lgg_1123$Case[nmf_lgg_1123$nmf_subtype == "ProteoNMF1"] )]
pid_NMF2 <- pid_NonProMES[which(pid_NonProMES %in% nmf_lgg_1123$Case[nmf_lgg_1123$nmf_subtype == "ProteoNMF2"] )]
pid_NMF4 <- pid_NonProMES[which(pid_NonProMES %in% nmf_lgg_1123$Case[nmf_lgg_1123$nmf_subtype == "ProteoNMF4"] )]

cl.OS <- cl.lgg.p[, c("PATIENT_ID", "OS_MONTHS", "OS_STATUS")]

cl.OS.compare <- cl.OS
cl.OS.compare$Censor <- ifelse(cl.OS.compare$OS_STATUS == "1:DECEASED", 1, 0)
cl.OS.compare$Expre <- ifelse(cl.OS.compare$PATIENT_ID %in% pid_ProMES, "NMF3", 
                              ifelse(cl.OS.compare$PATIENT_ID %in% pid_NMF2, "NMF2", 
                                     ifelse(cl.OS.compare$PATIENT_ID %in% pid_NMF1, "NMF1", 
                                            ifelse(cl.OS.compare$PATIENT_ID %in% pid_NMF4, "NMF4", "Other"))))
cl.OS.compare <- cl.OS.compare[, c("OS_MONTHS", "Censor", "Expre")]
colnames(cl.OS.compare) <- c("OS", "Censor", "Expre")
cl.OS.compare$OS <- as.numeric(cl.OS.compare$OS)
cl.OS.compare$Censor <- as.numeric(cl.OS.compare$Censor)
Pth_color <- c('#33a02c','#e31a1c','#ff7f00','#1f78b4')

p <- ggsurvplot(
  fit = survfit(Surv(OS, Censor) ~ Expre, data = cl.OS.compare %>% filter(Expre != "Other")), size = .5,
  risk.table = F,pval = F,
  ggtheme = custom_theme(), 
  xlab = "Years from diagnosis", ylab = "Survival Probability",
  
  legend = "none",legend.title = "",font.legend = c(10, "bold", "black"), 
  
  break.x.by = 24, xscale = "m_y",break.y.by = 0.2,axes.offset = T,
  font.x = c(14, "plain", "black"),font.y = c(14, "plain", "black"),font.tickslab = c(12, "bold", "black"), 
  palette = Pth_color
)
p

ggsave(paste0("./figures/Fig2/1204_OS_tcga_vali.pdf"), width = 7.5, height = 6, units = 'cm', dpi = 600)

my.Surv <- with(cl.OS.compare,Surv(time = OS, event = Censor == 1))
surv.fit = survdiff(my.Surv ~ Expre, data = cl.OS.compare)
surv.fit$pvalue

tt1 <- cl.OS.compare[cl.OS.compare$Expre %in% c("NMF2", "NMF3"), ]
my.Surv <- with(tt1,Surv(time = OS, event = Censor == 1))
surv.fit = survdiff(my.Surv ~ Expre, data = tt1)
surv.fit$pvalue

tt2 <- cl.OS.compare[cl.OS.compare$Expre %in% c("NMF1", "NMF3", "NMF4"), ]
tt2$compare <- ifelse(tt2$Expre == "NMF3", "NMF3", "Other")
my.Surv <- with(tt2,Surv(time = OS, event = Censor == 1))
surv.fit = survdiff(my.Surv ~ compare, data = tt2)
surv.fit$pvalue

cl.OS.compare$Expre2 <- cl.OS.compare$Expre
cl.OS.compare$Expre2[cl.OS.compare$Expre %in% c("NMF1", "NMF4")] <- "others"
Pth_color <- c('#e31a1c','#ff7f00', "#298c70")

p <- ggsurvplot(
  fit = survfit(Surv(OS, Censor) ~ Expre2, data = cl.OS.compare %>% filter(Expre != "Other")), size = .5,
  risk.table = F,pval = F,
  ggtheme = custom_theme(), 
  xlab = "Years from diagnosis", ylab = "Survival Probability",
  
  legend = "none",legend.title = "",font.legend = c(10, "bold", "black"), 
  
  break.x.by = 24, xscale = "m_y",break.y.by = 0.2,axes.offset = T,
  font.x = c(14, "plain", "black"),font.y = c(14, "plain", "black"),font.tickslab = c(12, "bold", "black"), 
  palette = Pth_color
)
p
ggsave(paste0("./figures/Fig2/1130_OS_tcga_3line.pdf"), width = 7.5, height = 6, units = 'cm', dpi = 600)


## 2F CGGA 325 693 & combined ---- 

### cgga325 ---- 
nmf_cgga325 <- read_delim("../Fig_codes/data/Fig2/1204_nmf_cgga325.txt", delim = "\t")

OS_cgga325_all <- nmf_cgga325[nmf_cgga325$nmf_subtype != "Mixed", c("OS", "Censor", "nmf_subtype")]
colnames(OS_cgga325_all) <- c("OS", "Censor", "Expre")
OS_cgga325_all$OS <- as.numeric(OS_cgga325_all$OS)
OS_cgga325_all$Censor <- as.numeric(OS_cgga325_all$Censor)
Pth_color <- c('#33a02c','#e31a1c','#ff7f00','#1f78b4')

my.Surv <- with(OS_cgga325_all,Surv(time = OS, event = Censor == 1))
surv.fit = survdiff(my.Surv ~ Expre, data = OS_cgga325_all)
surv.fit$pvalue

tt1 <- OS_cgga325_all[OS_cgga325_all$Expre %in% c("ProteoNMF2", "ProteoNMF3"), ]
my.Surv <- with(tt1,Surv(time = OS, event = Censor == 1))
surv.fit = survdiff(my.Surv ~ Expre, data = tt1)
surv.fit$pvalue

tt2 <- OS_cgga325_all[OS_cgga325_all$Expre %in% c("ProteoNMF1", "ProteoNMF3", "ProteoNMF4"), ]
tt2$compare <- ifelse(tt2$Expre == "ProteoNMF3", "ProteoNMF3", "Other")
my.Surv <- with(tt2,Surv(time = OS, event = Censor == 1))
surv.fit = survdiff(my.Surv ~ compare, data = tt2)
surv.fit$pvalue

p <- ggsurvplot(
  fit = survfit(Surv(OS, Censor) ~ Expre, data = OS_cgga325_all), size = .5,
  risk.table = F,pval = F,
  ggtheme = custom_theme(), 
  xlab = "Years from diagnosis", ylab = "Survival Probability",
  
  legend = "none",legend.title = "",font.legend = c(10, "bold", "black"), 
  
  break.x.by = 365.25*2, xscale = "d_y",break.y.by = 0.2,axes.offset = T,
  font.x = c(14, "plain", "black"),font.y = c(14, "plain", "black"),font.tickslab = c(12, "bold", "black"), 
  palette = Pth_color
)
p
ggsave(paste0("./figures/Fig2/1204_OS_cgga325_vali.pdf"), width = 7.5, height = 6, units = 'cm', dpi = 600)

### cgga693 ----
nmf_cgga693 <- read_delim("../Fig_codes/data/Fig2/1204_nmf_cgga693.txt", delim = "\t")

OS_cgga693_all <- nmf_cgga693[nmf_cgga693$nmf_subtype != "Mixed", c("OS", "Censor", "nmf_subtype")]
colnames(OS_cgga693_all) <- c("OS", "Censor", "Expre")
OS_cgga693_all$OS <- as.numeric(OS_cgga693_all$OS)
OS_cgga693_all$Censor <- as.numeric(OS_cgga693_all$Censor)
Pth_color <- c('#33a02c','#e31a1c','#ff7f00','#1f78b4')

my.Surv <- with(OS_cgga693_all,Surv(time = OS, event = Censor == 1))
surv.fit = survdiff(my.Surv ~ Expre, data = OS_cgga693_all)
surv.fit$pvalue

tt1 <- OS_cgga693_all[OS_cgga693_all$Expre %in% c("ProteoNMF2", "ProteoNMF3"), ]
my.Surv <- with(tt1,Surv(time = OS, event = Censor == 1))
surv.fit = survdiff(my.Surv ~ Expre, data = tt1)
surv.fit$pvalue

tt2 <- OS_cgga693_all[OS_cgga693_all$Expre %in% c("ProteoNMF1", "ProteoNMF3", "ProteoNMF4"), ]
tt2$compare <- ifelse(tt2$Expre == "ProteoNMF3", "ProteoNMF3", "Other")
my.Surv <- with(tt2,Surv(time = OS, event = Censor == 1))
surv.fit = survdiff(my.Surv ~ compare, data = tt2)
surv.fit$pvalue

p <- ggsurvplot(
  fit = survfit(Surv(OS, Censor) ~ Expre, data = OS_cgga693_all), size = .5,
  risk.table = F,pval = F,
  ggtheme = custom_theme(), 
  xlab = "Years from diagnosis", ylab = "Survival Probability",
  
  legend = "none",legend.title = "",font.legend = c(10, "bold", "black"), 
  
  break.x.by = 365.25*2, xscale = "d_y",break.y.by = 0.2,axes.offset = T,
  font.x = c(14, "plain", "black"),font.y = c(14, "plain", "black"),font.tickslab = c(12, "bold", "black"), 
  palette = Pth_color
)
p
ggsave(paste0("./figures/Fig2/1204_OS_cgga693_vali.pdf"), width = 7.5, height = 6, units = 'cm', dpi = 600)

### cgga combined ----
OS_cgga_all <- rbind(OS_cgga325_all, OS_cgga693_all)

#OS_cgga_all <- nmf_cgga[nmf_cgga$nmf_subtype != "Mixed", c("OS", "Censor", "nmf_subtype")]
colnames(OS_cgga_all) <- c("OS", "Censor", "Expre")
OS_cgga_all$OS <- as.numeric(OS_cgga_all$OS)
OS_cgga_all$Censor <- as.numeric(OS_cgga_all$Censor)
Pth_color <- c('#33a02c','#e31a1c','#ff7f00','#1f78b4')

my.Surv <- with(OS_cgga_all,Surv(time = OS, event = Censor == 1))
surv.fit = survdiff(my.Surv ~ Expre, data = OS_cgga_all)
surv.fit$pvalue

tt1 <- OS_cgga_all[OS_cgga_all$Expre %in% c("ProteoNMF2", "ProteoNMF3"), ]
my.Surv <- with(tt1,Surv(time = OS, event = Censor == 1))
surv.fit = survdiff(my.Surv ~ Expre, data = tt1)
surv.fit$pvalue

tt2 <- OS_cgga_all[OS_cgga_all$Expre %in% c("ProteoNMF1", "ProteoNMF3", "ProteoNMF4"), ]
tt2$compare <- ifelse(tt2$Expre == "ProteoNMF3", "ProteoNMF3", "Other")
my.Surv <- with(tt2,Surv(time = OS, event = Censor == 1))
surv.fit = survdiff(my.Surv ~ compare, data = tt2)
surv.fit$pvalue


p <- ggsurvplot(
  fit = survfit(Surv(OS, Censor) ~ Expre, data = OS_cgga_all), size = .5,
  risk.table = F,pval = F,
  ggtheme = custom_theme(), 
  xlab = "Years from diagnosis", ylab = "Survival Probability",
  
  legend = "none",legend.title = "",font.legend = c(10, "bold", "black"), 
  
  break.x.by = 365.25*2, xscale = "d_y",break.y.by = 0.2,axes.offset = T,
  font.x = c(14, "plain", "black"),font.y = c(14, "plain", "black"),font.tickslab = c(12, "bold", "black"), 
  palette = Pth_color
)
p
ggsave(paste0("./figures/Fig2/1128_OS_cgga_vali.pdf"), width = 7.5, height = 6, units = 'cm', dpi = 600)


OS_cgga_all$Expre2 <- OS_cgga_all$Expre
OS_cgga_all$Expre2[OS_cgga_all$Expre %in% c("ProteoNMF1", "ProteoNMF4")] <- "others"
Pth_color <- c("#298c70", '#e31a1c','#ff7f00')

p <- ggsurvplot(
  fit = survfit(Surv(OS, Censor) ~ Expre2, data = OS_cgga_all), size = .5,
  risk.table = F,pval = F,
  ggtheme = custom_theme(), 
  xlab = "Years from diagnosis", ylab = "Survival Probability",
  
  legend = "none",legend.title = "",font.legend = c(10, "bold", "black"), 
  
  break.x.by = 365.25*2, xscale = "d_y",break.y.by = 0.2,axes.offset = T,
  font.x = c(14, "plain", "black"),font.y = c(14, "plain", "black"),font.tickslab = c(12, "bold", "black"), 
  palette = Pth_color
)
p
ggsave(paste0("./figures/Fig2/1130_OS_cgga_vali_3lines.pdf"), width = 7.5, height = 6, units = 'cm', dpi = 600)



count.data <- data.frame(
  Var1 = paste0("ProteoNMF",1:4),
  Freq = c(26+58, 33+54, 12+20, 22+48)
)

count.data <- count.data[order(-count.data$Freq),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.25*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

order <- c(1:nrow(count.data))
ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=1,stat = "identity", color = "white",alpha=0.85) +
  coord_polar("y", start = 0)+
  #geom_text(aes(x=1.1,y = lab.ypos, label = paste0(Freq," (",round(percent,digits = 2),")") ), color = "black")+
  #scale_fill_manual(name=NULL,values=gg_color_hue(7)) +
  scale_fill_manual(values = c('#33a02c','#e31a1c','#ff7f00','#1f78b4')) + 
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

ggsave("figures/Fig2/1128_pie_cgga_both.pdf", width = 5.98, height = 4.72)



# 2D - Evolutional changes of RNA signature scores ---- 
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

nmf_IR_1123 <- readxl::read_xlsx("../Fig_codes/data/Fig5/1124_proteomics_subgroups_IR_compare.xlsx")

plt_IR <- nmf_IR_1123[, seq(1, 14)]

ht_cl <- plt_IR[, c("Patient_ID", "nmf_subtype_I", "nmf_subtype_R")]
ht_input <- t(plt_IR[, c(seq(2, 5), seq(9,12))])
colnames(ht_input) <- ht_cl$Patient_ID

col_order <- ht_cl$Patient_ID[order(ht_cl$nmf_subtype_R)]
row_split <- data.frame(rname = rownames(ht_input), c(rep("I",4), rep("R",4)))
colnames(row_split) <- c("rname", "split")
ht_input_scale <- apply(ht_input, 1, scale) %>% t()
colnames(ht_input_scale) <- colnames(ht_input)

library(ComplexHeatmap)
library(dendextend)
library(circlize)

col_fun = colorRamp2(seq(-1.5, 1.5, 0.1), colorRampPalette(c('#29419E','#698DC9',"grey95",'#F67172','#DC2B18'))(31))
ht_input <- ht_input[, ht_cl$Patient_ID]
df_wide <- ht_input %>% as.data.frame()
rownames(df_wide) <- c("ProMIX_I", "ProPPR_I", "ProMES_I", "ProNEU_I", "ProMIX_R", "ProPPR_R", "ProMES_R", "ProNEU_R")
df_wide_Trans <- df_wide
ww1 <- (df_wide_Trans[seq(1,4),]) #%>% apply(., 1,function(x)x/sum(x))
ww2 <- (df_wide_Trans[seq(5,8),]) #%>% apply(., 1,function(x)x/sum(x))

df_wide2 <- rbind(ww1,ww2) %>% as.data.frame()
df_long <- df_wide2 %>% rownames_to_column("Proteomics") %>% gather(Patient_ID, value, -Proteomics)

df_long$group <- sapply(df_long$Proteomics, function(x)substr(x, 1, 6))
df_long$IR <- sapply(df_long$Proteomics, function(x)substr(x, 8, 8))

dd <- df_long %>% group_by(Patient_ID, group) %>% summarise(diff = diff(value))

df_long <- merge(df_long, dd, by = c("Patient_ID", "group"))
plt_IRcomapre <- merge(df_long, nmf_IR_1123, by = "Patient_ID")

plt_IR_diff <- plt_IRcomapre #%>% filter(group == "ProMES") #%>% filter(abs(diff) > 0.2)
#plt_IR_nodiff <- plt_IRcomapre %>% filter(abs(diff) <= 0.2)

Fig5B <- ggplot() + theme_classic()
Fig5B <- Fig5B + geom_boxplot(data = plt_IR_diff, aes(x = IR, y = value),width = 0.4,size=0.5, alpha = 0.7,fill="transparent", outlier.colour = "white")
#Fig5B <- Fig5B + geom_line(data = plt_IR_nodiff, aes(x = IR, y = value,group=Patient_ID),size=0.25, color = "#dddddd", lty = "dashed")

#Fig5B <- Fig5B + geom_point(data = plt_IR_nodiff, aes(x = IR, y = value, color=nmf_subtype_I,group=Patient_ID),size=2, alpha = 0.7,shape=16) + 
#facet_wrap(factor(group, levels = c("ProMIX", "ProPPR", "ProMES", "ProNEU")) ~., ncol = 4)
Fig5B <- Fig5B + geom_line(data = plt_IR_diff, aes(x = IR, y = value,group=Patient_ID),size=0.25, color = "black",  lty = "dashed")
Fig5B <- Fig5B + geom_point(data = plt_IR_diff, aes(x = IR, y = value, fill=IR,group=Patient_ID),size=2, alpha = 0.7,stroke = .25,shape=21, color = "black") + 
  facet_grid(factor(group, levels = c("ProMIX", "ProPPR", "ProMES", "ProNEU")) ~ . ) 

Fig5B <- Fig5B + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                       plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                       legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                       axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                       axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

Fig5B <- Fig5B + scale_y_continuous(expand=c(0,0),limits=c(-3.4,3.4),breaks = seq(-4,4,1)) +
  scale_x_discrete(limits = c("R","I"))
Fig5B<- Fig5B +ylab("NES") +xlab(NULL) + coord_flip()
Fig5B <- Fig5B +   scale_fill_manual(values = c("I" = "#009BFF", "R" = "#FF5D49")) 
#Fig5B <- Fig5B + scale_color_manual(name=NULL,values =  c("ProteoNMF1" = "#33a02c", "ProteoNMF2" = '#e31a1c', "ProteoNMF3" = '#ff7f00', "ProteoNMF4" = '#1f78b4', "Mixed" = "grey"))
Fig5B
figure_2<-rbind(ggplotGrob(Fig5B),size="last")
ggsave(file="./figures/Fig2/1023_boxplot_compare_ProSig_scores.pdf", plot=figure_2,bg = 'white', width =12, height = 15, units = 'cm', dpi = 600)
ggsave(file="./figures/Fig2/1115_boxplot_compare_ProSig_scores.pdf", plot=figure_2,bg = 'white', width =9, height = 15, units = 'cm', dpi = 600)

library(ggpubr)
plt_IR_diff %>% 
  ggplot()+
  geom_boxplot(aes(x = IR, y = value),width = 0.5)+ 
  geom_line(aes(x = IR, y = value,group = Patient_ID), size=1, color='gray', alpha=0.6)+ 
  geom_point(aes(x = IR, y = value,color=nmf_subtype_I,group=Patient_ID),size=2,shape=16)+
  stat_compare_means(aes(x=IR, y = value), paired = T)+
  scale_fill_manual(values = c("#1f78b4", "#e31a1c"))+
  facet_grid(. ~ group)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "", y = "ProMES score", fill = "")

# 1115 2D alluvial plots ---- 
library(ggalluvial)
ddd_IR_nmf <- readxl::read_xlsx("data/Fig2/1115_proteomics_subgroups_IR_compare.xlsx")

library(ggalluvial)
library(tidyverse)

data<-ddd_IR_nmf %>% filter(nmf_subtype_I != "Mixed") %>% filter(nmf_subtype_R != "Mixed")
data$group2<-"AFM"
data$group2[which(data$nmf_subtype_R=="ProteoNMF2")]<-"PPR"
data$group2[which(data$nmf_subtype_R=="ProteoNMF3")]<-"IME"
data$group2[which(data$nmf_subtype_R=="ProteoNMF4")]<-"NEU"

plot_table<-as.data.frame(table(data$nmf_subtype_I,data$nmf_subtype_R,data$group2))
colnames(plot_table)<-c("Subtype.I","Subtype.R","group","Count")
ggplot(plot_table,
       aes(axis1 = Subtype.I,
           axis2 = Subtype.R,
           y = Count)) +
  geom_alluvium(aes(fill = group)) +
  #scale_fill_manual(values = c(AFM = "#b2df8a", PPR = "#fb9a99", IME = "#fdbf6f", NEU = "#a6cee3"))+
  scale_fill_manual(values = c(AFM = "#33a02c", PPR = "#e31a1c", IME = "#ff7f00", NEU = "#1f78b4"))+

  geom_stratum() +
  #geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
  scale_x_discrete(limits = c("Primary", "Recurrent"),
                   expand = c(.1, .1)) +
  theme_classic()+
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
ggsave(file="./figures/Fig2/1115_2D_alluvial_plot.pdf", plot=last_plot(),bg = 'white', width = 10, height = 9, units = 'cm', dpi = 600)

chiq<-chisq.test(table(data$nmf_subtype_I,data$nmf_subtype_R))
chiq

sum((data$nmf_subtype_I == "ProteoNMF1"))
sum((data$nmf_subtype_R == "ProteoNMF1"))

prop.test(x = c(19, 7), n = c(49, 49), alternative = "two.sided", correct = FALSE)

sum((data$nmf_subtype_I == "ProteoNMF2"))
sum((data$nmf_subtype_R == "ProteoNMF2"))
prop.test(x = c(13, 21), n = c(49, 49), alternative = "two.sided", correct = FALSE)

sum((ddd_IR_nmf$nmf_subtype_I_V2 == "ProteoNMF3"))
sum((ddd_IR_nmf$nmf_subtype_R_V2 == "ProteoNMF3"))
prop.test(x = c(4, 9), n = c(49, 49), alternative = "two.sided", correct = FALSE)

sum((data$nmf_subtype_I == "ProteoNMF4"))
sum((data$nmf_subtype_R == "ProteoNMF4"))
prop.test(x = c(13, 12), n = c(49, 49), alternative = "two.sided", correct = FALSE)

chi_square_result <- chisq.test(table(data$nmf_subtype_I == "ProteoNMF3",data$nmf_subtype_R == "ProteoNMF3"))
residuals(chi_square_result, type = "standardized")

chis <- chisq.test(data$nmf_subtype_I, data$nmf_subtype_R)
