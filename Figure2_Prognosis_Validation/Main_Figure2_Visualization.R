#### Visualization codes for Figure 2
#### Integrative analysis of RNA and protein validated the prognostic value of protein clustering in multiple cohorts
#### Author: Jihong TANG; Jiguang WANG

# Fig 2A - Schematic workflow for transcriptomics classifier development ----
## All the figures were manually created using PowerPoint software.

# Fig 2B - Volcano plot of IME signatures ---- 
library(tidyverse)
## prepare the published mes markers
suva_sc <- readxl::read_xlsx("data/0104_Glioma_CellState_Jan04.xlsx", sheet = "cell2019_GBM")
mes.suva.sc <- suva_sc$MES
mes.suva.sc <- mes.suva.sc[!duplicated(mes.suva.sc)]

suva_spatial <- readxl::read_xlsx("data/0104_Glioma_CellState_Jan04.xlsx", sheet = "Suva_spatial")
mes.suva.spa <- c(suva_spatial$MES, suva_spatial$MES.Hyp, suva_spatial$MES.Ast)
mes.suva.spa <- mes.suva.spa[!duplicated(mes.suva.spa)]

verhaak_V1 <- readxl::read_xlsx("data/0104_Glioma_CellState_Jan04.xlsx", sheet = "TCGA_V1")
verhaak_V2 <- readxl::read_xlsx("data/0104_Glioma_CellState_Jan04.xlsx", sheet = "TCGA_V2")
mes.verhaak <- c(verhaak_V1$Mesenchymal, verhaak_V2$Mesenchymal)
mes.verhaak <- mes.verhaak[!duplicated(mes.verhaak)]


gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}

plotTable_up <- plotTable[which(((plotTable$pn_fc >= 0.5 & plotTable$pn_pvalue < 0.05) | (plotTable$gn_fc >= 0.5 & plotTable$gn_pvalue < 0.05)) & (plotTable$gn_fc >= 0.5 & plotTable$pn_fc >= 0.5)), ]
plotTable_down <- plotTable[which(((plotTable$pn_fc <= -0.5 & plotTable$pn_pvalue < 0.05) | (plotTable$gn_fc <= -0.5 & plotTable$gn_pvalue < 0.05)) & (plotTable$gn_fc <= -0.5 & plotTable$pn_fc <= -0.5)), ]
plotTable_other <- plotTable[which(!(plotTable$pid %in% plotTable_up$pid | plotTable$pid %in% plotTable_down$pid)), ]

mes.other <- c(mes.suva.spa, mes.verhaak)
mes.other <- mes.other[!duplicated(mes.other)]
mes.final <- c(
  mes.suva.sc[mes.suva.sc %in% mes.verhaak],
  mes.suva.sc[mes.suva.sc %in% plotTable_up$pid],
  mes.suva.spa[mes.suva.spa %in% plotTable_up$pid],
  mes.verhaak[mes.verhaak %in% plotTable_up$pid]
)
mes.final <- mes.final[!duplicated(mes.final)]
mes.final <- c(mes.final, "LDHA", "SPP1", "MGST1", "NAMPT")

plotTable_mes_up <- plotTable_up[plotTable_up$pid %in% mes.final, ]
plotTable_mes_other <- plotTable_other[plotTable_other$pid %in% mes.final, ]

plotTable_sig1 <- plotTable_up[((plotTable_up$pn_fc > 1.2 & plotTable_up$gn_fc > 1) | plotTable_up$pid %in% c(
  "TYMP", "ORM1", "ORM2", "ENPP6", "ISG15",
  "S100A6", "SERPINA3", "SERPINA1", "IFITM3"
)) & plotTable_up$pid != "IGHV4-30-4", ]
plotTable_sig2 <- plotTable_down[(plotTable_down$pn_fc < -0.8 & plotTable_down$gn_fc < -0.5), ]

mySub2 <- ggplot() +
  theme_classic()
mySub2 <- mySub2 + geom_point(data = plotTable_other, aes(x = (pn_fc), y = gn_fc), fill = "grey", alpha = 0.5, size = 1.5, shape = 21, stroke = 0.2)
mySub2 <- mySub2 + geom_point(data = plotTable_up, aes(x = (pn_fc), y = gn_fc), fill = gg_color_hue(3)[1], alpha = 0.8, size = 2, shape = 21, stroke = 0.2)
mySub2 <- mySub2 + geom_point(data = plotTable_down, aes(x = (pn_fc), y = gn_fc), fill = gg_color_hue(3)[3], alpha = 0.8, size = 2, shape = 21, stroke = 0.2)
mySub2 <- mySub2 + geom_point(data = plotTable_mes_other, aes(x = (pn_fc), y = gn_fc), alpha = 0.8, size = 1.5, shape = 21, stroke = 0.5, color = "black")

mySub2 <- mySub2 + geom_point(data = plotTable_mes_up, aes(x = (pn_fc), y = gn_fc), alpha = 0.8, size = 2, shape = 21, stroke = 0.6, color = "black")

mySub2

mySub2 <- mySub2 + geom_vline(xintercept = 0.5, linetype = 2, color = "grey50") +
  geom_vline(xintercept = -0.5, linetype = 2, color = "grey50") +
  geom_hline(yintercept = 0.5, linetype = 2, color = "grey50") +
  geom_hline(yintercept = -0.5, linetype = 2, color = "grey50")
mySub2 <- mySub2 + ggrepel::geom_text_repel(data = plotTable_mes_up, aes(x = (pn_fc), y = gn_fc, label = pid), size = 2.5, color = "black", max.overlaps = 30)

mySub2

mySub2 <- mySub2 + xlab("Log2 FC(ProteNMF3 vs. Other)") + ylab("Log2 FC(ProteNMF3 vs. Other)") +
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 3.5), breaks = seq(-2, 4, 0.5)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-2.3, 2.3), breaks = seq(-6, 6, 0.5))
mySub2 <- mySub2 + theme(
  panel.background = element_rect(fill = "transparent", color = "transparent"), plot.margin = unit(c(2, 1, 0.5, 1), "lines"), plot.title = element_text(size = 14, vjust = 0.5, hjust = 0.5, face = "bold.italic"),
  legend.key.width = unit(1.5, "cm"), legend.key.height = unit(0.5, "cm"), legend.position = "top",
  legend.margin = margin(t = 0.1, r = 0.1, b = 0, l = 0.1, unit = "cm"), legend.text = element_text(size = 12), axis.text.y = element_text(size = 10, face = "plain", color = "black"),
  axis.text.x = element_text(size = 10, face = "plain", color = "black"), axis.title.x = element_text(size = 12, face = "plain", color = "black"), axis.title.y = element_text(size = 12, hjust = 0.5, vjust = 2, face = "plain", color = "black")
)
figure_4 <- rbind(ggplotGrob(mySub2), size = "last")
mySub2
ggsave(file = "./figures/Fig2B_scatter_pn_gn.pdf", plot = figure_4, bg = "white", width = 12, height = 9, units = "cm", dpi = 600)

# Fig 2C - Pie charts of classification results ----
library(tidyverse)
count.data <- data.frame(
  Var1 = paste0("ProteoNMF", 1:4),
  Freq = c(26 + 58, 33 + 54, 12 + 20, 22 + 48)
)

count.data <- count.data[order(-count.data$Freq), ]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.25 * Freq) %>%
  mutate(percent = Freq / sum(count.data$Freq))
count.data

order <- c(1:nrow(count.data))
ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size = 1, stat = "identity", color = "white", alpha = 0.85) +
  coord_polar("y", start = 0) +
  # geom_text(aes(x=1.1,y = lab.ypos, label = paste0(Freq," (",round(percent,digits = 2),")") ), color = "black")+
  # scale_fill_manual(name=NULL,values=gg_color_hue(7)) +
  scale_fill_manual(values = c("#33a02c", "#e31a1c", "#ff7f00", "#1f78b4")) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_blank(),
    axis.title = element_blank(),
    legend.title = element_blank(),
    # panel.background=element_rect(fill='transparent',color='black',size=1),plot.margin=unit(c(0.5,1,0.5,1),'lines'),
    legend.key.width = unit(0.6, "cm"), legend.key.height = unit(0.6, "cm"), legend.key.size = unit(5, "lines"), legend.key = element_rect(size = 0.1, color = NA),
    legend.position = "right", legend.text = element_text(size = 14, face = "bold.italic"), legend.margin = margin(t = 0.1, r = 0.1, b = 0, l = 0.1, unit = "cm")
  )
ggsave("figures/Fig2B_pie_cgga_both.pdf", width = 5.98, height = 4.72)

# Fig 2D - survival validation ---- 
custom_theme <- function() {
  theme_classic() %+replace%
    theme(
      text = 	element_text(size = 14, face = 'bold'),
      plot.title=element_text(hjust=0.5)
    )
}

### CGGA 325 693 & combined
nmf_cgga325 <- read_delim("./data/1204_nmf_cgga325.txt", delim = "\t")

OS_cgga325_all <- nmf_cgga325[nmf_cgga325$nmf_subtype != "Mixed", c("OS", "Censor", "nmf_subtype")]
colnames(OS_cgga325_all) <- c("OS", "Censor", "Expre")
OS_cgga325_all$OS <- as.numeric(OS_cgga325_all$OS)
OS_cgga325_all$Censor <- as.numeric(OS_cgga325_all$Censor)

nmf_cgga693 <- read_delim("./data/1204_nmf_cgga693.txt", delim = "\t")

OS_cgga693_all <- nmf_cgga693[nmf_cgga693$nmf_subtype != "Mixed", c("OS", "Censor", "nmf_subtype")]
colnames(OS_cgga693_all) <- c("OS", "Censor", "Expre")
OS_cgga693_all$OS <- as.numeric(OS_cgga693_all$OS)
OS_cgga693_all$Censor <- as.numeric(OS_cgga693_all$Censor)

OS_cgga_all <- rbind(OS_cgga325_all, OS_cgga693_all)

tt1 <- OS_cgga_all[OS_cgga_all$Expre %in% c("ProteoNMF2", "ProteoNMF3"), ]
my.Surv <- with(tt1,Surv(time = OS, event = Censor == 1))
surv.fit = survdiff(my.Surv ~ Expre, data = tt1)
surv.fit$pvalue

tt2 <- OS_cgga_all[OS_cgga_all$Expre %in% c("ProteoNMF1", "ProteoNMF3", "ProteoNMF4"), ]
tt2$compare <- ifelse(tt2$Expre == "ProteoNMF3", "ProteoNMF3", "Other")
my.Surv <- with(tt2,Surv(time = OS, event = Censor == 1))
surv.fit = survdiff(my.Surv ~ compare, data = tt2)
surv.fit$pvalue

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
ggsave(paste0("./figures/Fig2D_OS_cgga_vali_3lines.pdf"), width = 7.5, height = 6, units = 'cm', dpi = 600)

# Fig 2E - Forest plot of Cox HR results ----
nmf_cgga325 <- readxl::read_xlsx("./data/0310_MultiVariate_Prognosis.xlsx", sheet = "Vali_CGGA325")

OS_cgga325_all <- nmf_cgga325[nmf_cgga325$Protein_Cluster != "Mixed", c(
  "OS", "Censor", "Protein_Cluster", "Grade2021", "Radio_status (treated=1;un-treated=0)",
  "Chemo_status (TMZ treated=1;un-treated=0)", "Age", "Gender"
)]

colnames(OS_cgga325_all) <- c("OS", "Censor", "pCluster", "Grade2021", "Radiotherapy", "Chemotherapy", "Age", "Gender")
OS_cgga325_all$OS <- as.numeric(OS_cgga325_all$OS)
OS_cgga325_all$Censor <- as.numeric(OS_cgga325_all$Censor)
OS_cgga325_all$Radiotherapy <- as.numeric(OS_cgga325_all$Radiotherapy) %>% as.factor()
OS_cgga325_all$Chemotherapy <- as.numeric(OS_cgga325_all$Chemotherapy) %>% as.factor()
Pth_color <- c("#33a02c", "#e31a1c", "#ff7f00", "#1f78b4")

nmf_CGGA693 <- readxl::read_xlsx("./data/0310_MultiVariate_Prognosis.xlsx", sheet = "Vali_CGGA693")

OS_CGGA693_all <- nmf_CGGA693[nmf_CGGA693$Protein_Cluster != "Mixed", c(
  "OS", "Censor", "Protein_Cluster", "Grade2021", "Radio_status (treated=1;un-treated=0)",
  "Chemo_status (TMZ treated=1;un-treated=0)", "Age", "Gender"
)]

colnames(OS_CGGA693_all) <- c("OS", "Censor", "pCluster", "Grade2021", "Radiotherapy", "Chemotherapy", "Age", "Gender")
OS_CGGA693_all$OS <- as.numeric(OS_CGGA693_all$OS)
OS_CGGA693_all$Censor <- as.numeric(OS_CGGA693_all$Censor)
OS_CGGA693_all$Radiotherapy <- as.numeric(OS_CGGA693_all$Radiotherapy) %>% as.factor()
OS_CGGA693_all$Chemotherapy <- as.numeric(OS_CGGA693_all$Chemotherapy) %>% as.factor()
Pth_color <- c("#33a02c", "#e31a1c", "#ff7f00", "#1f78b4")

OS_CGGA_all <- rbind(OS_cgga325_all, OS_CGGA693_all)
OS_CGGA_all$ProteinCluster <- ifelse(OS_CGGA_all$pCluster == "IME", "IME",
  ifelse(OS_CGGA_all$pCluster == "PPR", "PPR", "AllOther")
)
OS_CGGA_all$ProteinCluster <- as.factor(OS_CGGA_all$ProteinCluster)

OS_CGGA_all <- as.data.frame(OS_CGGA_all)
res.cox.cgga <- coxph(Surv(OS, Censor) ~ ProteinCluster + Grade2021 + Radiotherapy + Chemotherapy + Age + Gender, data = OS_cgga_all)
res.cox.cgga

source("ggforest_plot_custom.R")
ggforest_JH(res.cox.cgga,
  fontsize = .8,
  noDigits = 2
)
ggsave("figures/Fig2E_coxph_cgga.pdf", width = 7, height = 6)

# Fig 2F - Evolutionary changes of RNA signature scores ----
nmf_IR_1123 <- readxl::read_xlsx("./data/1124_proteomics_subgroups_IR_compare.xlsx")

plt_IR <- nmf_IR_1123[, seq(1, 14)]

ht_cl <- plt_IR[, c("Patient_ID", "nmf_subtype_I", "nmf_subtype_R")]
ht_input <- t(plt_IR[, c(seq(2, 5), seq(8, 11))])
ht_input <- ht_input[, ht_cl$Patient_ID]

df_wide <- ht_input %>% as.data.frame()
rownames(df_wide) <- c("ProMIX_I", "ProPPR_I", "ProMES_I", "ProNEU_I", "ProMIX_R", "ProPPR_R", "ProMES_R", "ProNEU_R")
df_wide_Trans <- 2^df_wide
ww1 <- t(df_wide_Trans[seq(1, 4), ]) %>% apply(., 1, function(x) x / sum(x))
ww2 <- t(df_wide_Trans[seq(5, 8), ]) %>% apply(., 1, function(x) x / sum(x))

df_wide2 <- rbind(ww1, ww2) %>% as.data.frame()
df_long <- df_wide2 %>%
  rownames_to_column("Proteomics") %>%
  gather(Patient_ID, value, -Proteomics)

df_long$group <- sapply(df_long$Proteomics, function(x) substr(x, 1, 6))
df_long$IR <- sapply(df_long$Proteomics, function(x) substr(x, 8, 8))

dd <- df_long %>%
  group_by(Patient_ID, group) %>%
  summarise(diff = diff(value))

df_long <- merge(df_long, dd, by = c("Patient_ID", "group"))
plt_IRcomapre <- merge(df_long, nmf_IR_1123, by = "Patient_ID")

plt_IR_diff <- plt_IRcomapre %>% filter(abs(diff) > 0.2)
plt_IR_nodiff <- plt_IRcomapre %>% filter(abs(diff) <= 0.2)

Fig5B <- ggplot() +
  theme_classic()
Fig5B <- Fig5B + geom_boxplot(data = plt_IR_diff, aes(x = IR, y = value), width = 0.4, size = 0.5, alpha = 0.7, fill = "transparent", outlier.colour = "white")
Fig5B <- Fig5B + geom_line(data = plt_IR_nodiff, aes(x = IR, y = value, group = Patient_ID), size = 0.25, color = "#969696", lty = "dashed")

Fig5B <- Fig5B + geom_point(data = plt_IR_nodiff, aes(x = IR, y = value, color = IR, group = Patient_ID), size = 2, alpha = 0.7, shape = 16) +
  facet_wrap(factor(group, levels = c("ProMIX", "ProPPR", "ProMES", "ProNEU")) ~ ., ncol = 4)
Fig5B <- Fig5B + geom_line(data = plt_IR_diff, aes(x = IR, y = value, group = Patient_ID), size = 0.3, color = "black", lty = "dashed")
Fig5B <- Fig5B + geom_point(data = plt_IR_diff, aes(x = IR, y = value, fill = IR, group = Patient_ID), size = 2, alpha = 0.7, stroke = .25, shape = 21, color = "black") +
  facet_wrap(factor(group, levels = c("ProMIX", "ProPPR", "ProMES", "ProNEU")) ~ ., ncol = 4)

Fig5B <- Fig5B + theme(
  panel.background = element_rect(fill = "transparent", color = "transparent"), plot.margin = unit(c(2, 2, 2, 2), "lines"),
  plot.title = element_text(size = 34, vjust = 0.5, hjust = 0.5, face = "bold.italic", color = "transparent"), text = element_text(size = 14, face = "bold"),
  legend.key.width = unit(0.6, "cm"), legend.key.height = unit(0.6, "cm"), legend.position = "none", legend.text = element_text(size = 14, hjust = 0, face = "bold"),
  axis.text.x = element_text(size = 12, face = "bold", color = "black"), axis.text.y = element_text(size = 12, face = "bold", color = "black"),
  axis.title.x = element_text(size = 14, face = "plain", color = "black"), axis.title.y = element_text(size = 14, face = "plain", color = "black")
)

Fig5B <- Fig5B + scale_y_continuous(expand = c(0, 0), limits = c(0, 1.03), breaks = seq(0, 1, 0.2))
Fig5B <- Fig5B + ylab("Normalized score") + xlab(NULL)
Fig5B <- Fig5B + scale_fill_manual(values = c("I" = "#009BFF", "R" = "#FF5D49"))
Fig5B <- Fig5B + scale_color_manual(values = c("I" = "#009BFF", "R" = "#FF5D49"))
Fig5B
figure_2 <- rbind(ggplotGrob(Fig5B), size = "last")
ggsave(file = "./figures/Fig2F_boxplot_compare.pdf", plot = figure_2, bg = "white", width = 15, height = 9, units = "cm", dpi = 600)

# Fig 2G - Evolutional changes of protein clusters ---- 
library(ggalluvial)
library(tidyverse)
ddd_IR_nmf <- readxl::read_xlsx("data/1124_proteomics_subgroups_IR_compare.xlsx")

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
ggsave(file="./figures/Fig2G_alluvial_plot.pdf", plot=last_plot(),bg = 'white', width = 10, height = 9, units = 'cm', dpi = 600)