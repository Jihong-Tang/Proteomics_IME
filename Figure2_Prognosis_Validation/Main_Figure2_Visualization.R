#### Visualization codes for Main Figure 2; 2024-12-09
#### Integrative analysis of RNA and protein validated the prognostic value of protein clustering in multiple cohorts
#### Author: Jihong TANG; Jiguang WANG

# Fig 2A - Schematic workflow for transcriptomics classifier development ----
## All the figures were manually created using PowerPoint software.

# Fig 2B - Pie charts of classification results ----
count.data <- data.frame(
  Var1 = paste0("ProteoNMF", 1:4),
  Freq = c(85, 59, 32, 58)
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

ggsave("figures/Fig2B_pie_tcga.pdf", width = 5.98, height = 4.72)

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

library(tidyverse)

# Fig 2C - survival validation ---- 
custom_theme <- function() {
  theme_classic() %+replace%
    theme(
      text = 	element_text(size = 14, face = 'bold'),
      plot.title=element_text(hjust=0.5)
    )
}

### TCGA ---- 
library(tidyverse)
library(survival)
library(survminer)
nmf_lgg_1123 <- readxl::read_xlsx("data/nmf_lgg_res_1123.xlsx")
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
ggsave(paste0("./figures/Fig2C_OS_tcga_3line.pdf"), width = 7.5, height = 6, units = 'cm', dpi = 600)


### CGGA 325 693 & combined ---- 
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
ggsave(paste0("./figures/Fig2C_OS_cgga_vali_3lines.pdf"), width = 7.5, height = 6, units = 'cm', dpi = 600)

# Fig 2D - Evolutional changes of RNA signature scores ---- 
library(ggalluvial)
library(tidyverse)
ddd_IR_nmf <- readxl::read_xlsx("data/1115_proteomics_subgroups_IR_compare.xlsx")

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
ggsave(file="./figures/Fig2D_alluvial_plot.pdf", plot=last_plot(),bg = 'white', width = 10, height = 9, units = 'cm', dpi = 600)