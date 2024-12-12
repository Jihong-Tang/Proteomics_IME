#### Visualization codes for Main Figure 4; 2024-12-11
#### Histopathological imaging analysis identified enriched gemistocytic differentiation in the IDHm-IME gliomas
#### Author: Jihong TANG; Jiguang WANG

# Fig 4A; Fig4C; Fig 4H Histopathology images ---- 
## All the histopatholoy images are annotated by licensed pathologist

# Fig 4B Boxplot for fraction comparison ---- 
prop.test(x = c(6, 3), n = c(6, 28))

specie <- c(rep("ProMES", 2), rep("Non-ProMES", 2))
condition <- rep(c("GD", "Non-GD") , 2)
value <- c(6, 0, 3, 25)
data <- data.frame(specie,condition,value)

# Stacked + percent
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="fill", stat="identity", show.legend = F, width = .7) + 
  scale_y_continuous(expand = c(0,0), breaks = seq(0,1,0.2), labels = scales::percent) + 
  scale_x_discrete(limits = c("ProMES", "Non-ProMES")) + 
  scale_fill_manual(values = c("GD" = '#e41a1c',
                               "Non-GD" = '#377eb8')) +
  labs(x = "", y = "Percentage") + 
  theme_classic() + 
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_blank(),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

ggsave("figures/Fig4B_CGPA_GD_IMEcompare.pdf", width = 3, height = 4)

prop.test(x = c(23, 43), n = c(26, 146))

specie <- c(rep("ProMES", 2), rep("Non-ProMES", 2))
condition <- rep(c("GD", "Non-GD") , 2)
value <- c(23, 3, 43, 103)
data <- data.frame(specie,condition,value)

# Stacked + percent
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="fill", stat="identity", show.legend = F, width = .7) + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0,1,0.2), labels = scales::percent) + 
  scale_x_discrete(limits = c("ProMES", "Non-ProMES")) + 
  scale_fill_manual(values = c("GD" = '#e41a1c',
                               "Non-GD" = '#377eb8')) +
  labs(x = "", y = "Percentage") + 
  theme_classic() + 
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_blank(),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

ggsave("figures/Fig4B_TCGA_GD_IMEcompare.pdf", width = 3, height = 4)

# Fig 4D AI model pathologist compare ----
library(tidyverse)
df_compare <- readxl::read_xlsx("../Fig_codes/data/Fig4/GEM_Model_Percent_Correlation.xlsx")
lm <- lm(df_compare$Model_percent_0408 ~ df_compare$Human_percent, data = df_compare)
summary(lm)

ggplot(df_compare, aes(x = Human_percent, y = Model_percent_0408)) + 
  geom_point(size = 2.5, color = "#3c77af") + 
  geom_smooth(method = "lm", se = F, color = "black") + 
  labs(x = "Pathologist (%)", y = "AI (%)") + 
  theme_classic() + 
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
ggsave("figures/Fig4D_scatter_Patho_model_correlation.pdf", width = 3.3, height = 3)

scc.test <- cor.test(df_compare$Model_percent_0408, df_compare$Human_percent, method = "spearman")
scc.test

pcc.test <- cor.test(df_compare$Model_percent_0408, df_compare$Human_percent, method = "pearson")
pcc.test

# Fig 4E Classification AUC plots ---- 
library(tidyverse)
res_image <- read_delim("results/1021_LOOCV_image_classifier.txt")
cutoff <- 0.2

frac_to_prob <- function(frac) {
  diff <- frac - cutoff
  if (diff > 0) {
    return(0.5 + 0.5 * diff / (1 - cutoff))
  } else {
    return(0.5 + 0.5 * diff / (cutoff))
  }
}

res_image$pred <- sapply(res_image$median_ratio, frac_to_prob)
res_image$label <- ifelse(res_image$IDHmes == "Yes", 1, 0)


library(plotROC)
library(pROC)
library(PRROC)
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}

roc.loocv <- roc(res_image$label, res_image$pred)

roc.loocv$auc


plt.roc.loocv <- ggroc(roc.loocv, color = "#ff7f00", linetype = 1, size = 1)
plt.roc.loocv <- plt.roc.loocv + theme_classic() +
  theme(
    panel.background = element_rect(fill = "transparent", color = "black"), plot.margin = unit(c(1, 1, 1, 1), "lines"), plot.title = element_text(size = 24, vjust = 0.5, hjust = 0.5, face = "bold.italic"),
    text = element_text(size = 24, vjust = 1.4, hjust = 0.5, face = "bold"), legend.key.width = unit(1, "cm"), legend.key.height = unit(0.5, "cm"), legend.position = "bottom",
    legend.text = element_text(size = 12, face = "plain"), axis.text.y = element_text(size = 16, vjust = 0.5, hjust = 1, face = "plain", color = "black"), legend.title = element_text(size = 12, vjust = 0.5, hjust = 0, face = "plain"),
    axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust = 1, face = "plain", color = "black"), axis.title.x = element_text(size = 20, vjust = 0, hjust = 0.5, face = "plain", color = "black"),
    axis.title.y = element_text(size = 20, hjust = 0.5, vjust = 2, face = "plain", color = "black"),
    strip.text = element_text(size = 18, face = "bold", vjust = 0.5, hjust = 0.5), strip.background = element_rect(colour = "black", fill = gg_color_hue(3))
  )
plt.roc.loocv
ggsave("figures/Fig4E_Imaging_IME_classifier_ROC.pdf", width = 3, height = 3)

# Fig 4G Evolution fraction comparison ---- 
pair_all <- readxl::read_xlsx("data/0324_PA_Longitunidal_compare_AI_Model.xlsx")

cutoffs <- seq(0, 1, 0.01)

get_GD_ratio_IR <- function(cutoff) {
  Ini_mean <- nrow(pair_all[pair_all$mean_ini2 >= cutoff, ]) / nrow(pair_all)
  Ini_median <- nrow(pair_all[pair_all$median_Ini >= cutoff, ]) / nrow(pair_all)

  Rec_mean <- nrow(pair_all[pair_all$mean_rec2 >= cutoff, ]) / nrow(pair_all)
  Rec_median <- nrow(pair_all[pair_all$median_Rec >= cutoff, ]) / nrow(pair_all)

  out <- data.frame(cutoff, Ini_mean, Ini_median, Rec_mean, Rec_median)
}

res_GD_ratio_IR <- do.call(rbind, lapply(cutoffs, get_GD_ratio_IR))

ggplot(res_GD_ratio_IR) +
  geom_line(aes(x = cutoff, y = Ini_mean), size = 1, color = "#377eb8") +
  geom_line(aes(x = cutoff, y = Rec_mean), size = 1, color = "#e41a1c") +
  scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0.01, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 0.8)) +
  labs(x = "GD percentage cutoff", y = "Patient percentage") +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "transparent", color = "transparent"), plot.margin = unit(c(2, 2, 2, 2), "lines"),
    plot.title = element_text(size = 34, vjust = 0.5, hjust = 0.5, face = "bold.italic", color = "transparent"), text = element_text(size = 14, face = "bold"),
    legend.key.width = unit(0.6, "cm"), legend.key.height = unit(0.6, "cm"), legend.position = "none", legend.text = element_text(size = 14, hjust = 0, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"), axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_text(size = 14, face = "plain", color = "black"), axis.title.y = element_text(size = 14, face = "plain", color = "black")
  )
ggsave("figures/Fig4G_PA_cohort_fraction_updates.pdf", width = 5, height = 4)