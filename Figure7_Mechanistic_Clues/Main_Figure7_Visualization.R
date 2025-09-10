#### Visualization codes for Figure 7
#### Molecular mechanisms associated with poor prognosis in the IDHm-IME gliomas
#### Author: Jihong TANG; Jiguang WANG

# Fig 7A HeatMap dotplots for visualize important genes DE results ---- 
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

de_pns <- readxl::read_xlsx("data/Fig6_mechanism_candidates.xlsx")
de_pns <- de_pns[de_pns$Cohort == "CGPA", seq(1,5)]
colnames(de_pns) <- c("Gene", "pvalue", "fc", "qvalue", "Cohort")

de_gns <- readxl::read_xlsx("data/Fig6_mechanism_candidates.xlsx", sheet = "RNA")
de_gns_cgga <- de_gns[de_gns$Cohort == "CGGA",seq(1,5)]
colnames(de_gns_cgga) <- c("Gene", "pvalue", "fc", "qvalue", "Cohort")
de_gns_tcga <- de_gns[de_gns$Cohort == "TCGA", seq(1,5)]
colnames(de_gns_tcga) <- c("Gene", "pvalue", "fc", "qvalue", "Cohort")

gns <- c("IGHG1", "IGHG2", "IGHG3", "IGKC", "FCGR2A", "CCL2", "CCR2", "GBP1", "GBP2", "SERPINA3")

plt_pns <- de_pns[de_pns$Gene %in% gns,]
plt_gns_cgga <- de_gns_cgga[de_gns_cgga$Gene %in% gns,]
plt_gns_tcga <- de_gns_tcga[de_gns_tcga$Gene %in% gns,]
plt_all <- rbind(plt_pns, plt_gns_cgga, plt_gns_tcga)

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
        row_order = c("CGPA", "CGGA", "TCGA"),
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
pdf("figures/Fig7A_genes_heatmap.pdf", width = 6, height = 4)
draw(ht)
dev.off()

# Fig 7B scRNA cell clusters barplots ----
count.data <- data.frame(
  Var1 = c("3CD4T","8CD56bright_NK","7CD56dim_NK", "5CD8T", "1B", "6NKT",  "2PC","4Treg" ),
  Freq = c(274, 41, 23, 987, 293, 36, 84, 46)
)
count.data <- count.data[order(-count.data$Freq),]
count.data <- count.data %>%
  mutate(lab.ypos = cumsum(Freq) - 0.25*Freq) %>%
  mutate(percent = Freq/sum(count.data$Freq)) 
count.data

ggplot(count.data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, size=0.5,stat = "identity", color = "white",alpha=0.85) +
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
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
        legend.position='right',legend.text=element_text(size=14,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm')) 
ggsave("./figures/Fig7B_sc.lymph.IME.bar.pdf", width = 3, height = 4)

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
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.key.size = unit(5, 'lines'),legend.key = element_rect(size = 0.1, color = NA),
        legend.position='right',legend.text=element_text(size=14,face='bold.italic'),legend.margin=margin(t=0.1,r=0.1,b=0,l=0.1,unit='cm')) 
ggsave("./figures/Fig7B_sc.lymph.other.bar.pdf", width = 3, height = 4)

# Fig7D Ripley's Cross L-function ----
library(tidyverse)
env_L_nsim99 <- read_delim("data/0710_CrossL_plasma_IDH1_CD38_nsim99.txt")
plot <- env_L_nsim99
plot <- plot[plot$r <= 700, ]
ggplot(plot) +
  geom_vline(xintercept = 31.44531, linetype = "dashed", color = "#dddddd") +
  geom_line(aes(x = r, y = obs), color = "red", linewidth = 0.6) +
  geom_line(aes(x = r, y = theo), color = "black", linewidth = 0.6, linetype = "dashed") +
  geom_ribbon(aes(x = r, ymin = lo, ymax = hi), alpha = 0.1) +
  scale_x_continuous(limits = c(-20, 800), expand = c(0, 0), breaks = seq(0, 800, 200)) +
  scale_y_continuous(limits = c(-20, 800), expand = c(0, 0), breaks = seq(0, 800, 200)) +
  labs(
    x = "Distance (µm)",
    y = "L(r)"
  ) +
  theme_classic() +
  theme( 
    legend.position = "none",
    legend.key.width = unit(.5, "cm"), legend.key.height = unit(0.5, "cm"),
    legend.text = element_text(size = 12, face = "plain"),
    axis.text.y = element_text(size = 16, vjust = 0.5, hjust = 1, face = "plain", color = "black"),
    legend.title = element_text(size = 12, vjust = 0.5, hjust = 0, face = "plain"),
    axis.text.x = element_text(size = 16, angle = 0, face = "plain", color = "black"),
    axis.title.x = element_text(size = 20, vjust = 0, hjust = 0.5, face = "plain", color = "black"),
    axis.title.y = element_text(size = 20, hjust = 0.5, vjust = 2, face = "plain", color = "black")
  )
ggsave("figures/Fig7D_Lcross_Tumor_Plasma.pdf", width = 3.2, height = 3)