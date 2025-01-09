#### Visualization codes for Main Figure 5; 2024-12-12
#### Spatial dynamics of tumor cell states and lymphocyte infiltration in the IDHm-IME gliomas
#### Author: Jihong TANG; Jiguang WANG

# Fig 5A schematic workflow ---- 
## The schematic workflow for the spatial multi-omics experiment design
## are created using Biorender (https://www.biorender.com/)

# Fig 5B data illustration: H&E; Visium HD; CODEX ----
## Left panel: H&E images of the Visium HD capture region
## right panel: multi-staining images of the matched Visium HD capture region
## middle panel: unbiased clustering of the Visium data 
library(tidyverse)
library(Seurat)

load("results/0724_ST_HD_BANKSY02_IniRec_016um.Rdata")

Idents(obj_rec_16_BSK02) <- "BANKSY_snn_res.0.2"

pp <- SpatialDimPlot(obj_rec_16_BSK02,
  image.alpha = 0, cols = c(
    "0" = "#33a02c", "2" = "#ff7f00", "6" = "#b15928", "4" = "#1f78b4", "5" = "#ffff99",
    "1" = "#e31a1c", "3" = "#cab2d6"
  ),
  label = FALSE, crop = T, pt.size.factor = 6.5, shape = 22
)
pp
ggsave("figures/Fig5B_rec_16um_banksy02_res.0.2.c2.pdf", plot = pp, width = 5, height = 4)

obj_ini_16_BSK02$sn_cluster <- NA
obj_ini_16_BSK02$sn_cluster[obj_ini_16_BSK02$BANKSY_snn_res.0.2 %in% c(0, 6)] <- "SN2"
obj_ini_16_BSK02$sn_cluster[obj_ini_16_BSK02$BANKSY_snn_res.0.2 %in% c(1, 7)] <- "LowQ"
obj_ini_16_BSK02$sn_cluster[obj_ini_16_BSK02$BANKSY_snn_res.0.2 == 2] <- "SN5"
obj_ini_16_BSK02$sn_cluster[obj_ini_16_BSK02$BANKSY_snn_res.0.2 %in% c(3, 9)] <- "SN1"
obj_ini_16_BSK02$sn_cluster[obj_ini_16_BSK02$BANKSY_snn_res.0.2 %in% c(4, 10)] <- "SN8"
obj_ini_16_BSK02$sn_cluster[obj_ini_16_BSK02$BANKSY_snn_res.0.2 == 5] <- "SN6"
obj_ini_16_BSK02$sn_cluster[obj_ini_16_BSK02$BANKSY_snn_res.0.2 == 8] <- "SN9"

DimPlot(obj_ini_16_BSK02, group.by = "sn_cluster", label = T)

Idents(obj_ini_16_BSK02) <- "sn_cluster"
p <- SpatialDimPlot(obj_ini_16_BSK02,
  image.alpha = 0, cols = c(
    "SN1" = "#b2df8a", "SN8" = "#fccde5", "SN9" = "#8dd3c7", "LowQ" = "#dddddd", "SN6" = "#ffff99",
    "SN2" = "#fb9a99", "SN5" = "#cab2d6"
  ),
  label = FALSE, crop = T, pt.size.factor = 6.5, shape = 22
)

p
ggsave("figures/Fig5B_ini_16um_banksy02_res.0.2.pdf", plot = p, width = 5, height = 4)

# Fig 5C cell neighbourhood analysis ----

plt_ht_nei <- read_delim("data/Fig5/hd_neighborhood/rec_10_table.txt")
plt_ht_nei <- read_delim("data/Fig5/hd_spatial/")

pltt <- apply(plt_ht_nei, 2, scale)
rownames(pltt) <- paste0("Niche", seq(1:nrow(pltt)))

# #plttt <- log2(plt_ht_nei + 1)
# plttt <- apply(plt_ht_nei, 2, function(x) x-median(x))
# rownames(plttt) <- paste0("Niche", seq(1:nrow(plttt)))

library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(seq(2, -2, -0.1), colorRampPalette(c('#d01c8b','#f1b6da','#f7f7f7','#b8e186','#4dac26'))(41))
col_fun = colorRamp2(seq(-2, 2, 0.1), colorRampPalette(c('#29419E','#698DC9',"grey95",'#F67172','#DC2B18'))(41))
row_order <- c('MES-like', "G2M", "G1S", "AC-like", "NPC-like", "OC-like", "Lymph", "Myeloid", "Vascular", "Oligo")
ccl <- data.frame(cty = row_order, 
                  group = c(rep("0tmr", 6), rep("1tme",4)))
pltt <- pltt[, row_order]
ht_nei <- Heatmap(as.matrix(pltt), 
                        name = "PCC", cluster_rows = T, cluster_columns = F, show_column_names = T, show_row_names = T,
                        col = col_fun, #rect_gp = gpar(),
                        #column_names_side = "top", column_names_gp = gpar(fontsize = 12), column_names_rot = 90,
                        width = ncol(pltt) * unit(6, "mm"), height = nrow(pltt) * unit(8, "mm"),
                        column_split  = ccl$group, 
                        row_title = " ", 
                        row_dend_width = unit(5, "mm"),
                  row_dend_gp = gpar(lwd = 1.5, col = "black"),
                        #left_annotation = hl_IDHmut, right_annotation = path_anno, top_annotation = ht,
                        show_heatmap_legend = T, 
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.rect(x = x, y = y, width = width, height = height, 
                              gp = gpar(col = "#969696", fill = NA, lwd = 1.5))
                  },
                        border = T, border_gp = gpar(lwd = 1.5, col = "black"),
                  column_gap = unit(2, "mm"))
ht_nei

pdf("./figures/Fig5C_1130_cell_neighborhood_niches.pdf", width = 6, height = 7)
draw(ht_nei, heatmap_legend_side = "right", 
     annotation_legend_side = "right")
dev.off()

# Fig 5D Example patch visualization ----
## Up panel: H&E images of the selected example patch
## Bottom panel: cell type annotation of the selected example patch

# Fig 5E 
## dotplots for all ---- 
res_cor_cut0.05 <- do.call(rbind, lapply(cellstates, get_scc_res, 0.05))
# res_cor_cut0.02 <- do.call(rbind, lapply(cellstates, get_scc_res, 0.02))
orderID <- c("n_MES", "n_G2M", "n_G1S", "n_AC", "n_NPC", "n_OC", "n_Lym", "n_Mye", "n_Vas", "n_Olig")
plt_dot <- res_cor_cut0.05[!(res_cor_cut0.05$Cell_State %in% c("n_lowQ", "n_Unknown")), ]
plt_dot$dotSize <- ifelse(plt_dot$PCC_pvalue < 0.1, -log10(plt_dot$PCC_pvalue), 1)
library(tidyverse)
plot_2.data <- plt_dot[which(plt_dot$dotSize > 1), ]
plot_3.data <- plt_dot[which(plt_dot$dotSize == 1), ]
coMu.plot <- ggplot() +
  theme_classic()
coMu.plot <- coMu.plot + geom_vline(xintercept = unique(plt_dot$Cutoff), linetype = 2, color = "grey", size = 0.3)
coMu.plot <- coMu.plot + geom_hline(yintercept = unique(plt_dot$Cell_State), linetype = 2, color = "grey", size = 0.3)

coMu.plot <- coMu.plot + geom_point(data = plt_dot, aes(y = Cell_State, x = Cutoff, size = dotSize, fill = PCC), shape = 21, color = "white") + ylab(NULL) + xlab(NULL) + ggtitle(NULL)
coMu.plot
coMu.plot <- coMu.plot + geom_point(data = plot_2.data, aes(y = Cell_State, x = Cutoff, size = dotSize, fill = PCC), shape = 21) + ylab(NULL) + xlab(NULL) + ggtitle(NULL)
coMu.plot
coMu.plot <- coMu.plot + theme(
  panel.background = element_rect(fill = "transparent", color = "black"), plot.margin = unit(c(1, 4, 1, 1), "lines"), plot.title = element_text(size = 24, vjust = 0.5, hjust = 0.5, face = "bold.italic"),
  text = element_text(size = 24, vjust = 1.4, hjust = 0.5, face = "bold"), legend.key.width = unit(1, "cm"), legend.key.height = unit(0.5, "cm"), legend.position = "bottom",
  legend.text = element_text(size = 12, face = "plain"), axis.text.y = element_text(size = 16, vjust = 0.5, hjust = 1, face = "plain", color = "black"), legend.title = element_text(size = 12, vjust = 0.5, hjust = 0, face = "plain"),
  axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust = 1, face = "plain", color = "black"), axis.title.x = element_text(size = 20, vjust = 0, hjust = 0.5, face = "plain", color = "black"),
  axis.title.y = element_text(size = 20, hjust = 0.5, vjust = 2, face = "plain", color = "black"),
  strip.text = element_text(size = 18, face = "bold", vjust = 0.5, hjust = 0.5), strip.background = element_rect(colour = "black", fill = gg_color_hue(3))
)
coMu.plot
coMu.plot <- coMu.plot + scale_shape_manual(
  name = NULL, values = c(A_green = 15, B_red = 17, C_grey = 16),
  labels = c(A_green = "Mutual exclusion", B_red = "Co-occurrence", C_grey = "Not significant"),
  guide = guide_legend(override.aes = list(size = 4), nrow = 4), na.translate = F
)
coMu.plot <- coMu.plot + scale_size(name = "Fisher's test\n   p-value", breaks = c(1, 2, 4, 6, 8), labels = expression(0.01, 10^-2, 10^-4, 10^-6, 10^-8), guide = guide_legend(nrow = 2), range = c(3, 6))
coMu.plot <- coMu.plot + scale_fill_gradientn(name = "        Odds ratio\n(log10 transformed)", colours = c("#276419", "#4d9221", "#7fbc41", "white", "#de77ae", "#c51b7d", "#8e0152"), limits = c(-1, 1), breaks = seq(-1, 1, 0.5))
coMu.plot <- coMu.plot + scale_x_discrete(position = "bottom") + scale_y_discrete(limits = rev(orderID))
coMu.plot
plotAll <- rbind(ggplotGrob(coMu.plot), size = "first")
plotAll
ggsave(file = "./figures/Fig5E_cor_spatialDensity.pdf", plot = plotAll, bg = "white", width = 6, height = 12, units = "cm", dpi = 600)

## mes-like plots ----
res_bin_ave_IR %>%
  filter(Cell_State == "n_MES") %>%
  ggplot() +
  theme_classic() +
  geom_point(aes(x = nCount_Ave, y = bins), size = 2, alpha = 1, color = "#4daf4a") +
  stat_smooth(aes(x = nCount_Ave, y = bins), method = "lm", color = "black", se = F) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(
    panel.background = element_rect(fill = "transparent", color = "transparent"), plot.margin = unit(c(2, 2, 2, 2), "lines"),
    plot.title = element_text(size = 34, vjust = 0.5, hjust = 0.5, face = "bold.italic", color = "transparent"), text = element_text(size = 14, face = "bold"),
    legend.key.width = unit(0.6, "cm"), legend.key.height = unit(0.6, "cm"), legend.position = "right", legend.text = element_text(size = 14, hjust = 0, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"), axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_text(size = 14, face = "plain", color = "black"), axis.title.y = element_text(size = 14, face = "plain", color = "black")
  ) +
  xlab("Binned patch cell states density") +
  ylab("GD prediction score")
ggsave("figures/Fig5E_density_GDscore_n_MES.pdf", width = 3.43, height = 3)

## Lym-like plots ----
res_bin_ave_IR %>%
  filter(Cell_State == "n_Lym") %>%
  ggplot() +
  theme_classic() +
  geom_point(aes(x = nCount_Ave, y = bins), size = 2, alpha = 1, color = "#d53f88") +
  stat_smooth(aes(x = nCount_Ave, y = bins), method = "lm", color = "black", se = F) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(
    panel.background = element_rect(fill = "transparent", color = "transparent"), plot.margin = unit(c(2, 2, 2, 2), "lines"),
    plot.title = element_text(size = 34, vjust = 0.5, hjust = 0.5, face = "bold.italic", color = "transparent"), text = element_text(size = 14, face = "bold"),
    legend.key.width = unit(0.6, "cm"), legend.key.height = unit(0.6, "cm"), legend.position = "right", legend.text = element_text(size = 14, hjust = 0, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"), axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_text(size = 14, face = "plain", color = "black"), axis.title.y = element_text(size = 14, face = "plain", color = "black")
  ) +
  xlab("Binned patch cell states density") +
  ylab("GD prediction score")
ggsave("figures/Fig5E_density_GDscore_n_Lym.pdf", width = 3.43, height = 3)
