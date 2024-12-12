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

# Fig 5E 
