library(tidyverse)

# 
codex.lym <- readxl::read_xlsx("../Fig_codes/data/fig6/Ini_Lymphocyte_Manual_Count/0815_PLC_Lymphocyte_Manual_Count.xlsx", sheet = 2)

list.region <- unique(codex.lym$Region[codex.lym$`Tcell (CD3E)` != 0])

region <- list.region[1]

get.region <- function(region){
  
  cd3 <- codex.lym$`Tcell (CD3E)`[codex.lym$Region == region]
  cd4 <- codex.lym$CD4T[codex.lym$Region == region]
  cd8 <- codex.lym$CD8T[codex.lym$Region == region]
  cddn <- codex.lym$`Double Negative`[codex.lym$Region == region]
  
  CD4pd1 <- codex.lym$`CD4+PD1+`[codex.lym$Region == region]
  CD4tox <- codex.lym$`CD4+TOX+`[codex.lym$Region == region]
  CD4gzmk <- codex.lym$`CD4+GZMK+`[codex.lym$Region == region]
  CD4pd1tox <- codex.lym$`CD4+PD1+TOX+`[codex.lym$Region == region]
  CD4pd1gzmk <- codex.lym$`CD4+PD1+GZMK+`[codex.lym$Region == region]
  CD4gzmktox <- codex.lym$`CD4+GZMK+TOX+`[codex.lym$Region == region]
  CD4all <- codex.lym$`CD4+GZMK+TOX+PD1+`[codex.lym$Region == region]
  
  CD8pd1 <- codex.lym$`CD8+PD1+`[codex.lym$Region == region]
  CD8tox <- codex.lym$`CD8+TOX+`[codex.lym$Region == region]
  CD8gzmk <- codex.lym$`CD8+GZMK+`[codex.lym$Region == region]
  CD8pd1tox <- codex.lym$`CD8+PD1+TOX+`[codex.lym$Region == region]
  CD8pd1gzmk <- codex.lym$`CD8+PD1+GZMK+`[codex.lym$Region == region]
  CD8gzmktox <- codex.lym$`CD8+GZMK+TOX+`[codex.lym$Region == region]
  CD8all <- codex.lym$`CD8+GZMK+TOX+PD1+`[codex.lym$Region == region]
  
  df.region <- data.frame(region = rep(region, cd3), 
                          CD3 = rep(1, cd3),
                          CD4 = c(rep(1, cd4), rep(0, cd3 - cd4)),
                          CD8 = c(rep(0, cd4), rep(1, cd8), rep(0, cd3 - cd4 - cd8))
                          
  )
  pd1 <- c(rep(1, CD4all), rep(1, CD4pd1tox), rep(1, CD4pd1gzmk), rep(0, CD4gzmktox),
           rep(1, CD4pd1), rep(0, CD4tox), rep(0, CD4gzmk), 
           rep(0, cd4 - CD4all - CD4pd1 - CD4tox - CD4gzmk - CD4pd1tox - CD4pd1gzmk - CD4gzmktox),
           rep(1, CD8all), rep(1, CD8pd1tox), rep(1, CD8pd1gzmk), rep(0, CD8gzmktox),
           rep(1, CD8pd1), rep(0, CD8tox), rep(0, CD8gzmk),
           rep(0, cd8 - CD8all - CD8pd1 - CD8tox - CD8gzmk - CD8pd1tox - CD8pd1gzmk - CD8gzmktox),
           rep(0, cddn)
  ) 
  tox <- c(rep(1, CD4all), rep(1, CD4pd1tox), rep(0, CD4pd1gzmk), rep(1, CD4gzmktox),
           rep(0, CD4pd1), rep(1, CD4tox), rep(0, CD4gzmk), 
           rep(0, cd4 - CD4all - CD4pd1 - CD4tox - CD4gzmk - CD4pd1tox - CD4pd1gzmk - CD4gzmktox),
           rep(1, CD8all), rep(1, CD8pd1tox), rep(0, CD8pd1gzmk), rep(1, CD8gzmktox),
           rep(0, CD8pd1), rep(1, CD8tox), rep(0, CD8gzmk),
           rep(0, cd8 - CD8all - CD8pd1 - CD8tox - CD8gzmk - CD8pd1tox - CD8pd1gzmk - CD8gzmktox),
           rep(0, cddn)
  )
  gzmk <- c(
    rep(1, CD4all), rep(0, CD4pd1tox), rep(1, CD4pd1gzmk), rep(1, CD4gzmktox),
    rep(0, CD4pd1), rep(0, CD4tox), rep(1, CD4gzmk), 
    rep(0, cd4 - CD4all - CD4pd1 - CD4tox - CD4gzmk - CD4pd1tox - CD4pd1gzmk - CD4gzmktox),
    rep(1, CD8all), rep(0, CD8pd1tox), rep(1, CD8pd1gzmk), rep(1, CD8gzmktox),
    rep(0, CD8pd1), rep(0, CD8tox), rep(1, CD8gzmk),
    rep(0, cd8 - CD8all - CD8pd1 - CD8tox - CD8gzmk - CD8pd1tox - CD8pd1gzmk - CD8gzmktox),
    rep(0, cddn)
  )
  
  df.region$PD1 <- pd1
  df.region$TOX <- tox
  df.region$GZMK <- gzmk
  return(df.region)
}

df.region <- do.call(rbind, lapply(list.region, get.region))

rg.label <- codex.lym[, c("Region", "Group")]; colnames(rg.label) <- c("region", "group")
df.region <- merge(df.region, rg.label, by = "region")

df.region %>% group_by(group) %>% summarise(
  CD3 = sum(CD3),
  CD4 = sum(CD4),
  CD8 = sum(CD8),
  PD1 = sum(PD1),
  TOX = sum(TOX),
  GZMK = sum(GZMK))

## sub2 upset plots ---- 
library(ComplexHeatmap)
library(tidyverse)

### plt1 plc ---- 
up_plt1 <- df.region[df.region$group %in% c("PLC", "PLC_Near"),]
ctys <- c("PD1", "GZMK")
ctys <- colnames(up_plt1)[2:7]
up_plt1$cID <- paste0("cell", seq(1, nrow(up_plt1)))

list_ctys <- lapply(ctys, function(x){
  tmp = up_plt1$cID[up_plt1[,x] == 1]
  out = tmp[!is.na(tmp)]})
names(list_ctys) <- ctys
#color1<- c("CT" = "#e31a1c", "IT" = "#1f78b4",  "LE" ="#ff7f00", "MVP" = "#33a02c", "PAN" = "#6a3d9a")
m <- make_comb_mat(list_ctys, mode = "distinct")

UpSet(m)
ss = set_size(m)
cs = comb_size(m)
ht <- UpSet(m, comb_order = order(-comb_degree(m), -cs),
            comb_col = c("#885e80"), 
            top_annotation = HeatmapAnnotation(
              " " = anno_barplot(cs, 
                                 ylim = c(0, max(cs)*1.1),
                                 border = FALSE, 
                                 gp = gpar(fill = "#ee837a"), 
                                 height = unit(2, "cm")
              ), 
              annotation_name_side = "left", 
              annotation_name_rot = 90),
            left_annotation = rowAnnotation(
              set_name = anno_text(set_name(m), 
                                   location = .5, 
                                   just = "center",
                                   width = max_text_width(set_name(m)) + unit(4, "mm"))
            ), 
            right_annotation = rowAnnotation(
              "Cells" = anno_barplot(ss, 
                                     baseline = 0,
                                     axis_param = list(
                                       at = c(0, 45, 90, 135, 181),
                                       labels = c(0, 45, 90, 135, 180),
                                       labels_rot = 0),
                                     border = FALSE, 
                                     gp = gpar(fill = "#5caeae"), 
                                     width = unit(2, "cm")
              )
            ),
            show_row_names = F
)

pdf("figures/Fig6/0816_lym_upset_plc.pdf", width = 3.5, height = 2.3)
ht = draw(ht)
od = column_order(ht)
decorate_annotation(" ", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 10, col = "black"), rot = 45)
})
dev.off()

### plt2 vessel ----
up_plt1 <- df.region[df.region$group %in% c("Vessel"),]
ctys <- colnames(up_plt1)[2:7]
up_plt1$cID <- paste0("cell", seq(1, nrow(up_plt1)))
ctys <- c("PD1", "GZMK")

list_ctys <- lapply(ctys, function(x){
  tmp = up_plt1$cID[up_plt1[,x] == 1]
  out = tmp[!is.na(tmp)]})
names(list_ctys) <- ctys
#color1<- c("CT" = "#e31a1c", "IT" = "#1f78b4",  "LE" ="#ff7f00", "MVP" = "#33a02c", "PAN" = "#6a3d9a")
m <- make_comb_mat(list_ctys, mode = "distinct")

UpSet(m)
ss = set_size(m)
cs = comb_size(m)
ht <- UpSet(m, comb_order = order(-comb_degree(m), -cs),
            comb_col = c("#885e80"), 
            top_annotation = HeatmapAnnotation(
              " " = anno_barplot(cs, 
                                 ylim = c(0, max(cs)*1.1),
                                 border = FALSE, 
                                 gp = gpar(fill = "#ee837a"), 
                                 height = unit(2, "cm")
              ), 
              annotation_name_side = "left", 
              annotation_name_rot = 90),
            left_annotation = rowAnnotation(
              set_name = anno_text(set_name(m), 
                                   location = .5, 
                                   just = "center",
                                   width = max_text_width(set_name(m)) + unit(4, "mm"))
            ), 
            right_annotation = rowAnnotation(
              "Cells" = anno_barplot(ss, 
                                     baseline = 0,
                                     axis_param = list(
                                       at = c(0, 15, 30, 45),
                                       labels = c(0, 15, 30, 45),
                                       labels_rot = 0),
                                     border = FALSE, 
                                     gp = gpar(fill = "#5caeae"), 
                                     width = unit(2, "cm")
              )
            ),
            show_row_names = F
)

pdf("figures/Fig6/0816_lym_upset_vessel.pdf", width = 3.5, height = 2.3)
ht = draw(ht)
od = column_order(ht)
decorate_annotation(" ", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 10, col = "black"), rot = 45)
})
dev.off()


### plt3 plc_middle ---- 
up_plt1 <- df.region[df.region$group %in% c("PLC_Middle"),]
ctys <- colnames(up_plt1)[2:7]
up_plt1$cID <- paste0("cell", seq(1, nrow(up_plt1)))
ctys <- c("PD1", "GZMK")

list_ctys <- lapply(ctys, function(x){
  tmp = up_plt1$cID[up_plt1[,x] == 1]
  out = tmp[!is.na(tmp)]})
names(list_ctys) <- ctys
#color1<- c("CT" = "#e31a1c", "IT" = "#1f78b4",  "LE" ="#ff7f00", "MVP" = "#33a02c", "PAN" = "#6a3d9a")
m <- make_comb_mat(list_ctys, mode = "distinct")

UpSet(m)
ss = set_size(m)
cs = comb_size(m)
ht <- UpSet(m, comb_order = order(-comb_degree(m), -cs),
            comb_col = c("#885e80"), 
            top_annotation = HeatmapAnnotation(
              " " = anno_barplot(cs, 
                                 ylim = c(0, max(cs)*1.1),
                                 border = FALSE, 
                                 gp = gpar(fill = "#ee837a"), 
                                 height = unit(2, "cm")
              ), 
              annotation_name_side = "left", 
              annotation_name_rot = 90),
            left_annotation = rowAnnotation(
              set_name = anno_text(set_name(m), 
                                   location = .5, 
                                   just = "center",
                                   width = max_text_width(set_name(m)) + unit(4, "mm"))
            ), 
            right_annotation = rowAnnotation(
              "Cells" = anno_barplot(ss, 
                                     baseline = 0,
                                     axis_param = list(
                                       at = c(0, 15, 30, 45, 61),
                                       labels = c(0, 15, 30, 45, 60),
                                       labels_rot = 0),
                                     border = FALSE, 
                                     gp = gpar(fill = "#5caeae"), 
                                     width = unit(2, "cm")
              )
            ),
            show_row_names = F
)

pdf("figures/Fig6/0816_lym_upset_plc_middle.pdf", width = 3.5, height = 2.3)
ht = draw(ht)
od = column_order(ht)
decorate_annotation(" ", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 10, col = "black"), rot = 45)
})
dev.off()


### plt4 plc_far ---- 
up_plt1 <- df.region[df.region$group %in% c("PLC_Far"),]
ctys <- colnames(up_plt1)[2:7]
up_plt1$cID <- paste0("cell", seq(1, nrow(up_plt1)))
ctys <- c("PD1", "GZMK")

list_ctys <- lapply(ctys, function(x){
  tmp = up_plt1$cID[up_plt1[,x] == 1]
  out = tmp[!is.na(tmp)]})
names(list_ctys) <- ctys
#color1<- c("CT" = "#e31a1c", "IT" = "#1f78b4",  "LE" ="#ff7f00", "MVP" = "#33a02c", "PAN" = "#6a3d9a")
m <- make_comb_mat(list_ctys, mode = "distinct")

UpSet(m)
ss = set_size(m)
cs = comb_size(m)
ht <- UpSet(m, comb_order = order(-comb_degree(m), -cs),
            comb_col = c("#885e80"), 
            top_annotation = HeatmapAnnotation(
              " " = anno_barplot(cs, 
                                 ylim = c(0, max(cs)*1.1),
                                 border = FALSE, 
                                 gp = gpar(fill = "#ee837a"), 
                                 height = unit(2, "cm")
              ), 
              annotation_name_side = "left", 
              annotation_name_rot = 90),
            left_annotation = rowAnnotation(
              set_name = anno_text(set_name(m), 
                                   location = .5, 
                                   just = "center",
                                   width = max_text_width(set_name(m)) + unit(4, "mm"))
            ), 
            right_annotation = rowAnnotation(
              "Cells" = anno_barplot(ss, 
                                     baseline = 0,
                                     axis_param = list(
                                       at = c(0, 40, 80, 120, 161),
                                       labels = c(0, 40, 80, 120, 160),
                                       labels_rot = 0),
                                     border = FALSE, 
                                     gp = gpar(fill = "#5caeae"), 
                                     width = unit(2, "cm")
              )
            ),
            show_row_names = F
)

pdf("figures/Fig6/0816_lym_upset_plc_far.pdf", width = 3.5, height = 2.3)
ht = draw(ht)
od = column_order(ht)
decorate_annotation(" ", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 10, col = "black"), rot = 45)
})
dev.off()


## sub3 barplot for cd4 cd8 ratio ---- 

df.region %>% group_by(group) %>% summarise(
  CD3 = sum(CD3),
  CD4 = sum(CD4),
  CD8 = sum(CD8),
  PD1 = sum(PD1),
  TOX = sum(TOX),
  GZMK = sum(GZMK))

library(tidyverse)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}

specie <- c(rep("Near", 2), rep("Middle", 2), rep("Far",2))
condition <- rep(c("CD4", "CD8") , 3)
#value <- c(21, 520-21, 62, 283-62)

value <- c(97, 79, 26, 37, 47, 97)
data <- data.frame(specie,condition,value)

# Stacked + percent
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="fill", stat="identity", show.legend = F, width = .7) + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0,1,0.2), labels = scales::percent) + 
  scale_x_discrete(limits = c("Near", "Middle", "Far")) + 
  scale_fill_manual(values = c("CD4" = '#e41a1c', 
                               "CD8" = '#fffe54')) + 
  labs(x = "", y = "Percentage") + 
  theme_classic() + 
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.y=element_blank(),axis.text.x=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black')) + 
  coord_flip()

ggsave("figures/Fig6/0816_PLC_region_CD4CD8_comparison.pdf", width = 4, height = 2.5)
M <- as.table(rbind(c(97,26,47),c(176, 63, 144)))
chisq.test(M)


## sub3 barplot for gzmk pd1 ratio ---- 

df.region %>% group_by(group) %>% summarise(
  CD3 = sum(CD3),
  CD4 = sum(CD4),
  CD8 = sum(CD8),
  PD1 = sum(PD1),
  TOX = sum(TOX),
  GZMK = sum(GZMK))

library(tidyverse)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 150)[1:n]
}

specie <- c(rep("Near", 3), rep("Middle", 3), rep("Far",3))
condition <- rep(c("0PD1","1Both", "2GZMK") , 3)
#value <- c(21, 520-21, 62, 283-62)

value <- c(33,75,16,15,26,9,38,77,10)
data <- data.frame(specie,condition,value)

# Stacked + percent
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="fill", stat="identity", show.legend = F, width = .7) + 
  scale_y_continuous(expand = c(0,0),breaks = seq(0,1,0.2), labels = scales::percent) + 
  scale_x_discrete(limits = c("Near", "Middle", "Far")) + 
  scale_fill_manual(values = c("0PD1" = '#ffc0cb', 
                               '1Both' = "#c8e6f4",
                               "2GZMK" = '#75fbfd')) + 
  labs(x = "", y = "Percentage") + 
  theme_classic() + 
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.y=element_blank(),axis.text.x=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black')) + 
  coord_flip()

ggsave("figures/Fig5/1115_PLC_region_GZMK_PD1comparison.pdf", width = 4, height = 2.5)
M <- as.table(rbind(c(97,26,47),c(176, 63, 144)))
chisq.test(M)

# 5C cell neighborhood analysis ---- 
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
pdf("./figures/Fig3/1128_18MutA_tmr_TME_Interactions_PCC.pdf", width = 6, height = 7)
draw(ht_inter_pcc, heatmap_legend_side = "right", 
     annotation_legend_side = "right")
dev.off()
