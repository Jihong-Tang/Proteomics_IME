library(tidyverse)

# 1- DNA methylation heatmap ---- 

get_Subtype_DE <- function(sid, df_DE){
  df1 <- df_DE[sid, ]
  df2 <- df_DE[!(row.names(df_DE) %in% sid), ]
  do_DE <- function(i, df1, df2){
    tmp <- wilcox.test(df1[, i], df2[, i])
    pvalue <- tmp$p.value
    tstat <- tmp$statistic
    fc <- mean(df1[, i]) - mean(df2[, i])
    return(data.frame(pid = colnames(df1)[i], pvalue, fc, tstat))
  }
  res_DE <- do.call(rbind, lapply(1:ncol(df1), function(x)do_DE(x, df1, df2)))
  res_DE$fdr <- p.adjust(res_DE$pvalue, method = 'fdr')
  return(res_DE)
}

ssMwwGST <- function(geData, geneSet, nCores = 8){
  library(yaGST)
  library(doMC)
  
  means <- rowMeans(geData)
  sds <- apply(geData, 1, sd)
  
  registerDoMC(nCores)
  ans <- foreach(ss = 1:ncol(geData)) %dopar% {
    currentSample <- (geData[, ss] - means)/sds
    rankedList <- sort(currentSample, decreasing = T)
    
    aMwwGST <- lapply(geneSet, function(x) mwwGST(rankedList, geneSet = x, minLenGeneSet = 15, alternative = "two.sided", verbose = F))
    aMwwGST <- aMwwGST[sapply(aMwwGST, length) != 0]
    tmp_NES <- sapply(aMwwGST, function(x) x$log.pu)
    tmp_pValue <- sapply(aMwwGST, function(x) x$p.value)
    
    ans <- list(tmp_NES = tmp_NES, tmp_pValue = tmp_pValue)
    print(ss)
    return(ans)
  }
  NES <- sapply(ans, function(x) x$tmp_NES)
  pValue <- sapply(ans, function(x) x$tmp_pValue)
  colnames(NES) <- colnames(pValue) <- colnames(geData)
  
  FDR <- t(apply(pValue, 1, function(x) p.adjust(x, method = "fdr")))
  res <- list(NES = NES, pValue = pValue, FDR = FDR)
  return(res)
}

library(tidyverse)
cl_3C_IDHmut <- read_delim("../../1025_MutA_ProteoCluster_Characterize_RNA/data/nmf_subtypes_5k_clV2_noncodel_rank4_Aug21.txt")

df.methyl.sigs <- read_delim("../../241009_MultiOmics_Classifier_forIDHmes/results/1008_ProteoNMF_M450K_signatures_W.txt")

df_CGGA_methyl <- readRDS("../../241009_MultiOmics_Classifier_forIDHmes/data/CGGA_Methylation_Data_103samples_Filtered_1225.rds")
df_MutA_methyl <- df_CGGA_methyl[, colnames(df_CGGA_methyl) %in% cl_3C_IDHmut$Cohort_ID]

df_DE <- data.frame(t(df_MutA_methyl))
df_DE <- df_DE[row.names(df_DE) %in% cl_3C_IDHmut$Cohort_ID, ]
#demo_methyl <- df_DE[seq(1,3), ]
#write.csv(demo_methyl, "./results/demo_methyl.csv")
cl_DE <- cl_3C_IDHmut[cl_3C_IDHmut$Cohort_ID %in% row.names(df_DE), ]

sid_ProteoNMF1 <- cl_DE$Cohort_ID[cl_DE$NMF_Cluster == "ProteoNMF1"]
res_DE_ProteoNMF1 <- get_Subtype_DE(sid_ProteoNMF1, df_DE)

sid_ProteoNMF2 <- cl_DE$Cohort_ID[cl_DE$NMF_Cluster == "ProteoNMF2"]
res_DE_ProteoNMF2 <- get_Subtype_DE(sid_ProteoNMF2, df_DE)

sid_ProteoNMF3 <- cl_DE$Cohort_ID[cl_DE$NMF_Cluster == "ProteoNMF3"]
res_DE_ProteoNMF3 <- get_Subtype_DE(sid_ProteoNMF3, df_DE)

sid_ProteoNMF4 <- cl_DE$Cohort_ID[cl_DE$NMF_Cluster == "ProteoNMF4"]
res_DE_ProteoNMF4 <- get_Subtype_DE(sid_ProteoNMF4, df_DE)

save(res_DE_ProteoNMF1_methyl, res_DE_ProteoNMF2_methyl, 
     res_DE_ProteoNMF3_methyl, res_DE_ProteoNMF4_methyl, file = "results/1009_res_DE_ProteoNMF_methylation_signatures.Rdata")
load("../../241009_MultiOmics_Classifier_forIDHmes/results/1009_res_DE_ProteoNMF_methylation_signatures.Rdata")

top <- 750

methyl_Proteo1_W <- res_DE_ProteoNMF1_methyl$pid[order(res_DE_ProteoNMF1_methyl$fc, decreasing = T)][1:top]
methyl_Proteo2_W <- res_DE_ProteoNMF2_methyl$pid[order(res_DE_ProteoNMF2_methyl$fc, decreasing = T)][1:top]
methyl_Proteo3_W <- res_DE_ProteoNMF3_methyl$pid[order(res_DE_ProteoNMF3_methyl$fc, decreasing = T)][1:top]
methyl_Proteo4_W <- res_DE_ProteoNMF4_methyl$pid[order(res_DE_ProteoNMF4_methyl$fc, decreasing = T)][1:top]

methyl_Proteo1_W <- res_DE_ProteoNMF1_methyl$pid[order(res_DE_ProteoNMF1_methyl$tstat, decreasing = T)][1:top]
methyl_Proteo2_W <- res_DE_ProteoNMF2_methyl$pid[order(res_DE_ProteoNMF2_methyl$tstat, decreasing = T)][1:top]
methyl_Proteo3_W <- res_DE_ProteoNMF3_methyl$pid[order(res_DE_ProteoNMF3_methyl$tstat, decreasing = T)][1:top]
methyl_Proteo4_W <- res_DE_ProteoNMF4_methyl$pid[order(res_DE_ProteoNMF4_methyl$tstat, decreasing = T)][1:top]


#list_cpgs <- c(df.methyl.sigs$ProteoNMF1, df.methyl.sigs$ProteoNMF2, df.methyl.sigs$ProteoNMF3, df.methyl.sigs$ProteoNMF4)
list_cpgs <- c(methyl_Proteo1_W, methyl_Proteo2_W, methyl_Proteo3_W, methyl_Proteo4_W)

plt_mutA <- df_MutA_methyl[list_cpgs, ]
plt_cl <- cl_3C_IDHmut[cl_3C_IDHmut$Cohort_ID %in% colnames(plt_mutA), c("Cohort_ID", "NMF_Cluster")]

plt_mutA <- plt_mutA[, plt_cl$Cohort_ID]
plt_mutA_scale <- apply(plt_mutA, 1, scale) %>% t()
colnames(plt_mutA_scale) <- colnames(plt_mutA)
sample_order_MutA <- plt_cl$Cohort_ID[order(plt_cl$NMF_Cluster, decreasing = F)]

library(ComplexHeatmap)
library(dendextend)
library(circlize)
col_fun = colorRamp2(seq(-1, 1, 0.1), colorRampPalette(c('#29419E','#698DC9',"grey95",'#F67172','#DC2B18'))(21))
top_anno <- HeatmapAnnotation(NMF_Cluster = plt_cl$NMF_Cluster,
                              
                              col = list(
                                NMF_Cluster = c("ProteoNMF1" = "#33a02c", "ProteoNMF2" = '#e31a1c', "ProteoNMF3" = '#ff7f00', "ProteoNMF4" = '#1f78b4')
                              ),
                              annotation_label = c("Proteomics subgroup"), 
                              annotation_name_gp= gpar(fontsize = 8),
                              simple_anno_size = unit(.25, "cm"),
                              annotation_legend_param = list(direction = "horizontal"))

ht_methyl <- Heatmap(as.matrix(plt_mutA_scale), 
                 name = "DNA methylation", 
                 height = unit(4, "cm"), width = unit(6, "cm"),
                 column_split = plt_cl$NMF_Cluster, column_order = sample_order_MutA,
                 show_row_names = F, show_row_dend = F, top_annotation = top_anno,
                 show_column_names = F,
                 row_names_gp = gpar(fontsize = 10),
                 column_names_gp = gpar(fontsize = 10),
                 col = col_fun, 
                 cluster_rows = F, cluster_columns = F, column_title =' ',
                 use_raster = F
)                             
ht_methyl                             
pdf("./figures/Fig6/1130_methyl_heatmap.pdf", width = 5, height = 4)
draw(ht_methyl, heatmap_legend_side = "right", 
     annotation_legend_side = "bottom")
dev.off()


