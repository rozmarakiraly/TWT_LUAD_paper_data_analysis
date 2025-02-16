library(tidyverse)
library(MSnSet.utils)

# ..................................................................... ####
# read in data ####

# first, read the protein expression table
setwd("..")
mq_df <- read.csv("data/MQ_input_data.csv", row.names = "Protein.IDs")
# log transform data
mq_df <- log2(mq_df)

# check the distribution of your data
boxplot(mq_df)

# second, read the phenotype data table (all info/groupings of the samples)
phenotype_df <- read.csv("data/phenotype_input.csv", row.names = "X")

# next, create an MSnSet object (ecol should specify columns with values)
all_data <- readMSnSet2(mq_df, ecol = names(mq_df)[sapply(mq_df, is.numeric)])

# finally, add phenotype information to the object
pData(all_data) <- phenotype_df

# ..................................................................... ####
# create PCA plot ####
# it adds 50% confidence ellipse
phenotype_df$Color <- phenotype_df$Group
plot_pca(all_data, phenotype = "Group")
?plot_pca
# ..................................................................... ####
# create Heatmap ####
complex_heatmap(all_data, heatmap_type = "sample_correlation",
                anno_column = "Group",
                anno_column_colors = list(Group = c("#414DBE", "#BEB241")),
                anno_row_colors = list(Group = c("#414DBE", "#BEB241")),
                heatmap_args = list(col = c("white","#d7b5d8","#253494")))

# ..................................................................... ####
# limma ####
# ... using lm to check whether a variable has effect on the data (batch, age etc.) ####
lm_batches <- limma_gen(all_data, model.str = "~ Batch", coef.str = "Batch")
# covariates can be added to the model (e.g. batch effect on an other variable)

# ... limma on imputed data ####
library(imputeLCMD)
imputed_data <- impute(all_data, method = "QRILC")
sum(is.na(exprs(imputed_data)))

# moderated t-tests for comparing normal and tumor
NvsT_df_imputed <- limma_a_b(eset = imputed_data, model.str = "~ Group", coef.str = "Group")

# plot histogram of p-values, it should NOT be a uniform distribution
hist(NvsT_df_imputed$P.Value)

# number of significant differences
sum(NvsT_df_imputed$adj.P.Val < 0.05, na.rm = T)

# ... plot heatmap for significant proteins ####
# create a new MSnSet with only the data for the heatmap, just to make things easier
filtered_mq_df <- mq_df[row.names(mq_df) %in% row.names(NvsT_df)[NvsT_df$adj.P.Val < 0.05 & !is.na(NvsT_df$adj.P.Val)],]
filtered_mq_df <- as.data.frame(t(scale(t(filtered_mq_df))))
filtered_mq_df[is.na(filtered_mq_df)] <- 0
significant_data <- readMSnSet2(filtered_mq_df, ecol = names(mq_df)[sapply(mq_df, is.numeric)])
pData(significant_data) <- phenotype_df

complex_heatmap(significant_data, heatmap_type = "expression",
                heatmap_args = list(col = circlize::colorRamp2(
                  breaks = c(min(exprs(significant_data), na.rm = TRUE), 
                             -1.5, 1.5, # add color limits
                             max(exprs(significant_data), na.rm = TRUE)), 
                  colors = c("orange", "orange", "darkblue", "darkblue"))),
                anno_column = c("Group","Grade"),
                anno_column_colors = list(Group = c("blue", "red")))

# ..................................................................... ####
# create a volcano plot ####
NvsT_df_imputed %>%
  mutate(point_color = case_when(
    adj.P.Val < 0.05 & logFC < 0 ~ "down",
    adj.P.Val < 0.05 & logFC > 0 ~ "up",
    TRUE ~ "NS")) -> NvsT_df

v1 <- plot_volcano(df = NvsT_df, logFC = "logFC", 
                   pvals = "adj.P.Val", sig_threshold = 0.05,
                   point_args = list(mapping = aes(color = point_color)))

v1 + scale_color_manual(values = c("#5555ff", "red3", "lightgrey"), 
                        breaks = c("down", "up", "NS")) +
  theme(legend.position = "none") # do not show legend

# ..................................................................... ####
# GSEA and related plots ####
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# gsea_input includes the protein name (transformed to ENSEMBL) and B columns from NvsT_df_imputed
# ID matching was done using mapIDs()
data_for_gsea <- read.csv(file = "gsea_input.csv", header = F)

gsea_input <- setNames(data_for_gsea$V2, test_data_gsea$V1)

system.time(
  gseGO_res <- gseGO(
    geneList = gsea_input,
    ont = "BP",
    OrgDb = "org.Hs.eg.db",
    keyType = "ENSEMBL",
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 500,
    eps = 1e-10,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE,
    seed = FALSE,
    nPermSimple = 100000,
    by = "fgsea"
  )
)
gsea_results <- gseGO_res@result

gsea_simp <- clusterProfiler::simplify(gseGO_res)

gsea_simp <- pairwise_termsim(gsea_simp)

gsea_simp_results <- gsea_simp@result

color_scale <- scale_colour_gradient2(low = "blue", high = "red", midpoint = 0)

png(filename = "dotplot.png", 
    type = "cairo",
    units = "in", 
    width = 6,
    height = 7,
    res = 300)
dotplot(gsea_simp,
        showCategory = 20, # number of top GOs to include
        size = "Count",
        font.size = 8,
        label_format = 40,
        color = "NES",
        # foldChange = gsea_input,
        title = "Top 20 Gene Sets") + 
  theme(legend.text = element_text(size = 6), legend.title = element_text(size = 6), legend.key.size = unit(0.25, "in")) +
  color_scale
# coord_flip() +
# theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

png(filename = "treeplot.png", 
    type = "cairo",
    units = "in", 
    width = 10,
    height = 8,
    res = 300)
treeplot(gsea_simp, # foldChange = geneList,
         showCategory = 20,
         color = "NES",
         cex_category = 1,
         split = T,
         # fontsize = 5,
         cluster.params = list(method = "ward.D", n = 8, color = NULL, # group colors
                               label_words_n = 4, # number of words to display on tip
                               label_format = 20), # text box for labels
         hilight.params = list(hilight = T, align = "right"),
         clusterPanel.params = list(clusterPanel = "dotplot", # the type of plot to make
                                    pie = "equal", legend_n = 3,
                                    hexpand = 5),
         offset.params = list(bar_tree = rel(2),
                              tiplab = rel(1),
                              extend = 0.3, # cluster label line length
                              hexpand = 0.2)) + # size of dotplot / heatmap area
  color_scale
dev.off()

png(filename = "emapplot.png", 
    type = "cairo",
    units = "in", 
    width = 12,
    height = 12,
    res = 300)
emapplot(gsea_simp,
         showCategory = 75,
         color = "NES",
         shadowtext = T,
         repel = T,
         custom_node_colors = color_scale,
         node_label = "category", # other option is group for clusters
         layout.params = list(layout = NULL, coords = NULL),
         edge.params = list(show = TRUE, min = 0.2),
         cex.params = list(category_node = 1.5, # relative node size
                           category_label = 0.75, # label size
                           line = 1.5), # edge thickness
         hilight.params = list(category = c(), # name of nodes to highlight
                               alpha_hilight = 1,
                               alpha_no_hilight = 0.5),
         cluster.params = list(cluster = FALSE, method = stats::kmeans, n = NULL, legend = FALSE,
                               label_style = "shadowtext", label_words_n = 4, label_format = 30))
dev.off()

# write.csv(gsea_simp_results, file = "gsea_results_top_75_gs.csv")