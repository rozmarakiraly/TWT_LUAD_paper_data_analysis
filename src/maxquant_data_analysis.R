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
# in this part you calculate which proteins fall into each category
# and create a new column the specifies this
NvsT_df_imputed%>%
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