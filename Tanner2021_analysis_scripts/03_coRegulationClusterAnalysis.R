# coRegulation Cluster Analysis
# calculate and Visualize co-regulation
# MH -> Tanner2021
#######################################

# load packages/dependencies
library(data.table)
library(pheatmap)
library(corrplot)
library(gplots)
library(Biobase)
library(factoextra)
library(NbClust)

# load data
cRes = fread("differentialExpressionResults_combined_proteinlevel.csv")
poi = fread("proteins_of_interest_070521.tsv")
poi[, protein_of_interest:=TRUE]
names(poi) = gsub(" ", "_", names(poi))

# select proteins regulated significantly in at least 1 comparison
diffprots_n1 = cRes[abs(log2_fold_change_protein)>=1 &
    p_value_BHadj_protein <= 0.05, unique(Protein.Group)]

# Which proteins of interest are among them?
poi[which(poi$Protein.Group %in% diffprots_n1)]

# factorize compartment and comparison columns for corrent ordering in plots
cRes[, compartment:=factor(compartment, levels = c("Lung", "BALF"))]
cRes[, comparison:=factor(comparison, levels = c("V/B", "BTH/B",
"DEX/B", "TH/B"))]

# obtain FC profiles for all these proteins
fcs = dcast(cRes[Protein.Group %in% diffprots_n1], Protein.Group~compartment+comparison, value.var = "log2_fold_change_protein")
fcs.m = as.matrix(fcs[, 2:ncol(fcs), with = F])
row.names(fcs.m) = fcs$Protein.Group
fcs.m[is.na(fcs.m)] = 0

# Visualize response profiles in Heatmap
hm = pheatmap(fcs.m,
         cluster_cols = F,
         show_rownames = F,
         main = "LTr04 Fibrosis model response profiles",
         # cutree_rows = 9,
         color = colorRampPalette(c("blue","white","red"))(100),
         gaps_col = c(4))
dev.copy(pdf, "03-01 Fibrosis model regulation profile heatmap.pdf", width = 4, height = 8)
dev.off()

# Determine the optimal number of Protein clusters via 2 methods
set.seed(123)
fviz_nbclust(fcs.m, FUNcluster = hcut, nboot = 20, method = "gap_stat")+
  labs(subtitle = "Gap statistic method")
ggsave("03-02_cluster_number_determination_gap_stat_method.pdf", height = 3, width = 3)

fviz_nbclust(fcs.m, FUNcluster = hcut, nboot = 20, method = "silhouette")+
  labs(subtitle = "Silhouette score")
ggsave("03-02_cluster_number_determination_gap_silhouette_method.pdf", height = 3, width = 3)

# Gap statistical method chosen, N = 8 clusters

## Calculate clusters and highlight in heatmap

# prepare annotations
ann_row = data.table("rn" = rownames(fcs.m))
ann_row = merge(ann_row, poi[, .(Protein.Group, protein_of_interest, GO_term_association)],
                by.x = "rn", sort = F, by.y = "Protein.Group", all.x = T, all.y = F)
ann_row[is.na(ann_row)] = 0
ann_row = as.data.frame(ann_row)
rownames(ann_row) = ann_row$rn
ann_row$rn = NULL
str(ann_row)

ann_col = data.frame("rn" = colnames(fcs.m))
ann_col$compartment = sapply(ann_col$rn, FUN = function(x){unlist(strsplit(x, split = "_"))[1]})
ann_col$comparison = sapply(ann_col$rn, FUN = function(x){unlist(strsplit(x, split = "_"))[2]})
rownames(ann_col) = ann_col$rn
ann_col$rn = NULL

# get cluster membership from clustering denrogram:
ann_row$cluster_id_euclid=factor(cutree(hm$tree_row, k=8))

# re-draw heatmap with these annotations
pheatmap(fcs.m,
              cluster_cols = F,
              show_rownames = F,
              annotation_row = ann_row,
              main = "LTr04 Fibrosis model response profiles",
              cutree_rows = 8,
              color = colorRampPalette(c("blue","white","red"))(100),
              gaps_col = c(4))

dev.copy(pdf, "03-03 LTr04 Fibrosis model regulation profile heatmap with clusters annotated.pdf", width = 7, height = 8)
dev.off()

# Plot protein fold-change line graphs per clusters
ann_row$rn = row.names(ann_row)
cRes = merge(cRes, ann_row, by.x = "Protein.Group", by.y = "rn")

cRes[, compartment_comparison:=paste(compartment,comparison, sep = "_"), .(compartment,comparison)]
unique(cRes$compartment_comparison)
cRes[, compartment_comparison:=factor(compartment_comparison, levels = 
                                          c("Lung_V/B", "Lung_BTH/B", "Lung_DEX/B", "Lung_TH/B",
                                            "BALF_V/B", "BALF_BTH/B", "BALF_DEX/B", "BALF_TH/B"))]

ggplot(cRes, aes(compartment_comparison, log2_fold_change_protein, color = cluster_id_euclid)) +
  geom_line(aes(group = Protein.Group), alpha = 0.6, size = 0.1) +
  facet_wrap(~paste("Cluster",cluster_id_euclid)) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept = 0, col = "black", lty = 2) +
  geom_vline(xintercept = 4.5, lty = 2) +
  theme(legend.position = "none")
  # geom_line(data = cRes[protein_of_interest.y == TRUE],
  #          aes(group = Protein.Group), alpha = 1, col = "red", size =2)
ggsave("03-04 LTr04 Fibrosis model regulation profile clusters fc profiles.pdf", width = 6, height = 6)

# Write out data.table with the results
cRes.clustRes = merge(fcs.m, ann_row, by=0)
cRes.clustRes$rn = NULL
fwrite(cRes.clustRes, "03-05_Clustering_result_table.csv")


# Plot Similarity matrix
fcs.m.dist = dist(fcs.m, method = "euclidean")
fcs.m.dist.m = as.matrix(fcs.m.dist)

fcs.m.dist.m.sim = 1/(1+fcs.m.dist.m) #Convert from distance to similarity

corrplot(fcs.m.dist.m.sim, 
         tl.pos = "n",
         order = "hclust",
         addrect = 8,
         method = "color")
dev.copy(pdf, "03-06 LTr04 Fibrosis model response profile similarity corrplot.pdf", width = 12, height = 12)
dev.off()

hm2 = heatmap.2(fcs.m.dist.m.sim,
          main = "LTr04 Fibrosis model response profile similarity heatmap",
          trace = "none",
          col = colorRampPalette(colors = c("Blue","white","Red"))(100))

dev.copy(pdf, "03-07 LTr04 Fibrosis model response profile similarity heatmap.pdf", width = 12, height = 12)
dev.off()

# Clarify value distribution before filtering
hist(fcs.m.dist.m.sim)

# Convert correlation matrix to tabular format for 
# potential downstream exploration
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.table(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = (cormat)[ut]
  )
}

fcs.m.dist.m.sim.dt = flattenCorrMatrix(fcs.m.dist.m.sim)

# Export for visualization in e.g. Cytoscape
fwrite(fcs.m.dist.m.sim.dt, "03-08 coRegNetwork_Euklidean_n1_all.csv")
