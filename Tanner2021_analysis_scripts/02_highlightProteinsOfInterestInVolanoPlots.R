# Collect results and compare in heatmap and volcanos
# load packages
library(data.table)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(plotly)
library(crosstalk)

# collect results for protein differentials
cRes_lung = fread("differentialTestingResults_Lung.csv")
cRes_balf = fread("differentialTestingResults_BALF.csv")
cRes = rbind(cRes_lung, cRes_balf)

# Add proteins of interest provided by LT
poi = fread("proteins_of_interest_070521.tsv")

head(poi) #inspect
names(poi) = gsub(" ", "_", names(poi)) #replace whitespace in colnames
poi[, protein_of_interest:=TRUE]

# merge into differential testing results for highlighting
cRes = merge(cRes, poi[, .(Protein.Group, protein_of_interest, GO_term_association)], by = "Protein.Group", all.x = TRUE)

# set comparison factor levels to order plot panels as desired
cRes[, comparison:=factor(comparison, levels=c("V/B", "BTH/B", "DEX/B", "TH/B"))]

# refine protein identifiers
cRes[, n_proteins:=length(unlist(strsplit(Protein.Group, split = ";"))),Protein.Group]
cRes[, Protein.Names.Short:=gsub("\\(Bos;|_MOUSE", "", Protein.Names)]
cRes[, First.Protein.Name.Short:=unlist(strsplit(Protein.Names.Short, split = ";"))[1], Protein.Group]
cRes[, First.Protein.Name.Short:=unlist(strsplit(Protein.Names.Short, split = ";"))[1], Protein.Group]

# Remove residual precursor-level scores
cRes[, p_value:=NULL]
cRes[, log2_fold_change:=NULL]
cRes = unique(cRes)

fwrite(cRes, file = "differentialExpressionResults_combined_proteinlevel.csv")

# Make Volcano plots to compare
significance_threshold_p = 0.05
significance_threshold_fc = 2

# Plot selected volcanos
########################
# No highlight
volcanos = ggplot(cRes, aes(log2_fold_change_protein, -log10(p_value_BHadj_protein),
                            text = Protein.Names)) +
  geom_point(pch = 21, color = "grey", fill = "grey") +
  facet_wrap(~paste(compartment,comparison), nrow = 2) +
  theme_bw() +
  theme(legend.position = "right") +
  # scale_color_manual(values = c("darkgrey","red")) +
  geom_hline(yintercept = -log10(significance_threshold_p), col = "darkgrey", lty = 2) +
  geom_vline(xintercept = c(-log2(significance_threshold_fc),
                            log2(significance_threshold_fc)), col = "darkgrey", lty = 2) +
  ylab("-log10(BH-adj. p-value)") +
  xlab("log2(Fold-change)")
plot(volcanos)

ggsave("02_Volcanos_nohighlight.pdf", height = 7, width = 13)


# With proteins of interest highlighted
volcanos + geom_point(data = cRes[!is.na(GO_term_association)],
                      aes(pch = protein_of_interest,
                          col = GO_term_association), pch = 19, size = 2, alpha = 0.5) +
  scale_color_brewer(palette = "Spectral")
ggsave("02_Volcanos_highlight_allPOI.pdf", height = 7, width = 16)

volcanos + geom_point(data = cRes[!is.na(GO_term_association)],
                      aes(pch = protein_of_interest,
                          fill = GO_term_association,
                          col = GO_term_association), pch =21) +
  scale_fill_brewer(palette = "Spectral") +
  geom_text_repel(data = cRes[!is.na(GO_term_association)],
                  aes(label = Protein.Names, col = GO_term_association),
                  nudge_y = 5, max.iter = 20000, max.overlaps = 50)
  
ggsave("02_Volcanos_highlight_allPOI_labeled.pdf", height = 7, width = 16)

# Graphs too busy, do one by one:
for (term in unique(poi$GO_term_association)){
  print(term)
  volcanos + geom_point(data = cRes[GO_term_association == term],
                        aes(pch = protein_of_interest,
                            fill = GO_term_association,
                            col = GO_term_association), pch =21) +
    scale_fill_brewer(palette = "Spectral") +
    geom_text_repel(data = cRes[GO_term_association == term],
                    aes(label = Protein.Names, col = GO_term_association),
                    nudge_y = 5, max.iter = 20000, max.overlaps = 50)
  ggsave(paste0("02_Volcanos_highlight_",term,"_labeled.pdf"), height = 7, width = 16)
  
}

# write out lists for functional enrichment analysis in PantherDB and StringDB ( -> LT)
dir.create("02_regulatedProteinListsForFunctionalAnalyses")
setwd("02_regulatedProteinListsForFunctionalAnalyses")

# Write out lists for over/under-representation testing
up_regulated = cRes[p_value_BHadj_protein <= 0.05 & log2_fold_change_protein >= 1, unique(Protein.Group), .(compartment, comparison)]
down_regulated = cRes[p_value_BHadj_protein <= 0.05 & log2_fold_change_protein <= -1, unique(Protein.Group), .(compartment, comparison)]

for(comp in unique(cRes$compartment)){
  for (compar in unique(cRes$comparison)){
    proteins_up = up_regulated[compartment == comp & comparison == compar, unique(V1)]
    proteins_down = down_regulated[compartment == comp & comparison == compar, unique(V1)]
    
    write.table(proteins_up, file = paste(comp, gsub("\\/","vs", compar), "upregulated_grouped.txt", sep = "_"),
                sep = "\t", row.names = F, quote = F)
    write.table(proteins_down, file = paste(comp, gsub("\\/","vs", compar), "downregulated_grouped.txt", sep = "_"),
                sep = "\t", row.names = F, quote = F)
    
    proteins_up_disaggr = unique(unlist(lapply(proteins_up, function(x){strsplit(x, split = ";")})))
    proteins_down_disaggr = unique(unlist(lapply(proteins_down, function(x){strsplit(x, split = ";")})))
    
    write.table(proteins_up_disaggr, file = paste(comp, gsub("\\/","vs", compar), "upregulated_groupsResolved.txt", sep = "_"),
                sep = "\t", row.names = F, quote = F)
    write.table(proteins_down_disaggr, file = paste(comp, gsub("\\/","vs", compar), "downregulated_groupsResolved.txt", sep = "_"),
                sep = "\t", row.names = F, quote = F)
  }
}

# reference lists of all detected proteins per compartment
all_proteins_detected = cRes[, unique(Protein.Group)]
write.table(all_proteins_detected,
            file = "all_proteins_detected_grouped.txt",
            sep = "\t", row.names = F, quote = F)
all_proteins_detected = cRes[, unique(Protein.Group)]
write.table(unique(unlist(lapply(all_proteins_detected, function(x){strsplit(x, split = ";")}))),
            file = "all_proteins_detected_groupsResolved.txt",
            sep = "\t", row.names = F, quote = F)

write.table(cRes[compartment == "Lung", unique(Protein.Group)],
            file = "Lung_all_proteins_detected_grouped.txt",
            sep = "\t", row.names = F, quote = F)
all_proteins_detected = cRes[, unique(Protein.Group)]
write.table(unique(unlist(lapply(cRes[compartment == "Lung", unique(Protein.Group)], function(x){strsplit(x, split = ";")}))),
            file = "Lung_all_proteins_detected_groupsResolved.txt",
            sep = "\t", row.names = F, quote = F)

write.table(cRes[compartment == "BALF", unique(Protein.Group)],
            file = "BALF_all_proteins_detected_grouped.txt",
            sep = "\t", row.names = F, quote = F)
all_proteins_detected = cRes[, unique(Protein.Group)]
write.table(unique(unlist(lapply(cRes[compartment == "BALF", unique(Protein.Group)], function(x){strsplit(x, split = ";")}))),
            file = "BALF_all_proteins_detected_groupsResolved.txt",
            sep = "\t", row.names = F, quote = F)

setwd("../")

# Done