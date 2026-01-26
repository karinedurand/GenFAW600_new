#!/bin/R

setwd("/home/karine/Documents/Obsidian_Vault/Spodo/GenFAW_122025/TREEMIX/janvier_seed/")
source("plotting_funcs.R")
prefix="treemix.m"
library(RColorBrewer)
library(R.utils)
pdf("Besttreem1.pdf",width=15, height=7)
par(mfrow=c(2,3))
for(edge in 3:3){
  plot_tree(cex=0.8,paste0(prefix,".",edge))
  title(paste(edge,"edges"))
}
dev.off()


source("plotting_funcs.R")

# Remplace cette ligne par ton stem exact (sans .treeout.gz ni rien d'autre)
stem <- "treemix.m0.rep1"   # ← À ADAPTER : ton meilleur m et replicate

# Plot et sauvegarde en PDF
pdf("Best_TreeMix_m0.pdf", width = 15, height = 8)

plot_tree(stem, cex = 0.8)   # Affiche l’arbre avec la(les) migration(s)

dev.off()

# Chargement des fonctions TreeMix
source("plotting_funcs.R")

# Ouverture du PDF (une page avec 7 arbres en 1 ligne)
pdf("TreeMix_m0_to_m5_rep1.pdf", width = 15, height = 7)  # Large pour 7 arbres côte à côte

# Configuration : 3 ligne, 2 colonnes
par(mfrow=c(2,3))

# Boucle pour m = 0 à 5, avec replicate 1
for (m in 0:5) {
  stem <- paste0("treemix.m", m, ".rep1")  # Stem exact : treemix.m0.rep1, treemix.m1.rep1, etc.
  
  plot_tree(stem, cex = 0.8)  # Plot l'arbre (lit automatiquement .treeout.gz)
  
  title(main = paste(m, "edge", if(m > 1) "s" else ""), 
        cex.main = 1.8, col.main = "black")
}

dev.off()

cat("Tous les arbres (m=0 à m=6, replicate 1) sauvegardés dans : TreeMix_m0_to_m6_rep1.pdf\n")

# Message de confirmation
cat("Arbre sauvegardé dans Best_TreeMix_final.pdf\n")




LL_ <- read_delim("LL.csv", delim = ",",  escape_double = FALSE, col_names = FALSE,    trim_ws = TRUE)
# Plot cross-validation data
LL <- ggplot(LL_ , aes(x = X1, y = X2)) + geom_point() + 
  xlab("m") + ylab("Likelihood") + theme_bw()
LL

# Save cross-validation plot to PDF
pdf("LL.pdf", width = 5, height = 5)
LL
dev.off()

# Stem du meilleur run
stem_best <- paste0("treemix_optm_runs/treemix.m1.rep", best_rep_m1)

# Plot l’arbre
plot_tree(stem_best, plot_mig = TRUE)

# Sauvegarde en PDF
pdf("TreeMix_m1_best.pdf", width = 12, height = 8)
plot_tree(stem_best)
dev.off()

# Plot les residuals (pour vérifier le fit)
pdf("Residuals_m1_best.pdf")
plot_resid(stem_best)
dev.off()