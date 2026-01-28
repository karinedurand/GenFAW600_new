#!/usr/bin/env Rscript
#install.packages("BiocManager", 
#                # lib = "/storage/simple/users/durandk/R/x86_64-pc-linux-gnu-library/4.3/",
#                 dependencies = TRUE, repos = "https://cloud.r-project.org")  

# Installer LEA via BiocManager (paquet Bioconductor)
#BiocManager::install("LEA", 
#                   #  lib = "/storage/simple/users/durandk/R/x86_64-pc-linux-gnu-library/4.3/", 
#                     update = TRUE, 
#                     ask = FALSE, repos = "https://cloud.r-project.org")




################################################################################################################################
#PLOT
###############################################################################################################################

# Charger le projet SNMF
project <- load.snmfProject("WholeGenome_biallelic_max5missing_pruned.snmfProject")



# Définir la plage de K à tester
Ks <- 2:7

# Calculer la cross-entropy moyenne pour chaque K
ce <- sapply(Ks, function(k) mean(cross.entropy(project, K = k)))

# Identifier le meilleur K
best_K <- Ks[which.min(ce)]
#cat("Meilleur K :", best_K, "\n")

# ===============================
# Sauvegarder la courbe cross-entropy
# ===============================
pdf("cross_entropy.pdf")
plot(Ks, ce, type = "b", pch = 19, col = "blue",
     xlab = "K (nombre de clusters)", ylab = "Cross-Entropy",
     main = "Critère d'entropie croisée (SNMF)")
dev.off()

# Sélectionner le meilleur run pour ce K
best_run <- which.min(cross.entropy(project, K = best_K))

metadata  <- read_csv("merged_data_pop_sNMF.csv",          col_names = TRUE)


##############################################################################################
# NEW EXPORT MATRIX - V
##############################################################################################

# Sélection du meilleur run pour K=5 (ou best_K si tu l'as déjà déterminé)
best_K <- 5
best_run <- which.min(cross.entropy(project, K = best_K))

# Extraction de la Q-matrix pour le meilleur run
Q_matrix_K5 <- Q(project, K = best_K, run = best_run)

# Vérification cruciale : même nombre de lignes (et même ordre !)
stopifnot(nrow(Q_matrix_K5) == nrow(metadata))
cat("Ordre vérifié : ", nrow(Q_matrix_K5), "individus dans Q-matrix et metadata\n")

# Conversion en data.frame avec noms de colonnes clairs
Q_df <- as.data.frame(Q_matrix_K5)
colnames(Q_df) <- paste0("Cluster_", 1:best_K)

# # Ajout des colonnes utiles de metadata (sans risque de désordre)
# # Adaptez les noms de colonnes selon votre fichier metadata
# Q_df$ID <- metadata$ID       # ou metadata$SampleID, metadata$ID, etc.
# Q_df$Region     <- metadata$REGION
# Q_df$GEO_LOC_NAME <- metadata$GEO_LOC_NAME    # si cette colonne existe
# # Ajoutez d'autres colonnes si nécessaire : Q_df$Sex <- metadata$SEX, etc.

# # Réorganiser les colonnes pour plus de lisibilité : métadonnées d'abord
# desired_order <- c("ID", "Region", "GEO_LOC_NAME", 
#                    paste0("Cluster_", 1:best_K))
# # Garder seulement les colonnes qui existent
# desired_order <- desired_order[desired_order %in% colnames(Q_df)]

# Q_df <- Q_df[, desired_order]

# # Sauvegarde en CSV avec un nom explicite
# csv_filename <- paste0("Q_matrix_best_K", best_K, "_with_metadata.csv")
# write.csv(Q_df, file = csv_filename, row.names = FALSE, quote = TRUE)

# cat("Q-matrix avec métadonnées sauvegardée avec succès dans :\n")
# cat("   →", csv_filename, "\n")
# cat("   → Nombre d'individus :", nrow(Q_df), "\n")
# cat("   → Colonnes :", paste(colnames(Q_df), collapse = ", "), "\n")
##############################################################################################
##POP ordered plot 
#################################

# Ajouter les assignations Q_matrix dans les métadonnées
metadata$index <- 1:nrow(metadata)  # garder la position originale
metadata$plot_order <- factor(metadata$plot_order)  # s'assurer que pop est bien un facteur

# Réordonner les échantillons par population
order_idx <- order(metadata$plot_order)  # indices triés par pop
Q_matrix_ordered <- Q_matrix_K5[order_idx, ]
metadata_ordered <- metadata[order_idx, ]
write.csv(Q_matrix_ordered ,"Q_matrix_ordered_newpop.csv")
write.csv(metadata_ordered,"metadata_ordered_newpop.csv")
# Récupérer les labels à afficher sous le graphique (ici la pop)
labels <- metadata_ordered$GEO_LOC_NAME...17
Order_plot <- metadata_ordered$plot_order
# Sauvegarde du barplot en PDF
pdf("admixture_barplot_ordered_country.pdf", width = 7.1, height = 6.7, pointsize = 10)
barplot(t(Q_matrix_ordered),
        col = rainbow(5),
       border = NA,
        space = 0,
        main = paste("Structure génétique - K =", best_K),
        names.arg = labels,    # labels de population
        las = 2,               # labels verticaux
        cex.names = 0.7)       # taille des noms
dev.off()


##############################################################################################
## Barplots pour tous les K avec le même ordre
##############################################################################################

pdf("admixture_barplots_allK_ordered.pdf", width = 10, height = 5)

for (k in Ks) {
  # Sélectionner le meilleur run pour ce K
  best_run <- which.min(cross.entropy(project, K = k))
  
  # Extraire la Q-matrix
  Q_matrix <- Q(project, K = k, run = best_run)
  
  # Réordonner selon le même ordre d'individus
  Q_matrix_ordered <- Q_matrix[order_idx, ]
  
  # Barplot ordonné
  barplot(t(Q_matrix_ordered),
          col = rainbow(k),
          border = NA,
          space = 0,
          xlab = "Individus",
          ylab = "Proportions d'ancêtres",
          main = paste("Structure génétique - K =", k),
          names.arg = labels,
          las = 2,
          cex.names = 0.5)
}

dev.off()
##############################################################################################
## Sauvegarder les Q-matrix ordonnées pour chaque K selon l’ordre des REGIONS
##############################################################################################



# Boucle sur tous les K
for (k in Ks) {
  # Sélectionner le meilleur run pour ce K
  best_run <- which.min(cross.entropy(project, K = k))
  
  # Extraire la Q-matrix correspondante
  Q_matrix <- Q(project, K = k, run = best_run)
  
 
  # Sauvegarder la Q-matrix ordonnée
  file_name <- paste0("Q_matrix_K", k, ".csv")
  write.csv(Q_matrix_ordered, file_name, row.names = FALSE)
  
  cat("✅ Q-matrix ordonnée sauvegardée pour K =", k, "→", file_name, "\n")
}

