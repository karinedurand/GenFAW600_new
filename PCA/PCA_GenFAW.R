#===============================
#  LIBRAIRIES
#===============================
library(ggplot2)
library(tidyverse)
library(readr)
library(ggrepel)   # utilisé UNIQUEMENT pour le Whole genome

#===============================
#  PARAMÈTRES GÉNÉRAUX
#===============================
  setwd("/home/karine/Documents/Obsidian_Vault/Spodo/GenFAW_122025/PCA/")

# Palette couleur par pays
country_colors <- c(
  "Argentina"     = "#ffbb78",
  "Benin"         = "#2ca02c",
  "Brazil"        = "#ff7f0e",
  "China"         = "#e377c2",
  "French_Guiana" = "#17becf",
  "Ghana"         = "green",
  "Guadeloupe"    = "lightblue",
  "India"         = "#bcbd22",
  "Kenya"         = "#d62728",
  "Malawi"        = "purple",
  "Malaysia"      = "#8c564b",
  "Mexico"        = "black",
  "Puerto_Rico"   = "blue",
  "Rwanda"        = "red",
  "Sudan"         = "#7f7f7f",
  "USA"           = "#1f77b4",
  "Zambia"        = "#9467bd"
)

# Formes pour HOST
host_shapes <- c(21, 22, 4, 24, 25)

#===============================
#  METADATA
#===============================
  metadata <- read_delim(
    "~/Documents/Obsidian_Vault/Spodo/GenFAW_122025/LROH/metadata_whole_11122025.csv",
    delim = ",",
    trim_ws = TRUE
  )

metadata$HOST <- factor(
  metadata$HOST,
  levels = c("corn", "grasses", "missing", "rice", "sugarcane")
)

#=========================================================
#  WHOLE GENOME PCA (AVEC ggrepel)
#=========================================================
v_whole <- read.table("WholeGenome_biallelic_max80missing_pruned_PCA.eigenvec")
p_whole <- read.table("WholeGenome_biallelic_max80missing_pruned_PCA.eigenval")

vs_whole <- cbind(metadata, v_whole)
write.csv(vs_whole, "PCA_Whole.csv", row.names = FALSE)

# Variance expliquée
pc1_var <- round(p_whole$V1[1] * 100 / sum(p_whole$V1), 2)
pc2_var <- round(p_whole$V1[2] * 100 / sum(p_whole$V1), 2)

#📌 Plot PC1–PC2 (Whole genome, avec labels)
PC1PC2_whole <- ggplot(
  vs_whole,
  aes(x = V3, y = V4, color = GEO_LOC_NAME, shape = HOST)
) +
  geom_point(size = 2, stroke = 0.5) +
  scale_shape_manual(values = host_shapes) +
  scale_color_manual(values = country_colors) +
  xlab(paste0("PC1 (", pc1_var, "%)")) +
  ylab(paste0("PC2 (", pc2_var, "%)")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = c(0.8, 0.05),
    legend.justification = c(0, 0),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  
  ## 🔹 Labels des populations (ggrepel UNIQUEMENT ici)
  ggrepel::geom_text_repel(
    data = aggregate(cbind(V3, V4) ~ GEO_LOC_NAME, vs_whole, median),
    aes(x = V3, y = V4, label = GEO_LOC_NAME),
    size = 3.5,
    color = "black",
    box.padding = 0.3,
    max.overlaps = Inf,
    inherit.aes = FALSE
  )
PC1PC2_whole 
pdf("Whole_PC1_PC2_cluster.pdf", width = 10, height = 8)
PC1PC2_whole
dev.off()
#=========================================================
#  WHOLE GENOME PCA 
#=========================================================
v_whole <- read.table("WholeGenome_biallelic_max80missing_pruned_PCA.eigenvec")
p_whole <- read.table("WholeGenome_biallelic_max80missing_pruned_PCA.eigenval")

vs_whole <- cbind(metadata, v_whole)
write.csv(vs_whole, "PCA_Whole.csv", row.names = FALSE)

pc1_var_whole <- round(p_whole$V1[1] * 100 / sum(p_whole$V1), 2)
pc2_var_whole <- round(p_whole$V1[2] * 100 / sum(p_whole$V1), 2)

PC1PC2_whole <- ggplot(
  vs_whole,
  aes(x = V3, y = V4, color = GEO_LOC_NAME, shape = HOST)
) +
  geom_point(size = 2, stroke = 0.5) +
  scale_shape_manual(values = host_shapes) +
  scale_color_manual(values = country_colors) +
  xlab(paste0("PC1 (", pc1_var_whole, "%)")) +
  ylab(paste0("PC2 (", pc2_var_whole, "%)")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = c(0.82, 0.1),
    legend.justification = c(0, 0),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )
pdf("Whole_PC1_PC2.pdf", width = 10, height = 8)
PC1PC2_whole
dev.off()






#=========================================================
#  AUTOSOME PCA 
#=========================================================
v_auto <- read.table("Autosome_biallelic_max80missing_pruned_PCA.eigenvec")
p_auto <- read.table("Autosome_biallelic_max80missing_pruned_PCA.eigenval")

vs_auto <- cbind(metadata, v_auto)
write.csv(vs_auto, "PCA_Autosome.csv", row.names = FALSE)

pc1_var_auto <- round(p_auto$V1[1] * 100 / sum(p_auto$V1), 2)
pc2_var_auto <- round(p_auto$V1[2] * 100 / sum(p_auto$V1), 2)

PC1PC2_auto <- ggplot(
  vs_auto,
  aes(x = V3, y = V4, color = GEO_LOC_NAME, shape = HOST)
) +
  geom_point(size = 2, stroke = 0.5) +
  scale_shape_manual(values = host_shapes) +
  scale_color_manual(values = country_colors) +
  xlab(paste0("PC1 (", pc1_var_auto, "%)")) +
  ylab(paste0("PC2 (", pc2_var_auto, "%)")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = c(0.82, 0.1),
    legend.justification = c(0, 0),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )
pdf("Autosome_PC1_PC2.pdf", width = 10, height = 8)
PC1PC2_auto
dev.off()


#=========================================================
#  Z PCA 
#=========================================================
v_Z<- read.table("Z_biallelic_max80missing_pruned_PCA.eigenvec")
p_Z<- read.table("Z_biallelic_max80missing_pruned_PCA.eigenval")

vs_Z<- cbind(metadata, v_auto)
write.csv(vs_auto, "PCA_Z.csv", row.names = FALSE)

pc1_var_Z<- round(p_Z$V1[1] * 100 / sum(p_Z$V1), 2)
pc2_var_Z<- round(p_Z$V1[2] * 100 / sum(p_Z$V1), 2)

#📌 Plot PC1–PC2 (Autosomes, sans labels)
PC1PC2_Z<- ggplot(
  vs_Z,
  aes(x = V3, y = V4, color = GEO_LOC_NAME, shape = HOST)
) +
  geom_point(size = 2, stroke = 0.5) +
  scale_shape_manual(values = host_shapes) +
  scale_color_manual(values = country_colors) +
  xlab(paste0("PC1 (", pc1_var_Z, "%)")) +
  ylab(paste0("PC2 (", pc2_var_Z, "%)")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = c(0.823, 0.1),
    legend.justification = c(0, 0),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

pdf("Z_PC1_PC2.pdf", width = 10, height = 8)
PC1PC2_Z
dev.off()

