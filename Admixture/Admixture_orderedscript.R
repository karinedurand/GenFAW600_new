library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

#https://connor-french.github.io/intro-pop-structure-r/
setwd("/home/karine/Documents/Obsidian_Vault/Spodo/GenFAW_122025/Admixture/Admixture_sansmaf_geno_dossierok_meme_fichier_treemix_pca/")
width_a4 <- 9      # Largeur A4 en pouces (portrait)
height_a4_5 <- 9 / 5  # 1/5 de la hauteur A4 en pouces

# Read cross-validation data
ce <- read.csv("CV.csv", sep="\t",header=FALSE )

# Plot cross-validation data
f.ce <- ggplot(ce , aes(x = V4, y = V5)) + geom_point() + 
  xlab("K") + ylab("CE") + theme_bw()
f.ce

# Save cross-validation plot to PDF
pdf("CV.pdf", width = 5, height = 5)
f.ce
dev.off()
pops  <- read.csv("/home/karine/Documents/Obsidian_Vault/Spodo/GenFAW_122025//Admixture/metadata_whole_11122025_pop.csv", sep=","  )

###Now time for plotting! First, you need the Q-matrix, which contains the ancestry proportions. You’re extracting the Q matrix for the best run of K=3. Each column is a population and each row is an individual, so let’s add column names to make that more apparent.
q_mat <- read_table("admix_geno0.05.2.Q",   col_names = FALSE)
new_pop <- read.csv("new_pops_sNMF.csv", sep="," )
merged_data <- merge(pops, new_pop, by = "ID", all = TRUE)
#write.csv( merged_data,"newpop_metadata.csv")
merged_data_popsNMF <- merged_data %>% 
  arrange(order)   # ordre croissant par défaut
write.csv(merged_data,"merged_data.csv")

colnames(q_mat) <- paste0("P", 1:6)

head(q_mat)

q_df <- q_mat %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(individual = pops$ID,
         region = pops$X1,
         order = pops$order)
q_df
#write.csv(q_df, "q_dfQ2.csv", row.names=TRUE)
q_df_long <- q_df %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q") 

q_df_long
q_df_2 <- q_df_long %>% 
  # arrange the data set by the plot order indicated in Prates et al.
  arrange(order) %>% 
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  mutate(individual = forcats::fct_inorder(factor(individual)))


q_palette <- c("#fde725", "#35b779", "#440154","red","orange")


k2<- q_df_2 %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  scale_fill_manual(values = q_palette) +  # Appliquer la palette personnalisée
  theme_minimal() +
  # Ajustements esthétiques
  theme(
    panel.spacing.x = unit(0, "lines"),
    axis.line = element_blank(),
    #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    axis.text.x = element_blank(),  # <--- Supprime les étiquettes de l'axe x,
    strip.background = element_rect(fill = "transparent", color = "black"),
    panel.background = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )
k2
width_a4 <- 9      # Largeur A4 en pouces (portrait)
height_a4_5 <- 9 / 5  # 1/5 de la hauteur A4 en pouces
ggsave(
  filename = "K2.pdf",  # Nom du fichier PDF
  plot = k2,                  # Objet ggplot à sauvegarder
  width = width_a4,                # Largeur A4
  height = height_a4_5,            # Hauteur 1/5 de A4
  units = "in",                    # Unités en pouces
  device = "pdf"                   # Format de sortie : PDF
)


####################################################################################################################""
q_mat <- read_table("admix_geno0.05.3.Q",   col_names = FALSE)

colnames(q_mat) <- paste0("P", 1:5)
head(q_mat)

q_df <- q_mat %>% 
  as_tibble() %>% 
  mutate(individual = pops$ID,
         region = pops$X1,
         order = pops$order)
q_df
#write.csv(q_df, "q_dfQ2.csv", row.names=TRUE)
q_df_long <- q_df %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q") 

q_df_long
q_df_3 <- q_df_long %>% 
  # arrange the data set by the plot order indicated in Prates et al.
  arrange(order) %>% 
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  mutate(individual = forcats::fct_inorder(factor(individual)))


q_palette <- c("#fde725", "#35b779", "#440154","red","orange")


k3 <- q_df_3 %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  scale_fill_manual(values = q_palette) +  # Appliquer la palette personnalisée
  theme_minimal() +
  # Ajustements esthétiques
  theme(
    panel.spacing.x = unit(0, "lines"),
    axis.line = element_blank(),
    #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    axis.text.x = element_blank(),  # <--- Supprime les étiquettes de l'axe x,
    strip.background = element_rect(fill = "transparent", color = "black"),
    panel.background = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )


k3
ggsave(
  filename = "K3.pdf",  # Nom du fichier PDF
  plot = k3,                  # Objet ggplot à sauvegarder
  width = width_a4,                # Largeur A4
  height = height_a4_5,            # Hauteur 1/5 de A4
  units = "in",                    # Unités en pouces
  device = "pdf"                   # Format de sortie : PDF
)
######################################################################################################################"
q_mat <- read_table("admix_geno0.05.4.Q",   col_names = FALSE)

colnames(q_mat) <- paste0("P", 1:5)

q_df <- q_mat %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(individual = merged_data$ID,
         region = merged_data$new_pop,
         order = pops$order)
q_df


q_palette <- c("#fde725", "#35b779", "#440154","red","orange")

#----------------------------------------------
ordre_regions <- c("Nat_C", "Nat_R","Mex", "Zam", "Inv_west","Inv_east")  # ← complète et ordonne comme tu veux

q_df_long <- q_df_long %>%
  mutate(
    region = factor(region, levels = ordre_regions),           # ordre des groupes
    pop    = factor(pop, levels = paste0("P", 1:6)),           # P1 en bas, P6 en haut
    individual = forcats::fct_inorder(individual)              # ordre interne conservé
  ) %>%
  arrange(region, order) %>%                                   # tri par région, puis par ton order interne
  mutate(individual = forcats::fct_inorder(individual))        # fixe l'ordre pour ggplot
q_palette <- c("red", "#fde725", "orange", "#35b779", "blue", "#440154")
k4 <- ggplot(q_df_long, aes(x = individual, y = q, fill = pop)) +
  geom_col(width = 1) +
  scale_fill_manual(values = q_palette) +
  facet_grid(~ region, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    panel.grid = element_blank(),
    strip.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none",          # ← Supprime complètement la légende
    panel.spacing.x = unit(0.3, "lines")
  ) +
  ylab("Ancestry proportion") +        # ← Axe y en anglais
  labs(title = NULL)                   # ← Au cas où il y aurait un titre, on l'enlève

# Afficher le graphique
print(k4)
# Sauvegarder (recommandé pour bien voir le résultat)
ggsave(
  filename = "new_pops_K6.pdf",  # Nom du fichier PDF
  plot = k6,                  # Objet ggplot à sauvegarder
  width = width_a4,                # Largeur A4
  height = height_a4_5,            # Hauteur 1/5 de A4
  units = "in",                    # Unités en pouces
  device = "pdf"                   # Format de sortie : PDF
)
#--------------------------------------------
#write.csv(q_df, "q_dfQ2.csv", row.names=TRUE)
q_df_long <- q_df %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q") 

q_df_long
q_df_4 <- q_df_long %>% 
  # arrange the data set by the plot order indicated in Prates et al.
  arrange(order) %>% 
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  mutate(individual = forcats::fct_inorder(factor(individual)))

k4<-q_df_4 %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  scale_fill_manual(values = q_palette) +  # Appliquer la palette personnalisée
  theme_minimal() +
  # Ajustements esthétiques
  theme(
    panel.spacing.x = unit(0, "lines"),
    axis.line = element_blank(),
    #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    axis.text.x = element_blank(),  # <--- Supprime les étiquettes de l'axe x,
    strip.background = element_rect(fill = "transparent", color = "black"),
    panel.background = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )
k4

ggsave(
  filename = "K4.pdf",  # Nom du fichier PDF
  plot = k4,                  # Objet ggplot à sauvegarder
  width = width_a4,                # Largeur A4
  height = height_a4_5,            # Hauteur 1/5 de A4
  units = "in",                    # Unités en pouces
  device = "pdf"                   # Format de sortie : PDF
)
##########################################################################################################################

q_mat <- read_table("admix_geno0.05.5.Q",   col_names = FALSE)

colnames(q_mat) <- paste0("P", 1:5)
q_df <- q_mat %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(individual = pops$ID,
         region = pops$X1,
         order = pops$order)
q_df

#write.csv(q_df, "q_dfQ2.csv", row.names=TRUE)
q_df_long <- q_df %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q") 

q_df_long
q_df_5 <- q_df_long %>% 
  # arrange the data set by the plot order indicated in Prates et al.
  arrange(order) %>% 
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  mutate(individual = forcats::fct_inorder(factor(individual)))


q_palette <- c("#fde725", "#35b779", "#440154","red","orange")


# Sauvegarde du barplot en PDF

k5<-q_df_5 %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  scale_fill_manual(values = q_palette) +  # Appliquer la palette personnalisée
  theme_minimal() +
  # Ajustements esthétiques
  theme(
    panel.spacing.x = unit(0, "lines"),
    axis.line = element_blank(),
    #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    axis.text.x = element_blank(),  # <--- Supprime les étiquettes de l'axe x,
    strip.background = element_rect(fill = "transparent", color = "black"),
    panel.background = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )
k5
ggsave(
  filename = "K5.pdf",  # Nom du fichier PDF
  plot = k5,                  # Objet ggplot à sauvegarder
  width = width_a4,                # Largeur A4
  height = height_a4_5,            # Hauteur 1/5 de A4
  units = "in",                    # Unités en pouces
  device = "pdf"                   # Format de sortie : PDF
)
###########################################################################################################"
q_mat <- read_table("admix_geno0.05.6.Q",   col_names = FALSE)

colnames(q_mat) <- paste0("P", 1:6)
q_df <- q_mat %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(ID = merged_data$ID,
         region = merged_data$new_pop,
         order = pops$order)
q_df
write.csv(q_df ,"K6.csv ")
K6_metadata<- merge(q_df,merged_data , by = "ID", all = TRUE)
write.csv(K6_metadata, "K6_withmetadata.csv", row.names=TRUE)
q_df_long <- q_df %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q") 

q_df_long

#-----------
ordre_regions <- c("Nat_C", "Nat_R","Mex", "Zam", "Inv_west","Inv_east")  # ← complète et ordonne comme tu veux

q_df_long <- q_df_long %>%
  mutate(
    region = factor(region, levels = ordre_regions),           # ordre des groupes
    pop    = factor(pop, levels = paste0("P", 1:6)),           # P1 en bas, P6 en haut
    ID = forcats::fct_inorder(ID)              # ordre interne conservé
  ) %>%
  arrange(region, order) %>%                                   # tri par région, puis par ton order interne
  mutate(ID = forcats::fct_inorder(ID))        # fixe l'ordre pour ggplot
q_palette <- c("red", "#fde725", "orange", "#35b779", "blue", "#440154")
k6 <- ggplot(q_df_long, aes(x = ID, y = q, fill = pop)) +
  geom_col(width = 1) +
  scale_fill_manual(values = q_palette) +
  facet_grid(~ region, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    panel.grid = element_blank(),
    strip.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none",          # ← Supprime complètement la légende
    panel.spacing.x = unit(0.3, "lines")
  ) +
  ylab("Ancestry proportion") +        # ← Axe y en anglais
  labs(title = NULL)                   # ← Au cas où il y aurait un titre, on l'enlève

# Afficher le graphique
print(k6)

# Sauvegarder (recommandé pour bien voir le résultat)
ggsave(
  filename = "K6.pdf",  # Nom du fichier PDF
  plot = k6,                  # Objet ggplot à sauvegarder
  width = width_a4,                # Largeur A4
  height = height_a4_5,            # Hauteur 1/5 de A4
  units = "in",                    # Unités en pouces
  device = "pdf"                   # Format de sortie : PDF
)






q_df_6 <- q_df_long %>% 
  # arrange the data set by the plot order indicated in Prates et al.
  arrange(order) %>% 
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  mutate(ID = forcats::fct_inorder(factor(ID)))


q_palette <- c("red","#fde725", "orange","#35b779", "blue","#440154")


k6<-q_df_6 %>% 
  ggplot() +
  geom_col(aes(x = indivu, y = q, fill = pop)) +
  scale_fill_manual(values = q_palette) +  # Appliquer la palette personnalisée
  theme_minimal() +
  # Ajustements esthétiques
  theme(
    panel.spacing.x = unit(0, "lines"),
    axis.line = element_blank(),
    #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    axis.text.x = element_blank(),  # <--- Supprime les étiquettes de l'axe x,
    strip.background = element_rect(fill = "transparent", color = "black"),
    panel.background = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )
k6
ggsave(
  filename = "K6.pdf",  # Nom du fichier PDF
  plot = k6,                  # Objet ggplot à sauvegarder
  width = width_a4,                # Largeur A4
  height = height_a4_5,            # Hauteur 1/5 de A4
  units = "in",                    # Unités en pouces
  device = "pdf"                   # Format de sortie : PDF
)


library(gridExtra)
pdf("admixturek6.pdf")

common_theme <- theme(
  legend.position = "none",
  axis.text.x  = element_blank(),
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank()
)

grid.arrange(
  k2 + common_theme,
  k3 + common_theme,
  k4 + common_theme,
  k5 + common_theme,
  k6 + common_theme,
  ncol = 1
)

dev.off()
pdf("admixture_ALLK.pdf")
par(mfrow = c(3, 1), mar = c(2, 4, 2, 1), oma = c(0, 0, 2, 0))
k2 + theme(legend.position = "none")
k3 + theme(legend.position = "none") 
k4 + theme(legend.position = "none") 
k5 + theme(legend.position = "none") 
dev.off()
pdf("admixture_ALLK.pdf")
grid.arrange(k2 + theme(legend.position = "none"),
             k3 + theme(legend.position = "none"),
             k4 + theme(legend.position = "none"),
             k5 + theme(legend.position = "none") ,
             heights=c(1,1,1,1))
dev.off()


#########################################################################################################################################"
