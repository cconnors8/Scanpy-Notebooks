#This is the code for running nichenet on a pair of csvs, on receptor and one ligand
library(nichenetr)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggpubr)

# Set the input CSV file paths for DE genes and O-link proteomics hits
de_genes_csv <- "scores-SCVI-Lifespan.csv"  # Change this to your DE genes CSV file path
olink_hits_csv <- "Media_full_ttest_results.csv"  # Change this to your O-link proteomics hits CSV file path

# Set the number of genes and proteins to use based on significance threshold
n_de_genes <- 738  # Adjust as needed
n_olink_hits <- 961  # Adjust as needed

#AS A NBOTE< YOU ONLY USE THE HEAD ON ONE OF THESE (the one that is the receptor)

# Read the DE genes and O-link hits from CSV files
de_genes <- read_csv(de_genes_csv)$long_lived_n  # Assuming the DE genes are in a column named "gene"
de_genes <- head(de_genes, n_de_genes)  # Limit to the top n DE genes

olink_hits <- read_csv(olink_hits_csv)$Assay  # Assuming O-link hits are in a column named "protein"
#olink_hits <- head(olink_hits, n_olink_hits)  # Limit to the top n proteins

# Load NicheNet data: ligand-target matrix, ligand-receptor network, and weighted networks
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network <- readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))

expressed_genes_sender <- olink_hits
expressed_genes_receiver <- de_genes


ligands <- lr_network %>% pull(from) %>% unique()
expressed_ligands <- intersect(ligands,expressed_genes_sender)

receptors <- lr_network %>% pull(to) %>% unique()
expressed_receptors <- intersect(receptors,expressed_genes_receiver)

potential_ligands <-  lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% 
  pull(from) %>% unique()


geneset_oi <-  head(de_genes, n_de_genes)  %>% .[. %in% rownames(ligand_target_matrix)] 

length(geneset_oi)



background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)

ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

(ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>%
    mutate(rank = rank(desc(aupr_corrected))))

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>%
  arrange(-aupr_corrected) %>% pull(test_ligand)

best_upstream_ligands
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 200) %>% bind_rows()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = .25)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target_network <- make_heatmap_ggplot(vis_ligand_target, "ASC-Ligands", "OLink-Genes",
                                               color = "purple") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_ligand_target_network


ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Olink-ligands", x_name = "Receptors expressed by ASC",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))


vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)
p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized O-link-ligands", "Ligand activity",
                                     color = "darkorange", legend_title = "AUPR") + 
  theme(axis.text.x.top = element_blank())
p_ligand_aupr
