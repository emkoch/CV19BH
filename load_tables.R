library(tidyverse)
library(ggplot2)
library(ggsci)
library(lubridate)
library(corrplot)
library(data.table)
library(cowplot)
library(scales)
library(aod)
library(gaston)

source("genotype_functions.R")

problematic.sites <- read.table("problematic_sites_sarsCov2.vcf")

data_muts <- read_csv("seq_data_hq_with_mutations_final.csv")

genes <- c("nsp1", "nsp2", "nsp3", "nsp4", "3CL-PRO", "nsp6", "nsp7", "nsp8", "nsp9", "nsp10",
           "Pol", "Hel", "ExoN", "nsp15", "nsp16", "S glycoprotein", "ORF3a", "E", "M", "ORF6", "ORF7a",
           "ORF8", "N")
idNdS <- c(0.51, 0.79, 0.45, 0.44, 2.08, 0.58, 0.14, 0.32, 0.46, 1.99, 0.34, 0.37, 0.48, 0.57, 0.54, 0.6, 2.29, 0.15,
           0.51, 0.97, 1.43, 0.17, 0.81)
gene.lengths <- c(180, 638, 1945, 500, 306, 290, 83, 198, 113, 139, 932, 601, 527, 346, 298, 1273, 275, 75, 222, 61, 121, 121, 419)

col.names <- do.call(c, lapply(genes, function(gene.name){
  colnames(data_muts)[grepl(paste0("^", gene.name, "_"), colnames(data_muts))]
}))

col.names.kill <- col.names[colSums(data_muts[,col.names]) == 0]
data_muts <- data_muts[, ! colnames(data_muts) %in% col.names.kill]

mut_annotations <- read_csv("mutation_annotations_082721.csv")
sum(duplicated(mut_annotations$mutation))

mut_annos <- data.table::data.table(mut_annotations)
data.table::setkey(mut_annos, mutation)

mut_annos %>% dplyr::mutate(mut_name = paste(gene, mutation, sep="_")) %>% dplyr::filter(mut_name %in% colnames(data_muts)) -> mut_annos
mut_annos %>% dplyr::mutate(gene_position = paste(gene, `aa index`, sep=":")) %>% dplyr::filter(!is.na(gene)) -> mut_annos
position_counts <- table(mut_annos$gene_position)
multi_pos <- names(position_counts)[position_counts > 1]
multi_muts <- dplyr::filter(mut_annos, gene_position %in% multi_pos)

position_counts <- list()
for(position in unique(multi_muts$gene_position)){
  gene <- stringr::str_split(position, ":")[[1]][1]
  foo <- make.geno.mat(data_muts, gene, var.type = "coding", multi_muts %>% dplyr::filter(gene_position == position))
  position_counts[[position]] <- table(rowSums(foo))
}

co_occur <- sapply(position_counts, function(X) ("2" %in% names(X)) | ("3" %in% names(X)) | ("4" %in% names(X)))
co_occur_names <- names(position_counts)[co_occur]
co_occur_muts <- (multi_muts %>% dplyr::filter(gene_position %in% co_occur_names))$mut_name

data_muts_noco <- dplyr::select(data_muts, -co_occur_muts)
mut_annos_noco <- dplyr::filter(mut_annos, ! gene_position %in% co_occur_names)

keep.severe <- !is.na(data_muts_noco$Disease.Severity)
data_muts_severe <- data_muts_noco[keep.severe,]
data_muts_severe <- relabel_severe(data_muts_severe)

syn.genes <- lapply(genes, function(gene){
  filter.positions(make.geno.mat(dat = data_muts_severe, gene.names = gene, var.type = "synonymous", lookup.table = mut_annos_noco, kill.zeros = T),
                   problematic.sites$V2)
}); names(syn.genes) <- genes

ms.genes <- lapply(genes, function(gene){
  filter.positions(make.geno.mat(dat = data_muts_severe, gene.names = gene, var.type = "missense", lookup.table = mut_annos_noco, kill.zeros = T),
                   problematic.sites$V2)
}); names(ms.genes) <- genes

all.gene.ms <-  filter.positions(make.geno.mat(dat = data_muts_severe, gene.names = genes, var.type = "missense",
                              lookup.table = mut_annos_noco, kill.zeros = T), problematic.sites$V2)

all.gene.coding <- filter.positions(make.geno.mat(dat = data_muts_severe, gene.names = genes, var.type = "coding",
                                 lookup.table = mut_annos_noco, kill.zeros = T), problematic.sites$V2)

all.gene.syn <-  filter.positions(make.geno.mat(dat = data_muts_severe, gene.names = genes, var.type = "synonymous",
                               lookup.table = mut_annos_noco, kill.zeros = T), problematic.sites$V2)

pca_coding <- calc.pca(all.gene.coding)
pca_synonymous <- calc.pca(all.gene.syn)
pca_missense <- calc.pca(all.gene.ms)

data_muts_severe <- add.pcs(data_muts_severe, pca_coding, pca.id = "pc.coding.")
data_muts_severe <- add.pcs(data_muts_severe, pca_synonymous, pca.id = "pc.synonymous.")
data_muts_severe <- add.pcs(data_muts_severe, pca_missense, pca.id = "pc.missense.")

clades.other <- c("20C", "20D", "20E (EU1)", "20G", "20H (Beta, V2)", "21D (Eta)", "21E (Theta)")
data_muts_severe$clade_grouped <- data_muts_severe$clade
data_muts_severe[data_muts_severe$clade_grouped %in% clades.other,]$clade_grouped <- "Other"
data_muts_severe$clade_grouped <- factor(data_muts_severe$clade_grouped,
                                         levels=c("20A", "19A", "19B", "20B", "20I (Alpha, V1)", "Other"))

gene.colors <- ggsci::pal_igv()(length(genes)); names(gene.colors) <- genes
gene.colors.alt <- gene.colors
gene.colors.alt[gene.lengths>200] <- gene.colors[1:sum(gene.lengths>200)]
gene.colors.alt[gene.lengths<=200] <-  gene.colors[(sum(gene.lengths>200)+1):length(gene.lengths)]
