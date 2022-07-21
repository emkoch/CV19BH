library(dplyr)
# I fixed the ORF here to be in line with the AA cooridnates
# of ORF1ab given here: https://www.ncbi.nlm.nih.gov/gene/?term=YP_009725307.1
prot_coords <- read.table("unipChainCov2_altered.bed", sep="\t", header=FALSE, quote=NULL)
prot_coords <- prot_coords[-6,]
# Narrow the genome down to non-overlapping protein products
keep_prots <- c("nsp1", "nsp2", "nsp3", "nsp4", "3CL-PRO", "nsp6", "nsp7", "nsp8", "nsp9", "nsp10",
                "Pol", "Hel", "ExoN", "nsp15", "nsp16", "S glycoprotein", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b",
                "ORF8", "N")
prot_coords <- dplyr::filter(prot_coords, V4 %in% keep_prots)
ns_coords <- read.table("nextstrainGene.bed", sep="\t", header=F) %>% dplyr::filter(V4 != "ORF9b")
# Verify that all proteins within the uniprot table are in the same open reading frame as the nextstrain proteins
for(ii in 1:nrow(prot_coords)){
  # Is the first coord a start position in the nextstrain ORF where it falls
  pos1 <- prot_coords$V2[ii]
  ns_ii_1 <- max(which((ns_coords$V2 <= pos1) & (ns_coords$V3 >= pos1)))
  print((pos1 - ns_coords$V2[ns_ii_1]) %% 3)
  # Is the second coord an end position in the nextstrain ORF where it falls
  pos2 <- prot_coords$V3[ii]
  ns_ii_2 <- min(which((ns_coords$V2 <= pos2) & (ns_coords$V3 >= pos2)))
  print((pos2 - ns_coords$V2[ns_ii_2]) %% 3)
  print(c(prot_coords$V4[ii], ns_ii_1, ns_coords$V4[ns_ii_1], "-", ns_ii_2, ns_coords$V4[ns_ii_2]))
}
