source("load_tables.R")

## Calculate burden associations
## Include clades + PCs, or either alone, as covariates
## Each gene individually, then all together

ms.burden.severity <- lapply(genes, function(gene){
  print(gene)
  get.rv.sequence.severity(data_muts_severe, ms.genes[[gene]], all.gene.syn, use.clade = T)
}); names(ms.burden.severity) <- genes

ms.burden.severity.clade <- lapply(genes, function(gene){
  print(gene)
  get.rv.sequence.severity(data_muts_severe, ms.genes[[gene]], all.gene.syn, use.clade = T, use.PC = F)
}); names(ms.burden.severity.clade) <- genes

ms.burden.severity.PC <- lapply(genes, function(gene){
  print(gene)
  get.rv.sequence.severity(data_muts_severe, ms.genes[[gene]], all.gene.syn, use.clade = F, use.PC = T)
}); names(ms.burden.severity.PC) <- genes

ms.burden.all.severity <- get.rv.sequence.severity(data_muts_severe, all.gene.ms, all.gene.syn, use.clade=T)
syn.burden.all.severity <- get.rv.sequence.severity(data_muts_severe, all.gene.syn, use.clade=T)

## Plot The effects of burden scores from all missense and all synonymous

pp.burden.total.severity <- cowplot::plot_grid(
  rv.burden.plot(get.rv.sequence.severity(data_muts_severe, all.gene.ms,  use.clade=T), title = "All missense", xlimits=c(-0.5, 0.5)),
  rv.burden.plot(get.rv.sequence.severity(data_muts_severe, all.gene.syn,  use.clade=T), title = "All synonymous",xlimits=c(-0.5, 0.5) ))
pp.burden.total.severity
ggsave(plot=pp.burden.total.severity, filename="Plots/Final/Burden/syn_ms_comapre_pc_plus_clades_severity.pdf")

## Generate tables of p-values for max allele counts 10 and 100

p.table.severity <- data.frame(do.call(rbind, lapply(ms.burden.severity, FUN = function(X){
  c(tail(X$p.value, n=1),
    tail(X$p.value[X$count.max<=10], n=1))
})))
colnames(p.table.severity) <- c("100", "10")
p.table.severity$gene <- unlist(rownames(p.table.severity))
p.table.severity$gene.lengths <- gene.lengths
p.table.severity %>% tidyr::pivot_longer(cols=1:2, names_to = "max count", values_to="p-value") -> p.table.severity
p.table.severity$gene <- factor(p.table.severity$gene, levels=genes)

p.table.severity.clade <- data.frame(do.call(rbind, lapply(ms.burden.severity.clade, FUN = function(X){
  c(tail(X$p.value, n=1),
    tail(X$p.value[X$count.max<=10], n=1))
})))
colnames(p.table.severity.clade) <- c("100", "10")
p.table.severity.clade$gene <- unlist(rownames(p.table.severity.clade))
p.table.severity.clade$gene.lengths <- gene.lengths
p.table.severity.clade %>% tidyr::pivot_longer(cols=1:2, names_to = "max count", values_to="p-value") -> p.table.severity.clade
p.table.severity.clade$gene <- factor(p.table.severity.clade$gene, levels=genes)

p.table.severity.PC <- data.frame(do.call(rbind, lapply(ms.burden.severity.PC, FUN = function(X){
  c(tail(X$p.value, n=1),
    tail(X$p.value[X$count.max<=10], n=1))
})))
colnames(p.table.severity.PC) <- c("100", "10")
p.table.severity.PC$gene <- unlist(rownames(p.table.severity.PC))
p.table.severity.PC$gene.lengths <- gene.lengths
p.table.severity.PC %>% tidyr::pivot_longer(cols=1:2, names_to = "max count", values_to="p-value") -> p.table.severity.PC
p.table.severity.PC$gene <- factor(p.table.severity.PC$gene, levels=genes)

## Make p-value barplots

pp.pval.severity <- ggplot(p.table.severity %>% dplyr::filter(gene.lengths>200)) +
  geom_bar(aes(x=gene, y=-log10(`p-value`), fill=`max count`, color=gene), stat="identity", position="dodge", size=1) +
  theme_minimal_hgrid() +
  scale_color_igv() +
  scale_fill_grey() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2)) +
  ylab("-log10(p-value)")+ xlab(NULL) + guides(color=FALSE)  + 
  geom_hline(aes(yintercept=-log10(0.05)), lty="dotted") +
  geom_hline(aes(yintercept=-log10(0.05/12)))

pp.pval.severity

pp.pval.severity.clade <- ggplot(p.table.severity.clade %>% dplyr::filter(gene.lengths>200)) +
  geom_bar(aes(x=gene, y=-log10(`p-value`), fill=`max count`, color=gene), stat="identity", position="dodge", size=1) +
  theme_minimal_hgrid() +
  scale_color_igv() +
  scale_fill_grey() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2)) +
  ylab("-log10(p-value)")+ xlab(NULL) + guides(color=FALSE)  + 
  geom_hline(aes(yintercept=-log10(0.05)), lty="dotted") +
  geom_hline(aes(yintercept=-log10(0.05/12)))+ ggtitle("clade only")

pp.pval.severity.clade 

pp.pval.severity.PC <- ggplot(p.table.severity.PC %>% dplyr::filter(gene.lengths>200)) +
  geom_bar(aes(x=gene, y=-log10(`p-value`), fill=`max count`, color=gene), stat="identity", position="dodge", size=1) +
  theme_minimal_hgrid() +
  scale_color_igv() +
  scale_fill_grey() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2)) +
  ylab("-log10(p-value)")+ xlab(NULL) + guides(color=FALSE)  + 
  geom_hline(aes(yintercept=-log10(0.05)), lty="dotted") +
  geom_hline(aes(yintercept=-log10(0.05/12))) + ggtitle("PCs only")

pp.pval.severity.PC

cowplot::plot_grid(pp.pval.severity + ggtitle("Missense"), pp.pval.severity.clade, pp.pval.severity.PC, ncol=1, byrow=F)
ggsave("Plots/Final/Burden/p_value_comparison.pdf", width=6, height=10)

## Make burden score by max count plots for each gene

pp.set.ms.severity <- list()
for(ii in 1:length(genes)){
  pp <- rv.burden.plot(data.test = ms.burden.severity[[genes[ii]]], color = gene.colors.alt[ii], title = paste(genes[ii], "missense"))
  pp.set.ms.severity[[ii]] <- pp
  ggsave(filename = file.path("Plots/Final/Burden", paste0("severity_ms_", genes[ii], "_burden.pdf")), plot = pp, width=4, height=4)
}

ggsave(plot=cowplot::plot_grid(plotlist=c(pp.set.ms.severity[gene.lengths>200], 
                                          list(rv.burden.plot(ms.burden.all.severity, title = "All missense"))), ncol=4),
       filename = file.path("Plots/Final/Burden/severity_ms_grid_burden.pdf"), width = 16, height=18)

## Permute variants matching frequencies from each gene and calculate burden scores

ms.perms.severity <- lapply(genes, function(gene){
  match.rv.samps(data_muts_severe, geno.mat = ms.genes[[gene]], comparison.mat = all.gene.ms, match.mat = all.gene.syn, n.samps=1000)
})
names(ms.perms.severity) <- genes

saveRDS(ms.perms.severity, "ms_test_perms_severity.RDS")

ms.perms.severity <- readRDS("ms_test_perms_severity.RDS")
names(ms.perms.severity ) <- genes

## Make showing the distribution of association statistics from permuted burden scores

perm.plots.ms.severity <- lapply(genes, function(gene){
  print(gene)
  plot.p.perms(ms.burden.severity[[gene]], ms.perms.severity[[gene]], gene.name = paste(gene, "missense"), color = gene.colors.alt[[gene]])
}); names(perm.plots.ms.severity) <- genes

perm.plots.ms.severity.z <- lapply(genes, function(gene){
  plot.p.perms.z(ms.burden.severity[[gene]], ms.perms.severity[[gene]], gene.name = paste(gene, "missense"), color = gene.colors.alt[[gene]])
}); names(perm.plots.ms.severity.z) <- genes

ggsave(plot = cowplot::plot_grid(plotlist = perm.plots.ms.severity[gene.lengths>200], ncol=4),
       filename="Plots/Final/Burden/ms_perm_grid_burden_severity.pdf", width = 16, height=21)

ggsave(plot = cowplot::plot_grid(plotlist = perm.plots.ms.severity.z[gene.lengths>200], ncol=4),
       filename="Plots/Final/Burden/ms_perm_grid_burden_severity_z.pdf", width = 16, height=21)

for(gene in genes){
  pp <- plot.p.perms(ms.burden.severity[[gene]], ms.perms.severity[[gene]], gene.name = paste(gene, "missense"), color = gene.colors[[gene]])
  ggsave(filename = file.path("Plots/Final/Burden/", paste0("ms_", gene, "_perm_burden_severity.pdf")), plot = pp)
}

## Make percentile tables for different max count cutoffs

perc.table.10.severity <- do.call(rbind,
                             lapply(genes, function(gene){
                               result <- data.frame(
                                 percentile=get.perc.perms(ms.burden.severity[[gene]], ms.perms.severity[[gene]]),
                                 gene=gene,
                                 Outcome="Severity",
                                 count.max=ms.burden.severity[[gene]]$count.max,
                                 beta=ms.burden.severity[[gene]]$beta
                               )
                               result <- tail(dplyr::filter(result, count.max <= 10), n=1)
                               return(result)
                             })
)

perc.table.100.severity <- do.call(rbind,
                                  lapply(genes, function(gene){
                                    result <- data.frame(
                                      percentile=get.perc.perms(ms.burden.severity[[gene]], ms.perms.severity[[gene]]),
                                      gene=gene,
                                      Outcome="Severity",
                                      count.max=ms.burden.severity[[gene]]$count.max,
                                      beta=ms.burden.severity[[gene]]$beta
                                    )
                                    result <- tail(dplyr::filter(result, count.max <= 100), n=1)
                                    return(result)
                                  })
)

perc.table.100.severity$gene <- factor(perc.table.100.severity$gene, levels = genes)
perc.table.10.severity$gene <- factor(perc.table.10.severity$gene, levels = genes)

severity.perc.bars <- dplyr::bind_rows(dplyr::mutate(perc.table.10.severity, `max count`="10"),
                 dplyr::mutate(perc.table.100.severity, `max count`="100")) %>% 
  dplyr::mutate(gene.lengths=rep(gene.lengths, 2)) %>%
  dplyr::filter(gene.lengths>200) %>% 
  ggplot() +
  geom_bar(aes(x=gene, y=percentile-0.5, fill=`max count`, color=gene), stat="identity", position="dodge", size=1) +
  theme_minimal_hgrid() +
  scale_color_igv() +
  scale_fill_grey() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2)) +
  scale_y_continuous(breaks=c(0.5, 0.25, 0, -0.25, -0.5), labels=c("1\n+ risk", 0.75, 0.5, 0.25, "0\n- risk")) +
  ylab("Percentile among permutations") + xlab(NULL) + guides(color=FALSE)

severity.perc.bars

pp.pval.severity <- ggplot(p.table.severity %>% dplyr::filter(gene.lengths>200)) +
  geom_bar(aes(x=gene, y=-log10(`p-value`), fill=`max count`, color=gene), stat="identity", position="dodge", size=1) +
  theme_minimal_hgrid() +
  scale_color_igv() +
  scale_fill_grey() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2)) +
  ylab("-log10(p-value)")+ xlab(NULL) + guides(color=FALSE)  

pp.pval.severity

## Run AUC.R for this last plot

combo.grid.severity <- plot_grid(plot_grid(pp.burden.total.severity,
                                          plot_grid( pp.set.ms.severity[[11]], perm.plots.ms.severity[[11]]+ggtitle("")), labels=c("A", "B")),
                                plot_grid(plot_grid(pp.set.ms.severity[[16]], perm.plots.ms.severity[[16]]+ggtitle("")),
                                          pp.pval.severity, labels=c("C", "D")),
                                plot_grid(severity.perc.bars+guides(fill=FALSE), pp.pol, labels=c("E", "F")), nrow = 3)
combo.grid.severity
ggsave(plot=combo.grid.severity, filename="Plots/Final/Burden/combo_grid_severity.pdf", width=14, height=14)
