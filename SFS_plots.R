source("load_tables.R")

## Generate SFS rank and percentile plots

sfs.rank.ms <- lapply(genes, function(gene){
  print(gene)
  rank.sfs.plot(ms.genes[[gene]], all.gene.ms, gene.name=paste(gene, "missense"), reps = 200, color = gene.colors.alt[[gene]])
})
names(sfs.rank.ms) <- genes

ms.lower.perc <- lapply(genes, function(gene.name){
  print(gene.name)
  get.rank.lower.alt(ms.genes[[gene.name]], all.gene.ms, reps=1000, reps.2=2000, rep.mult = 1)
})
names(ms.lower.perc) <- genes

## Calculate percentiles among randomizations for derived allele burden, nucleotide diversity and rank sum percentile

ms.lower.pvals <- data.frame(do.call(rbind, lapply(ms.lower.perc, p.val.calc)))
ms.lower.pvals$analysis <- "MS/MS"

all.pvals <- do.call(rbind, list(ms.lower.pvals))
all.pvals$p.val.dab <- unlist(all.pvals$p.val.dab)
all.pvals$p.val.pi <- unlist(all.pvals$p.val.pi)
all.pvals$p.val.rps <- unlist(all.pvals$p.val.rps)
all.pvals$gene <- rownames(all.pvals)

{
ms.lower.pvals$nmut <- lapply(ms.genes, ncol)
ms.lower.pvals$prop.mut <- unlist(ms.lower.pvals$nmut)/(3*gene.lengths/1000)
ms.lower.pvals$gene <- rownames(ms.lower.pvals)
pp.ms <- ggplot(unlist.dframe(ms.lower.pvals)[gene.lengths>200,]) + geom_label(aes(x=p.val.dab, y=p.val.rps, label=gene,
                                                      color=prop.mut), alpha=0.5) +
  xlim(0, 1.1) + ylim(0, 1.1) +
  scale_x_sqrt(limits=c(0.0025, 1.3), breaks=c(0.01, 0.1, 0.2, 0.5, 1.0)) +
  scale_y_sqrt(limits=c(0.0025, 1.3), breaks=c(0.01, 0.1, 0.2, 0.5, 1.0)) +
  theme_minimal_grid() + scale_color_gradient(low = "#d01c8b", high="#4dac26") +
  xlab("Derived allele burden percentile") + ylab("RPSS") +
  labs(color = "MS muts / \n gene length (kb)")
pp.ms
ggsave(plot = pp.ms, filename = "Plots/Final/SFS/ms_sel_percentiles.pdf")

## Plot examples for particular genes and add the percentiles scatterplot

ms.grid <- plot_grid(plot_grid(sfs.rank.ms[["Pol"]],
                      lower.perc.plot(ms.lower.perc[[genes[11]]], color=gene.colors.alt[11], gene.name = NULL, add.p = FALSE),
            nrow = 2),
            plot_grid(sfs.rank.ms[["S glycoprotein"]],
                      lower.perc.plot(ms.lower.perc[[genes[16]]], color=gene.colors.alt[16], gene.name = NULL, add.p = FALSE),
                                 nrow = 2),
            plot_grid(sfs.rank.ms[["3CL-PRO"]],
                      lower.perc.plot(ms.lower.perc[[genes[5]]], color=gene.colors.alt[5], gene.name = NULL, add.p = FALSE),
                      nrow = 2),
            pp.ms, labels = c("A", "B", "C", "D"), label_size = 20)
ggsave(ms.grid, filename = "Plots/Final/SFS/ms_sel_grid.pdf", width=12, height = 8)
}

# saveRDS(ms.lower.perc, "ms_lower_perc.RDS")

# ms.lower.perc <- readRDS("ms_lower_perc.RDS")

## Save figures for the whole pile of SFS plots

rank.perc.plots.ms <- lapply(1:length(genes), function(ii){
  lower.perc.plot(ms.lower.perc[[genes[ii]]], color=gene.colors.alt[ii], gene.name = genes[ii], add.p = F)
})

cowplot::plot_grid(plotlist=rank.perc.plots.ms[gene.lengths>200], ncol=3)
ggsave("Plots/Final/SFS/ms_rank_percentile_grid.pdf", height=10, width=15)


for(ii in 1:length(genes)){
  pp <- lower.perc.plot(ms.lower.perc[[genes[ii]]], color=gene.colors.alt[ii], gene.name = genes[ii], add.p = F)
  ggsave(filename = file.path("Plots/Final/SFS", paste0("ms_", genes[ii], "_rank_percentile.pdf")), plot = pp)
}

for(ii in 1:length(genes)){
  pp <- lower.perc.plot(ms.lower.perc.100[[genes[ii]]], color=gene.colors.alt[ii], gene.name = genes[ii], add.p = F)
  ggsave(filename = file.path("Plots/Final/SFS", paste0("ms_", genes[ii], "_rank_percentile_100.pdf")), plot = pp)
}

for(ii in 1:length(genes)){
  ggsave(filename = file.path("Plots/Final/SFS", paste0("ms_", genes[ii], "_rank_sfs.pdf")), plot = sfs.rank.ms[[genes[ii]]])
}

cowplot::plot_grid(plotlist=sfs.rank.ms, ncol=3)
ggsave("Plots/Final/SFS/ms_rank_sfs_grid.pdf", height=20, width=20)


