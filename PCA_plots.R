source("load_tables.R")

pp.pca.coding <- make.pca.plot(data_muts_severe, "pc.coding.") + theme(legend.position = c(0.5,1),
                                                                                        legend.justification = c(0,1))
pp.pca.syn <- make.pca.plot(data_muts_severe, "pc.synonymous.") + ggtitle("Synonymous") + theme(legend.position = "none")
pp.pca.ms <- make.pca.plot(data_muts_severe, "pc.missense.") + ggtitle("Missense") + theme(legend.position = "none")

cowplot::plot_grid(pp.pca.coding, pp.pca.syn, pp.pca.ms, nrow = 1)
ggsave("Plots/Final/pca_coding_compare.pdf", width=14, height=5)
make.pca.plot(data_muts_severe, "pc.coding.")
ggsave("Plots/Final/PCA_clade_coding.pdf", width=7, height=5)

data_muts_severe$Collection.Range <- as.Date(cut(data_muts_severe$Collection.Date , breaks = "7 days"))

pp.clade.time <- ggplot(data_muts_severe, aes(x=Collection.Range, fill=clade)) + geom_bar() +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y") +
  scale_fill_igv() + theme_cowplot() + ylab("Case count") + xlab(NULL)

date.severity <- ggplot(data_muts_severe) +
  geom_jitter(aes(x=Collection.Date, y=`Severe Disease`), height=0.1) + theme_minimal_hgrid() +
  geom_smooth(aes(x=Collection.Date, y=`Severe Disease`), se=T, method="loess") +
  ylim(-.1, 1.1) + scale_y_continuous(name="Severe COVID19",breaks=c(0, 0.1,0.3, 0.5, 0.7, 0.9, 1),
                                      labels=c("No", 0.1, 0.3, 0.5, 0.7, 0.9,  "Yes")) +
  xlab(NULL)
date.severity

cowplot::plot_grid(pp.clade.time, plot_grid(pp.pca.coding + theme(legend.position = "none"), date.severity, labels=c("B", "C")), nrow=2,
                   labels=c("A", ""))
ggsave(filename="Plots/Final/basics_grid.pdf", width=10, height=7)

cowplot::plot_grid(plotlist=list(pp.pca.coding, pp.pca.syn, pp.pca.ms), nrow = 1)

plotlist_pcgrid <- list()
for(ii in 1:7){
  for(jj in 1:7){
    if(ii != jj){
      plotlist_pcgrid[[paste(ii, jj, sep="_")]] <- make.pca.plot(data_muts_severe, "pc.coding.", pc1 = ii, pc2 = jj, PC.label="PC") +
        theme(legend.position = "none")
    }
  }
}
pp.legend <- cowplot::get_legend(ggplot(data_muts_severe) + geom_point(aes(x=1, y=2, color=clade), size=2, alpha=0.75) +
                                   theme_cowplot() + scale_color_igv() + guides(color=guide_legend(ncol=2)))
pp.empty <- ggplot()+theme_void()
pp.grid <- plot_grid(plotlist_pcgrid[["1_2"]], pp.legend , pp.empty, pp.empty, pp.empty, pp.empty,
                     plotlist_pcgrid[["1_3"]], plotlist_pcgrid[["2_3"]],  pp.empty, pp.empty, pp.empty, pp.empty,
                     plotlist_pcgrid[["1_4"]], plotlist_pcgrid[["2_4"]], plotlist_pcgrid[["3_4"]], pp.empty, pp.empty,  pp.empty,
                     plotlist_pcgrid[["1_5"]], plotlist_pcgrid[["2_5"]], plotlist_pcgrid[["3_5"]], plotlist_pcgrid[["4_5"]], pp.empty, pp.empty,
                     plotlist_pcgrid[["1_6"]], plotlist_pcgrid[["2_6"]], plotlist_pcgrid[["3_6"]], plotlist_pcgrid[["4_6"]], plotlist_pcgrid[["5_6"]],  pp.empty,
                     plotlist_pcgrid[["1_7"]], plotlist_pcgrid[["2_7"]], plotlist_pcgrid[["3_7"]], 
                     plotlist_pcgrid[["4_7"]], plotlist_pcgrid[["5_7"]], plotlist_pcgrid[["6_7"]], 
                     ncol=6, nrow=6)
pp.grid
pp.grid <- plot_grid(plotlist_pcgrid[["1_2"]], pp.legend , pp.empty, pp.empty,
                     plotlist_pcgrid[["1_3"]], plotlist_pcgrid[["2_3"]],  pp.empty, pp.empty,
                     plotlist_pcgrid[["1_4"]], plotlist_pcgrid[["2_4"]], plotlist_pcgrid[["3_4"]], pp.empty,
                     plotlist_pcgrid[["1_5"]], plotlist_pcgrid[["2_5"]], plotlist_pcgrid[["3_5"]], plotlist_pcgrid[["4_5"]], ncol=4, nrow=4)
ggsave(filename = "Plots/Final/PCA_grid_coding.pdf", plot = pp.grid, height=12, width=12)

pp.loadings  <- lapply(1:7, function(X) plot.pca.loadings(all.gene.coding, pc = X))

cowplot::plot_grid(plotlist=pp.loadings, ncol=2)
ggsave("Plots/Final/PCs/loading_grid_all_coding.pdf", width=14, height=20)
