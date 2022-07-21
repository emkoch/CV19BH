source("load_tables.R")

## Grab the top 50 most common mutations
top.50 <- sort(unname(colSums(all.gene.ms)), decreasing = T, index.return=T)$ix[1:50]
top.50.names <- colnames(all.gene.ms)[top.50]

## Regress each high-frequency mutation on severity
model.fits.severe <- lapply(top.50, function(ii){
  data_muts_severe$single.mut <- all.gene.ms[,ii]
  model.fit <- broom::tidy(glm(`Severe Disease` ~
                                 age_imputed +
                                 Vaccination.Status..Processed. +
                                 `Other Comorbidities` +
                                 Diabetes +
                                 Hypertension +
                                 Ever.Smoke +
                                 Gender.num +
                                 day_imputed +
                                 Nationality_processed +
                                 region.processed +
                                 single.mut +
                                 clade_grouped +
                                 pc.coding.1 + pc.coding.2 + pc.coding.3 + pc.coding.4 + pc.coding.5 + pc.coding.6 + pc.coding.7,
                               family="binomial", data_muts_severe))
  return(model.fit)
})

## Generate QQ plot of common variant association tests
single.ests.severe <- do.call(rbind, top.ests <- lapply(model.fits.severe, function(X){
  dplyr::filter(X, term=="single.mut")
}))
single.ests.severe$mutation <- top.50.names

pdf("Plots/Final/PCs/qq_plot_common.pdf")
gaston::qqplot.pvalues(single.ests.severe$p.value[single.ests.severe$std.error<10], main=NULL)
dev.off()

## Generate plot of common variant associations with 95% CI
pp.single.severe.or <- single.ests.severe %>% dplyr::filter(std.error<10) %>%
  ggplot() +
  geom_point(aes(x=exp(estimate), y=substr(mutation, 6, 100))) +
  geom_errorbarh(aes(xmin=exp(estimate-1.96*std.error), xmax=exp(estimate+1.96*std.error), y=substr(mutation, 6, 100))) +
  theme_cowplot() + geom_vline(aes(xintercept=1), color="red") + ylab(NULL) + xlab("OR") + ggtitle("Common mutations") + 
  scale_x_log10()
pp.single.severe.or

## Fit a model with both clades and PCs
clade.fit.severe.model <- glm(`Severe Disease` ~
                               I(age_imputed/10) +
                               Vaccination.Status..Processed. +
                               `Other Comorbidities` +
                               Diabetes + 
                               Hypertension +
                               Ever.Smoke +
                               Gender.num+
                               day_imputed +
                               Nationality_processed +
                               region.processed +
                               clade_grouped +
                               pc.coding.1 + pc.coding.2 + pc.coding.3 + pc.coding.4 + pc.coding.5 + pc.coding.6 + pc.coding.7,
                             family="binomial", data_muts_severe)

## Wald tests for all caldes and all PCs
aod::wald.test(b=coef(clade.fit.severe.model), Sigma=vcov(clade.fit.severe.model), Terms=c(17:21))
aod::wald.test(b=coef(clade.fit.severe.model), Sigma=vcov(clade.fit.severe.model), Terms=c(22,28))

## Generate plots of clade and PC effects in the full model
clade.fit.severe <- broom::tidy(clade.fit.severe.model)

clade.fit.severe <- clade.fit.severe %>% mutate(or = exp(estimate),  # Odds ratio/gradient
                            var.diag = diag(vcov(clade.fit.severe.model)),  # Variance of each coefficient
                            or.se = sqrt(or^2 * var.diag),
                            or.low = exp(estimate-1.96*std.error),
                            or.high = exp(estimate+1.96*std.error))

clade.only.severe <- clade.fit.severe[grepl("clade", clade.fit.severe$term),] %>%
  tibble::add_row(term="clade_grouped20A", estimate=0, std.error=clade.fit.severe[1,]$std.error, statistic=0, p.value=0)

pp.clade.severe.or <- clade.only.severe %>% dplyr::filter(std.error<10) %>%
  ggplot() +
  geom_point(aes(x=exp(estimate), y=substr(term, 14, 100))) +
  geom_errorbarh(aes(xmin=exp(estimate-1.96*std.error), xmax=exp(estimate+1.96*std.error), y=substr(term, 14, 100))) +
  theme_cowplot() + ylab(NULL) + xlab("OR")+ geom_vline(aes(xintercept=1), color="red") + ggtitle("Nextclade designations")
pp.clade.severe.or + scale_x_log10()

pc.only.severe <- clade.fit.severe[grepl("coding", clade.fit.severe$term),]
pp.pc.severe <- pc.only.severe %>% dplyr::filter(std.error<10) %>%
  ggplot() +
  geom_point(aes(x=estimate, y=term)) +
  geom_errorbarh(aes(xmin=estimate-1.96*std.error, xmax=estimate+1.96*std.error, y=term)) +
  theme_cowplot() + ylab(NULL) + xlab("beta")+ geom_vline(aes(xintercept=0), color="red") + ggtitle("Principal Components")
pp.pc.severe

cowplot::plot_grid(pp.single.severe.or, cowplot::plot_grid(pp.clade.severe.or + scale_x_log10(), pp.pc.severe, nrow=2, 
                                                           labels=c("B", "C"), label_size=15), ncol=2,
                   labels=c("A", ""), label_size = 15)
ggsave("Plots/Final/PCs/individual_pc_grid_severity.pdf", width=10, height=8)


## Same analyses as above but use clades and PCs separately

clade.only.model <- glm(`Severe Disease` ~
                          age_imputed +
                          Vaccination.Status..Processed. +
                          `Other Comorbidities` +
                          Diabetes + 
                          Hypertension +
                          Ever.Smoke +
                          Gender.num+
                          day_imputed +
                          Nationality_processed +
                          region.processed +
                          clade_grouped,
                        family="binomial", data_muts_severe)
  
clade.only.fit <- broom::tidy(glm(`Severe Disease` ~
                                    age_imputed +
                                    Vaccination.Status..Processed. +
                                    `Other Comorbidities` +
                                    Diabetes + 
                                    Hypertension +
                                    Ever.Smoke +
                                    Gender.num+
                                    day_imputed +
                                    Nationality_processed +
                                    region.processed +
                                    clade_grouped,
                                  family="binomial", data_muts_severe))

pc.only.model <- glm(`Severe Disease` ~
      age_imputed +
      Vaccination.Status..Processed. +
      `Other Comorbidities` +
      Diabetes + 
      Hypertension +
      Ever.Smoke +
      Gender.num+
      day_imputed +
      Nationality_processed +
      region.processed +
      pc.coding.1 + pc.coding.2 + pc.coding.3 + pc.coding.4 + pc.coding.5 + pc.coding.6 + pc.coding.7,
    family="binomial", data_muts_severe)

pc.only.fit <- broom::tidy(glm(`Severe Disease` ~
                                 age_imputed +
                                 Vaccination.Status..Processed. +
                                 `Other Comorbidities` +
                                 Diabetes + 
                                 Hypertension +
                                 Ever.Smoke +
                                 Gender.num+
                                 day_imputed +
                                 Nationality_processed +
                                 region.processed +
                                 pc.coding.1 + pc.coding.2 + pc.coding.3 + pc.coding.4 + pc.coding.5 + pc.coding.6 + pc.coding.7,
                               family="binomial", data_muts_severe))


clade.only.only <- clade.only.fit[grepl("clade", clade.only.fit$term),] %>%
  tibble::add_row(term="clade_grouped20A", estimate=0, std.error=clade.fit.severe[1,]$std.error, statistic=0, p.value=0)
pp.clade.only <- clade.only.only %>% dplyr::filter(std.error<10) %>%
  ggplot() +
  geom_point(aes(x=exp(estimate), y=substr(term, 14, 100))) +
  geom_errorbarh(aes(xmin=exp(estimate-1.96*std.error), xmax=exp(estimate+1.96*std.error), y=substr(term, 14, 100))) +
  theme_cowplot() + ylab(NULL) + xlab("OR")+ geom_vline(aes(xintercept=1), color="red") + ggtitle("Nextclade designations only") + 
  scale_x_log10()
pp.clade.only
ggsave("Plots/Final/PCs/clade_only.pdf", width=5, height=4)

aod::wald.test(b=coef(clade.only.model), Sigma=vcov(clade.only.model), Terms=c(17, 20))
aod::wald.test(b=coef(pc.only.model), Sigma=vcov(pc.only.model), Terms=c(17:23))

pc.only.only <- pc.only.fit[grepl("coding", pc.only.fit$term),]
pp.pc.only <- pc.only.only %>% dplyr::filter(std.error<10) %>%
  ggplot() +
  geom_point(aes(x=exp(estimate), y=term)) +
  geom_errorbarh(aes(xmin=exp(estimate-1.96*std.error), xmax=exp(estimate+1.96*std.error), y=term)) +
  theme_cowplot() + ylab(NULL) + xlab("OR")+ geom_vline(aes(xintercept=1), color="red") + ggtitle("Principal Components only") + 
  scale_x_log10()
pp.pc.only
ggsave("Plots/Final/PCs/pc_only.pdf", width=5, height=4)

cowplot::plot_grid(pp.single.severe.or, cowplot::plot_grid(pp.clade.only, pp.pc.only, nrow=2, 
                                                           labels=c("B", "C"), label_size=15), ncol=2,
                   labels=c("A", ""), label_size = 15)
ggsave("Plots/Final/PCs/individual_pc_grid_severity_uncombined.pdf", width=10, height=8)

cowplot::plot_grid(pp.clade.only, pp.pc.only, nrow=2)
ggsave("Plots/Final/PCs/clade_pc_only.pdf", width=5, height=6)

