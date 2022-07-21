source("load_tables.R")
library(precrec)
library(rsample)

#data_boots <- rsample::bootstraps(data_muts_severe, 100)
#data_boots_2 <- rsample::bootstraps(data_muts_severe, 400)

#saveRDS(object=data_boots, "data_boots.RDS")
#saveRDS(object=data_boots_2, "data_boots_2.RDS")

data_boots <- readRDS("data_boots.RDS")
data_boots_2 <- readRDS("data_boots_2.RDS")

nongenetic.loo <- function(df){
  nongenetic_preds <- c()
  for(ii in 1:nrow(df)){
    data_tmp <- df[-ii,]
    fit_tmp <- glm(`Severe Disease` ~
                     age_imputed +
                     Vaccination.Status..Processed. +
                     `Other Comorbidities`+
                     Diabetes + 
                     Nationality_processed +
                     region.processed,
                   family="binomial", data_tmp)
    nongenetic_preds[ii] <- predict.glm(fit_tmp, df[ii,], type="response")
  }
  df$nongenetic_preds <- nongenetic_preds
  precrec_nongenetic <- evalmod(scores = df$nongenetic_preds, labels = df$`Severe Disease`)
  return(attr(precrec_nongenetic$rocs[[1]], "auc"))
}

auc.nongenetic <- lapply(data_boots$splits, function(boot){
  print("go")
  nongenetic.loo(data_muts_severe[boot$in_id,])
})

auc.nongenetic_2 <- lapply(data_boots_2$splits, function(boot){
  print("go")
  nongenetic.loo(data_muts_severe[boot$in_id,])
})

nongenetic_preds <- c()
for(ii in 1:nrow(data_muts_severe)){
  data_tmp <- data_muts_severe[-ii,]
  fit_tmp <- glm(`Severe Disease` ~
                   age_imputed +
                   Vaccination.Status..Processed. +
                   `Other Comorbidities`+
                   Diabetes + 
                   Nationality_processed +
                   region.processed,
                 family="binomial", data_tmp)
  nongenetic_preds[ii] <- predict.glm(fit_tmp, data_muts_severe[ii,], type="response")
}
data_muts_severe$nongenetic_preds <- nongenetic_preds

data_muts_severe$rv.count <- get.rv.burden(ms.genes$Pol, 10)

genetic_preds <- c()
for(ii in 1:nrow(data_muts_severe)){
  data_tmp <- data_muts_severe[-ii,]
  fit_tmp <- glm(`Severe Disease` ~
                   age_imputed +
                   Vaccination.Status..Processed. +
                   `Other Comorbidities`+
                   Diabetes + 
                   Nationality_processed +
                   region.processed + 
                   rv.count,
                 family="binomial", data_tmp)
  genetic_preds[ii] <- predict.glm(fit_tmp, data_muts_severe[ii,], type="response")
}
data_muts_severe$genetic_preds <- genetic_preds

genetic.loo <- function(df){
  genetic_preds <- c()
  for(ii in 1:nrow(df)){
    data_tmp <- df[-ii,]
    fit_tmp <- glm(`Severe Disease` ~
                     age_imputed +
                     Vaccination.Status..Processed. +
                     `Other Comorbidities`+
                     Diabetes + 
                     Nationality_processed +
                     region.processed + 
                     rv.count,
                   family="binomial", data_tmp)
    genetic_preds[ii] <- predict.glm(fit_tmp, df[ii,], type="response")
  }
  df$genetic_preds <- genetic_preds
  precrec_genetic <- evalmod(scores = df$genetic_preds, labels = df$`Severe Disease`)
  return(attr(precrec_genetic$rocs[[1]], "auc"))
}

genetic.loo.spike <- function(df){
  df$rv.count.spike <- get.rv.burden(ms.genes$`S glycoprotein`, 10)
  genetic_preds <- c()
  for(ii in 1:nrow(df)){
    data_tmp <- df[-ii,]
    fit_tmp <- glm(`Severe Disease` ~
                     age_imputed +
                     Vaccination.Status..Processed. +
                     `Other Comorbidities`+
                     Diabetes + 
                     Nationality_processed +
                     region.processed + 
                     rv.count.spike,
                   family="binomial", data_tmp)
    genetic_preds[ii] <- predict.glm(fit_tmp, df[ii,], type="response")
  }
  df$genetic_preds <- genetic_preds
  precrec_genetic <- evalmod(scores = df$genetic_preds, labels = df$`Severe Disease`)
  return(attr(precrec_genetic$rocs[[1]], "auc"))
}

genetic.loo.nsp2 <- function(df){
  df$rv.count.nsp2 <- get.rv.burden(ms.genes$nsp2, 10)
  genetic_preds <- c()
  for(ii in 1:nrow(df)){
    data_tmp <- df[-ii,]
    fit_tmp <- glm(`Severe Disease` ~
                     age_imputed +
                     Vaccination.Status..Processed. +
                     `Other Comorbidities`+
                     Diabetes + 
                     Nationality_processed +
                     region.processed + 
                     rv.count.nsp2,
                   family="binomial", data_tmp)
    genetic_preds[ii] <- predict.glm(fit_tmp, df[ii,], type="response")
  }
  df$genetic_preds <- genetic_preds
  precrec_genetic <- evalmod(scores = df$genetic_preds, labels = df$`Severe Disease`)
  return(attr(precrec_genetic$rocs[[1]], "auc"))
}

genetic.loo.pc1 <- function(df){
  genetic_preds <- c()
  for(ii in 1:nrow(df)){
    data_tmp <- df[-ii,]
    fit_tmp <- glm(`Severe Disease` ~
                     age_imputed +
                     Vaccination.Status..Processed. +
                     `Other Comorbidities`+
                     Diabetes + 
                     Nationality_processed +
                     region.processed + 
                     pc.coding.1,
                   family="binomial", data_tmp)
    genetic_preds[ii] <- predict.glm(fit_tmp, df[ii,], type="response")
  }
  df$genetic_preds <- genetic_preds
  precrec_genetic <- evalmod(scores = df$genetic_preds, labels = df$`Severe Disease`)
  return(attr(precrec_genetic$rocs[[1]], "auc"))
}

auc.genetic <- lapply(data_boots$splits, function(boot){
  genetic.loo(data_muts_severe[boot$in_id,])
})

auc.genetic_2 <- lapply(data_boots_2$splits, function(boot){
  genetic.loo(data_muts_severe[boot$in_id,])
})

auc.genetic.spike <- lapply(data_boots$splits, function(boot){
  genetic.loo.spike(data_muts_severe[boot$in_id,])
})

auc.genetic.spike_2 <- lapply(data_boots_2$splits, function(boot){
  genetic.loo.spike(data_muts_severe[boot$in_id,])
})

auc.genetic.nsp2 <- lapply(data_boots$splits, function(boot){
  genetic.loo.nsp2(data_muts_severe[boot$in_id,])
})

auc.genetic.nsp2_2 <- lapply(data_boots_2$splits, function(boot){
  genetic.loo.nsp2(data_muts_severe[boot$in_id,])
})

auc.genetic.pc1 <- lapply(data_boots$splits, function(boot){
  genetic.loo.pc1(data_muts_severe[boot$in_id,])
})

auc.genetic.pc1_2 <- lapply(data_boots_2$splits, function(boot){
  genetic.loo.pc1(data_muts_severe[boot$in_id,])
})

write.table(data.frame(auc_genetic=c(unlist(auc.genetic), unlist(auc.genetic_2)),
                       auc_nongenetic=c(unlist(auc.nongenetic), unlist(auc.nongenetic_2))), 
            file = "AUC_Pol.csv", sep=",", row.names = F)

Pol.table <- read.csv(file = "AUC_Pol.csv", sep = ",")

write.table(data.frame(auc_genetic=c(unlist(auc.genetic.spike), unlist(auc.genetic.spike_2)),
                       auc_nongenetic=c(unlist(auc.nongenetic), unlist(auc.nongenetic_2))), 
            file = "AUC_Spike.csv", sep=",", row.names = F)

write.table(data.frame(auc_genetic=c(unlist(auc.genetic.nsp2), unlist(auc.genetic.nsp2_2)),
                       auc_nongenetic=c(unlist(auc.nongenetic), unlist(auc.nongenetic_2))), 
            file = "AUC_nsp2.csv", sep=",", row.names = F)

write.table(data.frame(auc_genetic=c(unlist(auc.genetic.pc1), unlist(auc.genetic.pc1_2)),
                       auc_nongenetic=c(unlist(auc.nongenetic), unlist(auc.nongenetic_2))), 
            file = "AUC_pc1.csv", sep=",", row.names = F)


pp.spike <- ggplot(data=data.frame(delta.auc=c(unlist(auc.genetic.spike), unlist(auc.genetic.spike_2))-c(unlist(auc.nongenetic), unlist(auc.nongenetic_2))),
       mapping=aes(x=delta.auc)) + 
  geom_histogram(aes(y=..density..), bins=15, fill=gene.colors.alt["S glycoprotein"], alpha=0.7) +
  geom_density(size=2, color=gene.colors.alt["S glycoprotein"]) + 
  theme_cowplot() + 
  xlab(expression(paste(Delta, "AUC"))) + ylab("density") + 
  ggtitle("S glycoprotein: 95%CI=[-0.003, 0.003]") + theme(legend.position = "none") + 
  geom_vline(xintercept=spike.auc-attr(precrec_nongenetic$rocs[[1]], "auc"),
             color=gene.colors.alt["S glycoprotein"], lwd=2)
ggsave("Plots/Final/Spike_AUC.png")

spike.auc <- genetic.loo.spike(data_muts_severe)

pp.pol <- ggplot(data=data.frame(delta.auc=c(unlist(auc.genetic), unlist(auc.genetic_2))-c(unlist(auc.nongenetic), unlist(auc.nongenetic_2))),
       mapping=aes(x=delta.auc)) + 
  geom_histogram(aes(y=..density..), bins=15, fill=gene.colors.alt["Pol"], alpha=0.7) +
  geom_density(size=2, color=gene.colors.alt["Pol"]) + 
  theme_cowplot() + 
  xlab(expression(paste(Delta, "AUC"))) + ylab("density") + 
  ggtitle("Pol: 95%CI=[0.0012, 0.017]") + theme(legend.position = "none") + 
  geom_vline(xintercept=attr(precrec_genetic$rocs[[1]], "auc")-attr(precrec_nongenetic$rocs[[1]], "auc"),
             color=gene.colors.alt["Pol"], lwd=2)
ggsave("Plots/Final/Pol_AUC.png")

nsp2.auc <- genetic.loo.nsp2(data_muts_severe)

pp.nsp2 <- ggplot(data=data.frame(delta.auc=c(unlist(auc.genetic.nsp2), unlist(auc.genetic.nsp2_2))-c(unlist(auc.nongenetic), unlist(auc.nongenetic_2))),
       mapping=aes(x=delta.auc)) + 
  geom_histogram(aes(y=..density..), bins=15, fill=gene.colors.alt["nsp2"], alpha=0.7) +
  geom_density(size=2, color=gene.colors.alt["nsp2"]) + 
  theme_cowplot() + 
  xlab(expression(paste(Delta, "AUC"))) + ylab("density") + 
  ggtitle("nsp2: 95%CI=[-0.003, 0.004]") + theme(legend.position = "none") +
  geom_vline(xintercept=nsp2.auc-attr(precrec_nongenetic$rocs[[1]], "auc"),
             color=gene.colors.alt["nsp2"], lwd=2) 

pc1.auc <- genetic.loo.pc1(data_muts_severe)
pp.pc1 <- ggplot(data=data.frame(delta.auc=c(unlist(auc.genetic.pc1), unlist(auc.genetic.pc1_2))-c(unlist(auc.nongenetic), unlist(auc.nongenetic_2))),
       mapping=aes(x=delta.auc)) + 
  geom_histogram(aes(y=..density..), bins=15, fill="black", alpha=0.7) +
  geom_density(size=2, color="black") + 
  theme_cowplot() + 
  xlab(expression(paste(Delta, "AUC"))) + ylab("density") + 
  ggtitle("PC1: 95%CI=[-0.003, 0.01]") + theme(legend.position = "none") +
  geom_vline(xintercept=pc1.auc-attr(precrec_nongenetic$rocs[[1]], "auc"),
             color="black", lwd=2) 

cowplot::plot_grid(pp.pol, pp.spike, pp.nsp2, pp.pc1)
ggsave("Plots/Final/auc_bootstrap_grid.png", height=7, width=12)

quantile(c(unlist(auc.genetic), unlist(auc.genetic_2))-c(unlist(auc.nongenetic), unlist(auc.nongenetic_2)), c(0.025, 0.975))
quantile(c(unlist(auc.genetic.spike), unlist(auc.genetic.spike_2))-c(unlist(auc.nongenetic), unlist(auc.nongenetic_2)), c(0.025, 0.975))
quantile(c(unlist(auc.genetic.nsp2), unlist(auc.genetic.nsp2_2))-c(unlist(auc.nongenetic), unlist(auc.nongenetic_2)), c(0.025, 0.975))
quantile(c(unlist(auc.genetic.pc1), unlist(auc.genetic.pc1_2))-c(unlist(auc.nongenetic), unlist(auc.nongenetic_2)), c(0.025, 0.975))

