library(ggplot2)
library(cowplot)
library(ggsci)

unlist.dframe <- function(dat){
  for(name in names(dat)){
    dat[[name]] <- unlist(dat[[name]])
  }
  return(dat)
}

filter.positions <- function(geno.mat, filter.pos){
  var.pos <- as.integer(substring(sapply(strsplit(colnames(geno.mat), ";"), "[[", 2), 2, 
                                  nchar(sapply(strsplit(colnames(geno.mat), ";"), "[[", 2))-1))
  return(geno.mat[,! var.pos %in% filter.pos])
}

make.geno.mat <- function(dat, gene.names, var.type, lookup.table, kill.zeros = TRUE){
  col.names <- do.call(c, lapply(gene.names, function(gene.name){
    colnames(dat)[grepl(paste0("^", gene.name, "_"), colnames(dat))]
  }))
  if (length(col.names)==0)
    return(NULL)
  var.names <- sapply(stringr::str_split(col.names, "_"), "[", 2)
  var.annos <- lookup.table[var.names,]
  if (var.type == "synonymous"){
    cat("Making matrix for synonymous variants...\n")
    use.vars <- (var.annos$wt == var.annos$alt) & !is.na(var.annos$wt)
  } else if (var.type == "missense"){
    cat("Making matrix for missense variants...\n")
    use.vars <- (var.annos$wt != var.annos$alt) & !is.na(var.annos$wt) & (var.annos$alt != "*")
  } else if (var.type == "nonsense"){
    cat("Making matrix for nonsense variants...\n")
    use.vars <- (var.annos$alt == "*") & !is.na(var.annos$wt)
  } else if (var.type == "coding"){
    cat("Making matrix for all coding variants...\n")
    use.vars <- !is.na(var.annos$wt)
  } else{
    use.vars <- rep(TRUE, length(var.annos))
  }
  var.mat <- as.matrix(dat[,col.names[use.vars]])
  rownames(var.mat) <- rownames(dat[,col.names[use.vars]])
  var.info <- paste(paste0("gene:", var.annos[use.vars,]$gene), var.annos[use.vars,]$mutation,
                    var.annos[use.vars,]$wt, var.annos[use.vars,]$alt, sep=";")
  colnames(var.mat) <- var.info
  if( kill.zeros )
    var.mat <- var.mat[,colSums(var.mat) != 0]
  return(var.mat)
}

get.rv.counts <- function(geno.mat, max.count){
  allele.counts <- colSums(geno.mat)
  use.cols <- allele.counts <= max.count
  return(colSums(geno.mat[,use.cols]))
}

get.rv.burden <- function(geno.mat, max.count){
  if(is.vector(geno.mat)){
    allele.count <- sum(geno.mat)
    if(allele.count <= max.count){
      return(geno.mat)
    } else{
      return(rep(0, length(geno.mat)))
    }
  }
  allele.counts <- colSums(geno.mat)
  use.cols <- allele.counts <= max.count
  # print(c(max.count, Inf, unname(use.cols)))
  if(sum(use.cols)==1)
    return(geno.mat[,use.cols])
  return(rowSums(geno.mat[,use.cols]))
}

make.rv.samples <- function(geno.mat, num.vars, reps=1000){
  result <- list()
  total.vars <- ncol(geno.mat)
  for(ii in 1:reps){
    use.vars <- sample.int(total.vars, size=num.vars)
    result[[ii]] <- sort(unname(colSums(geno.mat[,use.vars])), decreasing=TRUE)
  }
  return(result)
}

rv.burden.plot <- function(data.test, color="black", title=NULL, xlimits=NULL){
  pp <- data.test %>% dplyr::filter(se<10) %>%
    ggplot() +
    geom_point(aes(x=beta, y=count.max), color=color) +
    geom_errorbarh(aes(xmin=beta-1.96*se, xmax=beta+1.96*se, y=count.max), color=color) +
    scale_y_log10(limits=c(0.5,120)) + theme_cowplot() +
    geom_vline(aes(xintercept=0)) +
    theme(legend.position = "none") +
    ylab("Maximum frequency in burden") +
    ggtitle(title)
  if(!is.null(xlimits))
    pp <- pp + xlim(xlimits)
  return(pp)
}

rank.sfs.plot <- function(gene.mat, comparison.mat, max.count=Inf, gene.name=NULL, reps=100, color="black"){
  gene.counts <- get.rv.counts(gene.mat, max.count=Inf)
  comparison.counts <- make.rv.samples(comparison.mat, num.vars = length(get.rv.counts(gene.mat, max.count=max.count)), reps = reps)
  pp <- ggplot()
  for(ii in 1:reps){
    df <- data.frame(comparison.counts[[ii]], 1:ncol(gene.mat))
    colnames(df) <- c("freq", "rank")
    pp <- pp + geom_line(data=df, aes(y=freq, x=rank), alpha=0.02)
  }
  pp <- pp + scale_y_log10() + scale_x_log10()
  df <- data.frame(freq=sort(unname(colSums(gene.mat)), decreasing=T), rank=1:ncol(gene.mat))
  pp <- pp + geom_line(data=df, aes(y=freq, x=rank), alpha=1, color=color, lwd=2) +
    geom_point(data=df, aes(y=freq, x=rank), alpha=1, color=color, size=3)+
    theme_cowplot() +
    ylab("Frequency") + ggtitle(gene.name)
  return(pp)
}

get.rank.percentiles <- function(gene.samps, gene.mat){
  samp.mat <- do.call(rbind, gene.samps)
  rank.sfs <- sort(unname(colSums(gene.mat)), decreasing=T)
  result <- rank.sfs*0
  for(ii in 1:ncol(gene.mat)){
    ecdf.tmp <- ecdf(samp.mat[,ii])
    result[ii] <- ecdf.tmp(rank.sfs[ii])
  }
  return(result)
}

get.rank.lower <- function(gene.mat, comparison.mat, reps=1000, reps.2=1000){
  samps.mat <- do.call(rbind, make.rv.samples(comparison.mat, ncol(gene.mat), reps=reps))
  rank.sfs <- sort(unname(colSums(gene.mat)), decreasing=T)
  result <- rank.sfs*0
  for(ii in 1:ncol(gene.mat)){
    if(ii==1){
      result[ii] <- sum(samps.mat[,ii] <= rank.sfs[ii])/length(samps.mat[,ii])
    } else{
      comparison.use <- samps.mat[samps.mat[,ii-1] >= rank.sfs[ii-1], ii] ## Only compare against samples with previous entry as high or higher
      if(rank.sfs[ii-1] %in% comparison.use) ## Delete previous entry so can't be resampled
        comparison.use <- comparison.use[-which(comparison.use==rank.sfs[ii-1])[1]]
      result[ii] <- sum(comparison.use <= rank.sfs[ii])/length(comparison.use)
    }
  }
  samps.mat.2 <- do.call(rbind, make.rv.samples(comparison.mat, ncol(gene.mat), reps=reps.2))
  perc.mat <- samps.mat.2*0
  for(jj in 1:reps.2){
    rank.sfs <-samps.mat.2[jj,]
    for(ii in 1:ncol(gene.mat)){
      if(ii==1){
        perc.mat[jj, ii] <- sum(samps.mat[,ii] <= rank.sfs[ii])/length(samps.mat[,ii])
      } else{
        comparison.use <- samps.mat[samps.mat[,ii-1] >= rank.sfs[ii-1], ii] ## Only compare against samples with previous entry as high or higher
        #if(rank.sfs[ii-1] %in% comparison.use) ## Delete previous entry so can't be resampled
        #  comparison.use <- comparison.use[-which(comparison.use==rank.sfs[ii-1])[1]]
        perc.mat[jj, ii] <- sum(comparison.use <= rank.sfs[ii])/length(comparison.use)
      }
    }
  }
  return(list(result, perc.mat))
}

get.rank.lower.alt <- function(gene.mat, comparison.mat, reps=1000, reps.2=1000, rep.mult=2, repl=FALSE, max.count=Inf){
  gene.mat <- gene.mat[,colSums(gene.mat) <= max.count]
  comparison.mat <- comparison.mat[,colSums(comparison.mat) <= max.count]

  samps.mat <- do.call(rbind, make.rv.samples(comparison.mat, ncol(gene.mat), reps=reps.2))

  rank.sfs <- sort(unname(colSums(gene.mat)), decreasing=T)
  len.sfs <- length(rank.sfs)
  var.counts <- colSums(comparison.mat)
  cond.samps.mat <- matrix(data=0, nrow=reps, ncol=len.sfs)

  result <- rank.sfs*0
  for(ii in 1:len.sfs){
    if(ii==1){
      comparison.vars <- var.counts
      len.vars <- length(comparison.vars)
      comparison.samps <- replicate(rep.mult*len.vars, max(sample(comparison.vars, len.sfs, replace = repl)) )
      #cond.samps.mat[,1] <- comparison.samps
      result[1] <- sum(comparison.samps <= rank.sfs[1])/(rep.mult*len.vars)
    } else{
      prev.entry <- rank.sfs[ii-1]
      prev.entry.count <- sum(rank.sfs[1:(ii-1)] == prev.entry)
      if(prev.entry == 1){
        result[ii] <- 1
      } else{
        # only want to sample where valid given previous observation
        comparison.vars <- var.counts[var.counts <= prev.entry]
        # remove entries each to number of times we've seen the previous observation
        if(sum(comparison.vars==prev.entry) == 0){
          comparison.vars <- comparison.vars
        } else if(sum(comparison.vars==prev.entry) <= prev.entry.count){
          comparison.vars <- comparison.vars[-which(comparison.vars==prev.entry)]
        } else{
          comparison.vars <- comparison.vars[-which(comparison.vars==prev.entry)[1:prev.entry.count]]
        }
        len.vars <- length(comparison.vars)
        comparison.samps <- replicate(rep.mult*len.vars, max(sample(comparison.vars, len.sfs-ii+1, replace = repl)))

        result[ii] <- sum(comparison.samps <= rank.sfs[ii])/(rep.mult*len.vars)
      }
    }
  }

  perc.mat <- samps.mat*0
  for(jj in 1:reps.2){
    rank.sfs <-samps.mat[jj,]
    for(ii in 1:len.sfs){
      if(ii==1){
        comparison.vars <- var.counts
        len.vars <- length(comparison.vars)
        comparison.samps <- replicate(rep.mult*len.vars, max(sample(comparison.vars, len.sfs, replace = repl)))
        perc.mat[jj, 1] <- sum(comparison.samps <= rank.sfs[1])/(rep.mult*len.vars)
      } else{
        prev.entry <- rank.sfs[ii-1]
        prev.entry.count <- sum(rank.sfs[1:(ii-1)] == prev.entry)
        if(prev.entry == 1){
          perc.mat[jj, ii] <- 1
        } else{
          # only want to sample where valid given previous observation
          comparison.vars <- var.counts[var.counts <= prev.entry]
          # remove entries each to number of times we've seen the previous observation
          if(sum(comparison.vars==prev.entry) == 0){
            comparison.vars <- comparison.vars
          } else if(sum(comparison.vars==prev.entry) <= prev.entry.count){
            comparison.vars <- comparison.vars[-which(comparison.vars==prev.entry)]
          } else{
            comparison.vars <- comparison.vars[-which(comparison.vars==prev.entry)[1:prev.entry.count]]
          }
          len.vars <- length(comparison.vars)
          comparison.samps <- replicate(rep.mult*len.vars, max(sample(comparison.vars, len.sfs-ii+1, replace = repl)))
          perc.mat[jj, ii] <- sum(comparison.samps <= rank.sfs[ii])/(rep.mult*len.vars)
        }
      }
    }
  }
  return(list(result, perc.mat, sort(unname(colSums(gene.mat)), decreasing=T), samps.mat, cond.samps.mat))
}

lower.perc.plot <- function(rank.lower.out, color="black", gene.name=NULL, n.lines=400, alpha=0.02, add.p=T){
  adj.mat <- rank.lower.out[[2]]
  for(ii in 1:nrow(rank.lower.out[[2]])){
    adj.mat[ii,which(rank.lower.out[[4]][ii,]==1)[-1]] <- 0
  }
  adj.rank <- rank.lower.out[[1]]
  adj.rank[which(rank.lower.out[[3]]==1)[-1]] <- 0
  p.val <- sum(rowSums(adj.mat, na.rm=T) <= sum(adj.rank))/nrow(adj.mat)
  pp <- ggplot()
  for(ii in 1:n.lines){
    keep <- adj.mat[ii,]>0
    df <- data.frame(x=(1:ncol(rank.lower.out[[2]]))[keep], y=rank.lower.out[[2]][ii,keep])
    pp <- pp + geom_line(data=df, aes(y=y, x=x), alpha=alpha)
  }
  keep <- adj.rank>0
  pp <- pp + geom_line(aes(y=rank.lower.out[[1]][keep], x=(1:length(rank.lower.out[[1]]))[keep]), color=color, lwd=2) +
    geom_point(aes(y=rank.lower.out[[1]][keep], x=(1:length(rank.lower.out[[1]]))[keep]), color=color, size=2) +
    theme_cowplot() +
    xlab("rank") + ylab("Rank count percentile") + ggtitle(paste(gene.name, "p =", p.val)) + scale_x_log10()
  if(add.p){
    pp <- pp + ggtitle(paste(gene.name, "p =", p.val))
  } else{
    pp <- pp + ggtitle(gene.name)
  }
  return(pp)
}

lower.count.plot <- function(rank.lower.out, color="black", gene.name=NULL, n.lines=400, alpha=0.02){
  pp <- ggplot()
  for(ii in 1:n.lines){
    df <- data.frame(x=(1:ncol(rank.lower.out[[5]])), y=rank.lower.out[[5]][ii,])
    pp <- pp + geom_line(data=df, aes(y=y, x=x), alpha=alpha)
  }
  pp <- pp + geom_line(aes(y=rank.lower.out[[3]], x=(1:length(rank.lower.out[[3]]))), color=color, lwd=2) +
    geom_point(aes(y=rank.lower.out[[3]], x=(1:length(rank.lower.out[[3]]))), color=color, size=2) +
    theme_cowplot() +
    xlab("rank") + ylab("Rank count percentile")+ scale_x_log10() + scale_y_log10()
  return(pp)
}

p.val.calc <- function(rank.lower.out, n.samps=988){
  ## Derived allele burden
  obs.val.dab <- sum(rank.lower.out[[3]])
  rand.vals.dab <- rowSums(rank.lower.out[[4]])
  ## pi
  obs.val.pi <- sum(rank.lower.out[[3]]/n.samps * (1-rank.lower.out[[3]]/n.samps))
  rand.vals.pi <- rowSums(rank.lower.out[[4]]/n.samps * (1-rank.lower.out[[4]]/n.samps))
  ## rank percentile sum
  adj.mat <- rank.lower.out[[2]]
  for(ii in 1:nrow(rank.lower.out[[2]])){
    adj.mat[ii,which(rank.lower.out[[4]][ii,]==1)[-1]] <- 0
  }
  adj.rank <- rank.lower.out[[1]]
  adj.rank[which(rank.lower.out[[3]]==1)[-1]] <- 0
  obs.val.rps <- sum(adj.rank)
  rand.vals.rps <- rowSums(adj.mat, na.rm=T)
  nn <- nrow(adj.mat)
  result <- list(p.val.dab = sum(rand.vals.dab <= obs.val.dab)/nn,
                 p.val.pi = sum(rand.vals.pi <= obs.val.pi)/nn,
                 p.val.rps = sum(rand.vals.rps <= obs.val.rps)/nn)
  return(result)
}

piNpiS.calc <- function(sfs.ms, sfs.syn, n.samps=988){
  ms.pi <- sum(sfs.ms/n.samps * (1-sfs.ms/n.samps))
  syn.pi <- sum(sfs.syn/n.samps * (1-sfs.syn/n.samps))
  return(ms.pi/syn.pi)
}

piNpiS.percentile <- function(d.ms, d.syn, n.samps=988){
  true.piNpiS <- piNpiS.calc(d.ms[[3]], d.syn[[3]], n.samps=n.samps)
  rand.piNpiS <- rowSums(d.ms[[4]]/n.samps * (1-d.ms[[4]]/n.samps)) /
    rowSums(d.syn[[4]]/n.samps * (1-d.syn[[4]]/n.samps))
  return(sum(rand.piNpiS <= true.piNpiS)/length(rand.piNpiS)) 
}

percentile.sfs.plot <- function(gene.mats, comparison.mat, gene.names=NULL, colors="black", reps=100, max.count=Inf){
  names(gene.mats) <- gene.names
  gene.counts <- lapply(gene.mats, function(gene.mat){
    get.rv.counts(gene.mat, max.count = max.count)
  }); names(gene.counts) <- gene.names

  gene.samps <- lapply(gene.mats, function(gene.mat){
    make.rv.samples(comparison.mat, num.vars = length(get.rv.counts(gene.mat, max.count=max.count)), reps = reps)
  }); names(gene.samps) <- gene.names

  gene.percentiles <- lapply(gene.names, function(gene.name){
    get.rank.percentiles(gene.samps[[gene.name]], gene.mats[[gene.name]])
  }); names(gene.percentiles) <- gene.names

  df.percentile <- data.frame(
    percentile=do.call(c, gene.percentiles),
    gene.name=do.call(c,lapply(gene.names, function(gene.name){
      rep(gene.name, ncol(gene.mats[[gene.name]]))
    })),
    rank=do.call(c, lapply(gene.names, function(gene.name){
      (1:ncol(gene.mats[[gene.name]]))/ncol(gene.mats[[gene.name]])
    }))
  )
  pp.perc <- ggplot(df.percentile, aes(x=rank, y=percentile, color=gene.name)) + geom_line(lwd=2, alpha=0.75) +
    theme_cowplot() + scale_color_igv() + labs(color="gene") + xlab("relative rank") +
    theme(legend.position = c(0.7, 0.5))
  return(list(pp.perc, df.percentile))
}

is_severe <- function(value) {
  if (is.na(value)) {
    return(NA)
  } else if (value == "Severe" || value == "Critical") {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

relabel_severe <- function(data) {
  new_severe_disease <- vector(mode="logical", length=nrow(data))
  for (i in 1:nrow(data)) {
    if (is_severe(data$`Disease.Severity`[i])) {
      new_severe_disease[i] <- 1;
    } else {
      new_severe_disease[i] <- 0;
    }
  }
  data["Severe Disease"] <- new_severe_disease
  return(data)
}

get.rv.sequence.severity <- function(seq_data, geno.mat, overall.mat=NULL, max.freq=100, use.clade=TRUE, use.PC=TRUE){
  rv.betas <- c()
  rv.ses <- c() 
  rv.ts <- c()
  rv.p <- c()
  seq_data_tmp <- seq_data
  freq.seq <- sort(unique(unname(colSums(geno.mat))))
  freq.seq <- freq.seq[freq.seq<=max.freq & freq.seq > 0]
  if(!is.null(overall.mat))
    seq_data_tmp$rv.count.syn.big <- get.rv.burden(overall.mat, max.freq)
  for(ii in freq.seq){
    seq_data_tmp$rv.count <- get.rv.burden(geno.mat, ii)
    if(!is.null(overall.mat))
      seq_data_tmp$rv.count.syn <- get.rv.burden(overall.mat, ii)
    if(use.clade){
      if(!is.null(overall.mat)){
        if(use.PC){
          model.fit <- glm(`Severe Disease` ~
                             age_imputed +
                             Vaccination.Status..Processed. +
                             `Other Comorbidities`+
                             Diabetes + 
                             Hypertension +
                             Ever.Smoke +
                             Gender.num +
                             day_imputed +
                             Nationality_processed +
                             region.processed +
                             rv.count +
                             rv.count.syn.big +
                             rv.count.syn +
                             clade_grouped +
                             pc.coding.1 + pc.coding.2 + pc.coding.3 + pc.coding.4 + pc.coding.5 + pc.coding.6 + pc.coding.7,
                           family="binomial", seq_data_tmp)
        } else{
          model.fit <- glm(`Severe Disease` ~
                             age_imputed +
                             Vaccination.Status..Processed. +
                             `Other Comorbidities`+
                             Diabetes + 
                             Hypertension +
                             Ever.Smoke +
                             Gender.num +
                             day_imputed +
                             Nationality_processed +
                             region.processed +
                             rv.count +
                             rv.count.syn.big +
                             rv.count.syn +
                             clade_grouped,
                           family="binomial", seq_data_tmp)
        }
    } else{
      if(use.PC){
        model.fit <- glm(`Severe Disease` ~
                           age_imputed +
                           Vaccination.Status..Processed. +
                           `Other Comorbidities`+
                           Diabetes +
                           Hypertension +
                           Ever.Smoke +
                           Gender.num +
                           day_imputed +
                           Nationality_processed +
                           region.processed +
                           rv.count +
                           clade_grouped +
                           pc.coding.1 + pc.coding.2 + pc.coding.3 + pc.coding.4 + pc.coding.5 + pc.coding.6 + pc.coding.7,
                         family="binomial", seq_data_tmp)
      } else{
        model.fit <- glm(`Severe Disease` ~
                           age_imputed +
                           Vaccination.Status..Processed. +
                           `Other Comorbidities`+
                           Diabetes +
                           Hypertension +
                           Ever.Smoke +
                           Gender.num +
                           day_imputed +
                           Nationality_processed +
                           region.processed +
                           rv.count +
                           clade_grouped,
                         family="binomial", seq_data_tmp)
      }
    }
    } else{
      if(!is.null(overall.mat)){
        if(use.PC){
          model.fit <- glm(`Severe Disease` ~
                             age_imputed +
                             Vaccination.Status..Processed. +
                             `Other Comorbidities`+
                             Diabetes +
                             Hypertension +
                             Ever.Smoke +
                             Gender.num +
                             day_imputed +
                             Nationality_processed +
                             region.processed +
                             rv.count +
                             rv.count.syn.big +
                             rv.count.syn +
                             pc.coding.1 + pc.coding.2 + pc.coding.3 + pc.coding.4 + pc.coding.5 + pc.coding.6 + pc.coding.7,
                           family="binomial", seq_data_tmp)
        } else{
          model.fit <- glm(`Severe Disease` ~
                             age_imputed +
                             Vaccination.Status..Processed. +
                             `Other Comorbidities`+
                             Diabetes +
                             Hypertension +
                             Ever.Smoke +
                             Gender.num +
                             day_imputed +
                             Nationality_processed +
                             region.processed +
                             rv.count +
                             rv.count.syn.big +
                             rv.count.syn,
                           family="binomial", seq_data_tmp)
        }
      } else{
        if(use.PC){
          model.fit <- glm(`Severe Disease` ~
                             age_imputed +
                             Vaccination.Status..Processed. +
                             `Other Comorbidities`+
                             Diabetes+
                             Hypertension +
                             Ever.Smoke +
                             Gender.num+
                             day_imputed +
                             Nationality_processed +
                             region.processed +
                             rv.count +
                             pc.coding.1 + pc.coding.2 + pc.coding.3 + pc.coding.4 + pc.coding.5 + pc.coding.6 + pc.coding.7,
                           family="binomial", seq_data_tmp)
        } else{
          model.fit <- glm(`Severe Disease` ~
                             age_imputed +
                             Vaccination.Status..Processed. +
                             `Other Comorbidities`+
                             Diabetes+
                             Hypertension +
                             Ever.Smoke +
                             Gender.num+
                             day_imputed +
                             Nationality_processed +
                             region.processed +
                             rv.count,
                           family="binomial", seq_data_tmp)
        }
      }
    }
    fit.summ <- summary(model.fit)
    rv.betas <- c(rv.betas, fit.summ$coefficients["rv.count",1])
    rv.ses <- c(rv.ses, fit.summ$coefficients["rv.count",2])
    rv.ts <- c(rv.ts, fit.summ$coefficients["rv.count",3])
    rv.p <- c(rv.p, fit.summ$coefficients["rv.count",4])
  }
  return(data.frame(beta=rv.betas, se=rv.ses, z.value=rv.ts, p.value=rv.p, count.max=freq.seq))
}

match.rv.samps <- function(seq_data, geno.mat, comparison.mat, match.mat=NULL, n.samps=100, max.freq=100){
  samp.betas <- c()
  samp.z <- c()
  samp.p <- c()
  samp.max <- c()

  seq_data_tmp <- seq_data
  # Use match.mat as a covariate for burden scores
  if(!is.null(match.mat))
    seq_data_tmp$rv.count.syn.big <- get.rv.burden(match.mat, max.freq)
  freq.seq <- sort(unique(unname(colSums(geno.mat))))
  freq.seq <- freq.seq[freq.seq<=max.freq & freq.seq > 0]
  print(freq.seq)
  prot.sfs <- table(unname(colSums(geno.mat)))
  sfs.entries <- as.integer(names(prot.sfs)) # entries are observed frequencies of mutations, not the number of times a frequency observed in the sample
  all.geno.counts <- unname(colSums(comparison.mat))
  for(ii in freq.seq){
    print(ii)
    if(!is.null(match.mat))
      seq_data_tmp$rv.count.syn <- get.rv.burden(match.mat, ii)
    sfs.use <- sfs.entries[sfs.entries <= ii]
    for(nn in 1:n.samps){
      geno.cols.use <- c()
      for(jj in 1:length(sfs.use)){
        sfs.entry <- sfs.use[jj]
        sfs.count <- unname(prot.sfs[jj])
        all.geno.cols <- which(all.geno.counts==sfs.entry) # get all cols matching the current frequency
        # Sample those columns given the count of that frequency in the current gene.
        if(length(all.geno.cols)==1){
          geno.cols.use <- c(geno.cols.use, all.geno.cols)
        } else{
          geno.cols.use <- c(geno.cols.use, sample(all.geno.cols, size=sfs.count))
        }
      }
      seq_data_tmp$rv.count <- get.rv.burden(comparison.mat[,geno.cols.use], ii)
      if(!is.null(match.mat)){
          model.fit <- glm(`Severe Disease` ~
                             age_imputed +
                             Vaccination.Status..Processed. +
                             `Other Comorbidities`+
                             Diabetes +
                             Hypertension +
                             Ever.Smoke +
                             Gender.num+
                             day_imputed +
                             Nationality_processed +
                             region.processed +
                             rv.count +
                             rv.count.syn.big +
                             rv.count.syn +
                             clade_grouped +
                             pc.coding.1 + pc.coding.2 + pc.coding.3 + pc.coding.4 + pc.coding.5 + pc.coding.6 + pc.coding.7,
                           family="binomial", seq_data_tmp)
        } else{
          model.fit <- glm(`Severe Disease` ~
                             age_imputed +
                             Vaccination.Status..Processed. +
                             `Other Comorbidities`+
                             Diabetes+
                             Hypertension +
                             Ever.Smoke +
                             Gender.num+
                             day_imputed +
                             Nationality_processed +
                             region.processed +
                             rv.count +
                             clade_grouped +
                             pc.coding.1 + pc.coding.2 + pc.coding.3 + pc.coding.4 + pc.coding.5 + pc.coding.6 + pc.coding.7,
                           family="binomial", seq_data_tmp)
      }
      fit.summ <- summary(model.fit)
      if(rownames(fit.summ$coefficients)[17]=="rv.count"){
        samp.betas <- c(samp.betas, fit.summ$coefficients["rv.count",1])
        samp.z <- c(samp.z, fit.summ$coefficients["rv.count",3])
        samp.p <- c(samp.p, fit.summ$coefficients["rv.count",4])
      } else{
        samp.betas <- c(samp.betas, NA)
        samp.z <- c(samp.z, NA)
        samp.p <- c(samp.p, NA)
      }
      samp.max <- c(samp.max, ii)
    }
  }
  return( data.frame(beta=samp.betas, z.value=samp.z, p.value=samp.p, max.count=samp.max) )
}

get.p.perm <- function(rv.tests, rv.perms){
  result <- c()
  for(ii in 1:nrow(rv.tests)){
    perms.tmp <- dplyr::filter(rv.perms, max.count==rv.tests$count.max[ii])
    ecdf.tmp <- ecdf(perms.tmp$beta)
    beta <- rv.tests$beta[ii]
    if(beta<0){
      result <- c(result, ecdf.tmp(beta) + 1-ecdf.tmp(-beta))
    } else{
      result <- c(result, ecdf.tmp(-beta) + 1-ecdf.tmp(beta))
    }
  }
  return(result)
}

get.p.perm.z <- function(rv.tests, rv.perms){
  result <- c()
  for(ii in 1:nrow(rv.tests)){
    perms.tmp <- dplyr::filter(rv.perms, max.count==rv.tests$count.max[ii])
    ecdf.tmp <- ecdf(perms.tmp$z.value)
    z.value <- rv.tests$z.value[ii]
    if(z.value<0){
      result <- c(result, ecdf.tmp(z.value) + 1-ecdf.tmp(-z.value))
    } else{
      result <- c(result, ecdf.tmp(-z.value) + 1-ecdf.tmp(z.value))
    }
  }
  return(result)
}

get.perc.perms <- function(rv.tests, rv.perms){
  result <- c()
  for(ii in 1:nrow(rv.tests)){
    perms.tmp <- dplyr::filter(rv.perms, max.count==rv.tests$count.max[ii])
    if(sum(!is.na(perms.tmp$beta))==0){
      result <- c(result, NA)
      next
    }
    ecdf.tmp <- ecdf(perms.tmp$beta)
    med.beta <- median(perms.tmp$beta, na.rm=T)
    beta <- rv.tests$beta[ii]
    if(rv.tests$se[ii]>10){
      result <- c(result, NA)
    } else{
      if(beta < med.beta){
        result <- c(result, ecdf.tmp(beta))
      } else{
        result <- c(result, ecdf.tmp(beta-1e-8))
      }
    }
  }
  return(result)
}

get.perc.perms.z <- function(rv.tests, rv.perms){
  result <- c()
  for(ii in 1:nrow(rv.tests)){
    perms.tmp <- dplyr::filter(rv.perms, max.count==rv.tests$count.max[ii])
    print(sum(!is.na(perms.tmp$z.value)))
    if(sum(!is.na(perms.tmp$z.value))==0){
      result <- c(result, NA)
      next
    }
    ecdf.tmp <- ecdf(perms.tmp$z.value)
    med.z <- median(perms.tmp$z.value, na.rm=T)
    z.value <- rv.tests$z.value[ii]
    if(rv.tests$se[ii]>10){
      result <- c(result, NA)
    } else{
      if(z.value < med.z){
        result <- c(result, ecdf.tmp(z.value))
      } else{
        result <- c(result, ecdf.tmp(z.value-1e-8))
      }
    }
  }
  return(result)
}

plot.p.perms <- function(rv.tests, rv.perms, gene.name=NULL, color="black", xlimits=NULL, text.point=2.7){
  rv.perms <- dplyr::filter(rv.perms, abs(beta)<5)
  rv.tests$perc.perm <- get.perc.perms(rv.tests, rv.perms)
  pp <- ggplot(rv.tests) +
    geom_violin(data=rv.perms,
                 aes(x=beta, y=as.factor(max.count), group=as.factor(max.count)), alpha=0.2) +
    geom_point(aes(x=beta, y=as.factor(count.max)), color=color, size=3) +
    geom_label(aes(x=text.point, y=as.factor(count.max), label=paste(round(perc.perm, digits=3))), color="blue", size=5) +
    theme_cowplot() +
    geom_vline(aes(xintercept=0)) +
    theme(legend.position = "none") + ylab("Maximum frequency in burden") +
    ggtitle(gene.name)
  if(is.null(xlimits)){
    pp <- pp + xlim(c(-3.5, 3.5))
  } else{
    pp <- pp + xlim(xlimits)
  }
  return(pp)
}

plot.p.perms.z <- function(rv.tests, rv.perms, gene.name=NULL, color="black", xlimits=NULL, text.point=4){
  rv.perms <- dplyr::filter(rv.perms, abs(beta)<5)
  rv.tests$perc.perm <- get.perc.perms.z(rv.tests, rv.perms)
  pp <- ggplot(rv.tests) +
    geom_violin(data=rv.perms,
                aes(x=z.value, y=as.factor(max.count), group=as.factor(max.count)), alpha=0.2) +
    geom_point(aes(x=z.value, y=as.factor(count.max)), color=color, size=3) +
    geom_label(aes(x=text.point, y=as.factor(count.max), label=paste(round(perc.perm, digits=3))), color="blue", size=5) +
    theme_cowplot() +
    geom_vline(aes(xintercept=0)) +
    theme(legend.position = "none") + ylab("Maximum frequency in burden") +
    ggtitle(gene.name)
  if(is.null(xlimits)){
    pp <- pp + xlim(c(-5.0, 5.0))
  } else{
    pp <- pp + xlim(xlimits)
  }
  return(pp)
}

calc.pca <- function(geno.mat, num.pcs = 7, scale. = FALSE){
  pca.out <- prcomp(geno.mat, scale.=scale.)
  indiv.vals <- scale(geno.mat, scale=scale.) %*% pca.out$rotation[,1:num.pcs]
  return(indiv.vals)
}

plot.pca.loadings <- function(geno.mat, pc = 1, scale. = FALSE){
  pca.out <- prcomp(geno.mat, scale. = scale.)
  pc.argsort <- sort(abs(pca.out$rotation[,pc]), index.return=T, decreasing = T)$ix
  pc.load <- data.frame(pc.loading=pca.out$rotation[pc.argsort, pc],
                         pc.rank=1:ncol(geno.mat),
                         mut.name=colnames(geno.mat)[pc.argsort])
  result <- ggplot(dplyr::filter(pc.load, pc.rank<15)) +
    geom_label(aes(x=pc.rank, y=16-pc.rank, size=abs(pc.loading), label=mut.name), alpha=0.5)  + theme_cowplot() +
    xlim(c(-2,16))+ theme(legend.position = "none") + xlab("Mutation loading rank") + ggtitle(paste0("PC", pc)) + ylab(NULL) +
    theme(axis.text.y = NULL) + scale_y_continuous(breaks=NULL)
  return(result)
}


add.pcs <- function(dat, indiv.pcs, pca.id = "pc."){
  for(ii in 1:ncol(indiv.pcs)){
    pc.name <- paste0(pca.id, ii)
    dat[[pc.name]] <- indiv.pcs[,ii]
  }
  return(dat)
}

make.pca.plot <- function(dat, pca.id, pc1=1, pc2=2, PC.label="Principal component"){
  pc1.name <- ggplot2::sym(paste0(pca.id, pc1))
  pc2.name <- ggplot2::sym(paste0(pca.id, pc2))
  result <- ggplot(dat) + geom_point(aes(x=!!pc1.name, y=!!pc2.name, color=clade), size=2, alpha=0.75) +
    theme_cowplot() + scale_color_igv() +
    xlab(paste(PC.label, pc1)) + ylab(paste(PC.label, pc2)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
  return(result)
}

get.pc2 <- function(geno.mat, num.pcs = 7, scale. = FALSE){
  pca.out <- prcomp(geno.mat, scale.=scale.)
  indiv.vals <- pca.out$rotation[,2]
  return(indiv.vals)
}
