args = commandArgs(trailingOnly = TRUE)
filename.count<- args[1]
filename.output<- args[2]

suppressMessages(library(data.table))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(rmutil))

# counts
cnt <- read.csv(filename.count)

# filter sites with more than 2 read mapping to error allele
cnt<- subset(cnt, errors <= 2)

# optimize dispersion and error
mtmp <- function(par, x, n, ge){
  
  # likelihood, conditional on genotype error
  term1<- 0.5 * dbetabinom(x, n, m= par[1], s= par[2])
  
  term2<- 0.5 * dbetabinom(n - x, n, m= par[1], s= par[2])
  
  err.like <- (ge)*(term1 + term2)
  
  # likelihood, conditional on no genotype error
  het.like = (1-ge) * dbetabinom(x, n, m= 0.5, s= par[2])
  
  # log likelihood
  ll<- -sum(log(err.like + het.like))
  
  # return
  return(ll)
}

# optimize
m0<- optim(par= c(1e-05,1), mtmp, x= cnt$ref.matches, n= cnt$N, ge= cnt$genotype.error,
           method="L-BFGS-B", lower = c(1e-10, 1e-05), upper = c(0.999, 100))

# get parameter
err = m0$par[1]

d = m0$par[2]

# likelihood function 
ll.fun <- function(par, x, n, ge, err, d) {
  
  # likelihood for first site
  allelic.imbalance <- par
  
  # likelihood, conditional on genotype error
  t1 = 0.5 * dbetabinom(x[1], n[1], m=err, s=d)
  t2 = 0.5 * dbetabinom(n[1]-x[1], n[1], m=err, s=d)
  er1 = ge[1] * (t1 + t2)
  
  # likelihood, conditional on no genotype error
  d1 = (1 - ge[1]) * dbetabinom(x[1], n[1], m =  0.5 + allelic.imbalance, s = d)
  
  # combined likelihood
  p1 = er1 + d1
  
  # for subsequent sites
  len = length(x)
  
  if(len > 1) {
    
    # likelihood given genotyping error
    ts1 <- 0.5 * dbetabinom(x[2:len], n[2:len], m=err , s=d)
    ts2 <- 0.5 * dbetabinom(n[2:len]- x[2:len], n[2:len], m=err, s=d)
    er2 <- (ge[2:len])*(ts1 + ts2)
    
    # consider all possible phasings with respect to the first SNP
    
    # precompute likelihoods of each subsequent SNP given 'in-phase' with first SNP
    snp.phase1.like <- (1-ge[2:len])*dbetabinom(x[2:len], n[2:len], m=0.5 + allelic.imbalance, s = d) + er2
    
    # precompute likelihoods of each subsequent SNP given 'out-of-phase' with first SNP
    snp.phase0.like <- (1-ge[2:len])*dbetabinom(x[2:len], n[2:len], m=0.5 - allelic.imbalance, s = d) + er2
    
    # create phase array
    phase1.like.array <- rep(NA, len)
    phase0.like.array <- rep(NA, len)
    
    # add likelihood for first site
    phase1.like.array[1] <- p1
    phase0.like.array[1] <- p1
    
    for(i in 2:len) {
      # prior SNP was either in-phase or out-of-phase with first SNP, consider
      # both mutually exclusive possibilities when computing combined likelihood of
      # all possible combinations of phases
      prev <- (0.5 * phase1.like.array[i-1]) + (0.5 * phase0.like.array[i-1])
      
      phase1.like.array[i] <- prev * snp.phase1.like[i-1]
      
      phase0.like.array[i] <- prev * snp.phase0.like[i-1]
    }
    
    # total likelihood is sum of last two elements
    l = -log(0.5*phase1.like.array[len] + 0.5*phase0.like.array[len])
  } else {
    l = -sum(log(p1))
  }
  
  return(l)
}


# get genes
genes= cnt$gene_id %>% unique %>% as.character

# function to run allele imbalance measurements
fun = function(i){
  
  # subset
  df= cnt[cnt$gene_id %in% genes[i],]
  
  # order
  df= df[order(df$start),]
  
  # number of het sites
  n.snp= nrow(df)
  
  # optimize
  m1 <- optimize(ll.fun, c(-0.5, 0.5), x= df$ref.matches, n= df$N, ge= df$genotype.error, err= err, d= d)
  
  # estimates of allelic.imbalance
  alt.ll <- m1$objective
  
  estimate <- m1$minimum
  
  # NULL hypothesis
  null.ll= ll.fun(par = 0, x= df$ref.matches, n = df$N, ge = df$genotype.error, err = err, d = d)
  
  # Likelihood ratio test
  lrt.stat <- 2 * (null.ll - alt.ll)
  
  pval <- pchisq(lrt.stat, df=1, lower.tail=F)
  
  result= data.frame(gene_id= genes[i], pval =pval, a= estimate)
  
  # return
  return(result)
  
}

my.list= llply(1:length(genes), fun, .progress = progress_text(char= "+"))

# combine list
res= do.call(rbind, my.list)

# correct p-values
res$fdr= p.adjust(res$pval, "fdr")

# order 
res= res[order(res$pval),]

# output results
write.csv(res, filename.output, row.names = F)



