args = commandArgs(trailingOnly = TRUE)
filename.count= args[1]
filename.region= args[2]
filename.output= args[3]
print(args)

# load libraries
suppressMessages(library(data.table))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(rmutil))
suppressMessages(library(rtracklayer))
suppressMessages(library(pbapply))
suppressMessages(library(parallel))

# select colnames
col.names= c("seqnames", "start", "ref.matches", "alt.matches", "ref", "alt")

# edit counts
cnt= read.table(filename.count, header=T, skip = 87)

colnames(cnt) = col.names

cnt$N= rowSums(cnt[, c("ref.matches", "alt.matches")])

# filter counts
cnt = subset(cnt, N >= 10)

# add end position and name
cnt$end= cnt$start

# function to estimate dispersion under assumption of no allele-imbalance
mtmp <- function(par, x, n){
  
  # likelihood given het sites 
  p= dbetabinom(x, n, m= 0.5, s= par)
  
  # maximize likelihood for being het
  -sum(log(p))
}

# estimate dispersion
print("estimating dispersion")

m0<- optimise(mtmp, c(1e-05, 100), x= cnt$ref.matches, n= cnt$N)

d<- m0$minimum

# convert to granges
gr1= cnt %>% makeGRangesFromDataFrame(keep.extra.columns = T)

# test regions
reg= import.bed(filename.region)

reg$name= paste0(seqnames(reg), "_" , start(reg), "_" , end(reg))

# function to estimate allele proportions
ll.new<- function(par, x, n, d){
  
  allelic.imbalance <- par
  
  # for first site
  p1= dbetabinom(x[1], n[1], m= 0.5 + allelic.imbalance, s= d)
  
  # for subsequent sites
  len = length(x)
  
  if(len > 1) {
    
    # precompute likelihoods of each subsequent SNP given 'in-phase' with first SNP
    snp.phase1.like <- dbetabinom(x[2:len], n[2:len], m=0.5 + allelic.imbalance, s = d)
    
    # precompute likelihoods of each subsequent SNP given 'out-of-phase' with first SNP
    snp.phase0.like <- dbetabinom(x[2:len], n[2:len], m=0.5 - allelic.imbalance, s = d)
    
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
  
  # return
  return(l)
}

# function to estimate allele imbalance

fun = function(i){
  
  # test region
  gr2= reg[i]
  
  # subset by overlaps
  dat= subsetByOverlaps(gr1, gr2, type = "any") %>% as.data.frame()
  
  if(nrow(dat)> 3){
    
    # sort by start
    dat=dat[order(dat$start),]
    
    rownames(dat)<- NULL
    
    m1<- optimize(ll.new, c(-0.49, 0.49), x= dat$ref.matches, n= dat$N, d= d)
    
    # result
    res= data.frame(seqnames= seqnames(gr2),start= start(gr2), end= end(gr2), name= gr2$name, a =  m1$minimum ,nsites= nrow(dat))
    
  } else{
    
    res<- NULL
    
  }
  res
}

my.list=llply(1:length(reg), fun, .progress = progress_text(char="+"))

res= rbindlist(my.list) %>% as.data.frame()

write.csv(res , filename.output, row.names = F)

