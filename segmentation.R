args = commandArgs(trailingOnly = TRUE)
filename.tumor= args[1]
filename.normal= args[2]
filename.output1= args[3]
filename.output2= args[4]

# load libraries
suppressMessages(library(data.table))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(DNAcopy))

# tumor
tm= fread(filename.tumor) %>% as.data.frame()

# normal
nr= fread(filename.normal) %>% as.data.frame()

# combine
df= merge(tm, nr, by= c("seqnames", "start", "end", "name"))

# change colnames
colnames(df)[5:8]= c("a_tumor", "nsites_tumor", "a_normal", "nsites_normal") 

# difference
df$delta_a= abs(df$a_tumor)- abs(df$a_normal)

write.csv(df, filename.output1)

# chr
chr= seq(1,22, by=1)

# perform binary segmentation
my.seg.fun= function(i){
  
  # subset
  x= df[df$seqnames %in% chr[i],]

  # create object
  CNA.object <- CNA(genomdat = x$delta_a, chrom = x$seqnames,
                    maploc = x$start,data.type="logratio",
                    sampleid=gsub(".csv", "",filename.tumor))
  
  # smooth
  smoothed.CNA.object <- smooth.CNA(CNA.object)
  
  # perform segmentation
  segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
  
  # parse obj
  seg= segment.smoothed.CNA.object$output
  
  seg$ID= gsub("\\.", "-", seg$ID)
  
  colnames(seg)<- c("sample", "seqnames","start", "end", "num.mark", "scna_score")
  
  seg
}

my.list= llply(1:length(chr), my.seg.fun, .progress = progress_text(char="+"))

# combine list
seg= rbindlist(my.list) %>% as.data.frame()

# output results
write.csv(seg, filename.output2, row.names = F)



