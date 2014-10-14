#!/usr/bin/Rscript

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'help', 'h', 0, "logical",
  'threads', 't', 1, "integer",
  'files', 'f', 1, "character",
  'output', 'o', 1, "character"
), byrow=TRUE, ncol=4);
library(getopt)
opt = getopt(spec);
# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}

if ( is.null(opt$threads)) { opt$threads = 1 }

library(reshape2)
library(EBSeq)

files <- read.table(opt$files, header=FALSE, sep=",", as.is=TRUE)

names(files) <- c("name", "type", "rep", "file")

data <- data.frame()

for (i in 1:nrow(files)) {
  name <- files[i,1]
  type <- files[i,2]
  rep <- files[i,3]
  file <- files[i,4]
  tmp <- read.table(file, sep="\t",header=TRUE, as.is=TRUE)
  tmp <- tmp[,c("target_id", "eff_counts")]
  tmp$eff_counts <- as.integer(tmp$eff_counts)
  tmp$variable <- paste(name, type, rep, sep="-")
  data <- rbind(data, tmp)
}

# unmelt data. acast returns matrix
# target_id is used as row names
# variable becomes column header
# value.var = "eff_counts" : the column eff_counts is used as the values in the new dataframe
eff_counts_data <- dcast(data, target_id ~ variable, value.var = "eff_counts")

anno <- eff_counts_data[,c(1)]                 # get names
eff_counts_data <- eff_counts_data[,-c(1)]     # remove names column
eff_counts_mat <- as.matrix(eff_counts_data)   # convert to matrix
rownames(eff_counts_mat) <- anno               # set matrix row names

f_mat <- as.matrix(files)

sizes <- MedianNorm(eff_counts_mat)
#sizes=QuantileNorm(counts, 0.75) # upper-quantile normalisation

conditions <- as.factor(rep(c("BS", "M"),each=nrow(files)/2))

ebOut <- EBTest(Data=eff_counts_mat, Conditions=conditions, sizeFactors=sizes, maxround=6)

posteriorProbs <- GetPPMat(ebOut)

probs <- data.frame(posteriorProbs)
names(probs) <- c("PPEE", "PPDE")

write.table(probs, file=paste(opt$output,"/",name,"_DE_ebseq.txt",sep=""), sep="\t",quote=F,row.names=T)