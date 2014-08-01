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
library(baySeq)

cl <- makeCluster(opt$threads)

files <- read.table(opt$files, header=FALSE, sep=",", as.is=TRUE)

names(files) <- c("name", "type", "rep", "file")
print(files)
print(opt$output)

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
replicates <- f_mat[,2]
nde <- rep(1,times = nrow(files))
de <- rep(c(1,2),times = 1, each=nrow(files)/2)
groups <- list(NDE = nde, DE = de)

CD <- new("countData", data = eff_counts_mat, replicates = replicates, groups = groups)

libsizes(CD) <- getLibsizes(CD)

CD <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = cl)
CD <- getLikelihoods.NB(CD, pET = 'BIC', cl = cl)

output_DE <- topCounts(CD, group = "DE", number = 100000, normaliseData=TRUE)

write.table(output_DE, file=paste(opt$output,"/",name,"_DE_bayseq.txt",sep=""), sep="\t", row.names=TRUE, quote=F)

