## TE window analyis script. 

setwd("/Volumes/HD-LXU3/Spring_2019_TTU_cluster/colaptes/tes_hi_c_repeatmasker")

## Input repeatmasker output that has had duplicate TE sequences removed. 

x <- read.table("c05_colaptes_arcs3.FINAL.reordered_2.fasta.no_dups.out", sep="", stringsAsFactors=F, fill=T, header=T)

# Turn off scientific notation for output with large numbers in sliding windows

options(scipen=999)

## Read in index file to get size of scaffolds 

require(GenomicRanges)

genome_index <- read.table("c05_colaptes_arcs3.FINAL.reordered_2.fasta.fai", sep="\t", stringsAsFactors=F)

genome_index <- genome_index[1:29,]

scaffold_sizes <- genome_index[,2]

scaffold_sizes

## Subset x to include only scaffolds of interest. In this case we are selecting all scaffolds that were syntenous with the T guttata genome, 
## but this can be modified. 

x.scaffolds.interest <- unique(x[,"query"])[1:29]

x.1 <- x[x$query %in% x.scaffolds.interest, ]

## Determine window size 

window_size <- 100000

# Loop through each scaffold 
# output headers = scaffold, start, end, prop TEs, prop CR1, prop ERV 

output <- c()
num_windows <- c()
# For each scaffold
for(a in 1:length(x.scaffolds.interest)) {
  options(scipen=999)
  num_windows[a] <-floor(scaffold_sizes[a]/window_size)
  start <- 1 
  end <- window_size
  a_TEs <- x.1[x.1[,5] == x.scaffolds.interest[a],]
  # For each window
  for (b in 1:num_windows[a]){
    b_TEs <- a_TEs[a_TEs[,6] >= start & a_TEs[,6] <= end | a_TEs[,7] >= start & a_TEs[,7] <= end,]
    if(nrow(b_TEs) > 0) {
      b_TEs[b_TEs[,7] > end, 7] <- end
      b_TEs[b_TEs[,6] < start, 6] <- start
      b_seq <- c()
      for(c in 1:nrow(b_TEs)) {
        b_seq <- c(b_seq, b_TEs[c,6]:b_TEs[c,7])
      }
      b_TE_count <- length(unique(b_seq)) / window_size
    } else {
      b_TE_count <- 0
    }
    output <- rbind(output, c(x.scaffolds.interest[a], start, end, b_TE_count))
    start <- start + window_size
    end <- end + window_size
  }
  print(x.scaffolds.interest[a]) 
}

output <- data.frame(Scaffold=as.character(output[,1]), Start=as.numeric(output[,2]), End=as.numeric(output[,3]), Prop_TE = as.numeric(output[,4]))

write.table(output, "colaptes_te_props.txt", sep = "\t", row.names = FALSE)



