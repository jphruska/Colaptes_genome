## TE window analyis script. 

setwd("/Volumes/HD-LXU3/Spring_2019_TTU_cluster/colaptes/gene_window_analysis/")

## Input txt file of headers from CDS that match scaffolds of interest. 

x <- read.table("cds_content_scaffolds_of_interest.txt", sep="\t", stringsAsFactors=F, fill=T, header=F)

## subtract cds enpoint from startpoint, and add difference as a column to data frame

x$V4 <- x$V3 - x$V2

# Turn off scientific notation for output with large numbers in sliding windows

options(scipen=999)


# sort, for each scaffold, the start and ends of the cds regions

x2 <-x %>% 
  arrange(V1,V3)

x2.colnames <- c("Scaffold", "CDS_start", "CDS_end", "length")

colnames(x2) <- x2.colnames


## Read in index file to get size of scaffolds 

require(GenomicRanges)


genome_index <- read.table("c05_colaptes_arcs3.FINAL.reordered_2.fasta.fai", sep="\t", stringsAsFactors=F)

genome_index <- genome_index[1:29,]

scaffold_sizes <- genome_index[,2]

scaffold_sizes

## Subset x to include only scaffolds of interest. In this case we are selecting all scaffolds that were syntenous with the T guttata genome, 
## but this can be modified. 

x.scaffolds.interest <- unique(x[,"V1"])[1:29]

x.1 <- x[x$query %in% x.scaffolds.interest, ]

## Determine window size 

window_size <- 100000

# Loop through each scaffold 
# output headers = scaffold, start, end, prop CDS 

output <- c()
num_windows <- c()
# For each scaffold
for(a in 1:length(x.scaffolds.interest)) {
  options(scipen=999)
  num_windows[a] <-floor(scaffold_sizes[a]/window_size)
  start <- 1 
  end <- window_size
  a_CDS <- x2[x2[,1] == x.scaffolds.interest[a],]
  # For each window
  for (b in 1:num_windows[a]){
    b_CDS <- a_CDS[a_CDS[,2] >= start & a_CDS[,2] <= end | a_CDS[,3] >= start & a_CDS[,3] <= end,]
    if(nrow(b_CDS) > 0) {
      b_CDS[b_CDS[,3] > end, 3] <- end
      b_CDS[b_CDS[,2] < start, 2] <- start
      b_seq <- c()
      for(c in 1:nrow(b_CDS)) {
        b_seq <- c(b_seq, b_CDS[c,2]:b_CDS[c,3])
      }
      b_CDS_count <- length(unique(b_seq)) / window_size
    } else {
      b_CDS_count <- 0
    }
    output <- rbind(output, c(x.scaffolds.interest[a], start, end, b_CDS_count))
    start <- start + window_size
    end <- end + window_size
  }
  print(x.scaffolds.interest[a]) 
}

output <- data.frame(Scaffold=as.character(output[,1]), Start=as.numeric(output[,2]), End=as.numeric(output[,3]), Prop_CDS = as.numeric(output[,4]))

write.table(output, "colaptes_cds_props.txt", sep = "\t", row.names = FALSE)


write.table(x2, "colaptes_cds_reordered.txt", sep = "\t", row.names = FALSE)
