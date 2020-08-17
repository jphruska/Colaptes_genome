x <- read.table("colaptes_repeatmodeler_genomic2.blast", sep="\t", stringsAsFactors=F)

# read in reference index to make sure bed file doesn't overlap edges
fai <- read.table("c05_colaptes_arcs3.fasta.fai", sep="\t", stringsAsFactors=F)

# minimum number of hits to keep element
min_hits <- 10

# maximum hits to keep 
max_hits <- 50

# define flank length to extract
flank_length <- 1000

# output directory
out_dir <- "colaptes_bed_extract"
dir.create(out_dir)

# remove hits less than 100 bp
x <- x[x[,7] >= 100,] 


# loop for each element
output <- c()
for(a in 1:length(unique(x[,1]))) {
  a_rep <- x[x[,1] == unique(x[,1])[a], ]
  # only continue if min_hits is met
  if(nrow(a_rep) >= min_hits) {
    a_rep <- a_rep[order(a_rep[,8]), ]
    if(nrow(a_rep) > max_hits) {
      a_rep <- a_rep[1:max_hits,]
    }
    a_name <- strsplit(a_rep[,1], "//")[[1]][1]
    
    # extract the "bed" part of the blast output and add the flanking length
    a_rep2 <- a_rep[,3:5]
    # sort the bed file start and end positions
    for(b in 1:nrow(a_rep2)) {
      b_rep <- a_rep2[b,]
      if(b_rep[1,2] > b_rep[1,3]) {
        b_rep1 <- b_rep[1,3]
        b_rep2 <- b_rep[1,2]
        a_rep2[b,2] <- b_rep1
        a_rep2[b,3] <- b_rep2
      }
    }
    a_rep2[,2] <- a_rep2[,2] - flank_length - 1 # minus the extra one for the zero-based bedtools
    a_rep2[,3] <- a_rep2[,3] + flank_length
    
    # loop to make sure that none of the new bed file is beyond the edges of scaffolds
    for(b in 1:nrow(a_rep2)) {
      b_fai <- fai[fai[,1] == a_rep2[b,1],]
      if(a_rep2[b,2] < 0) {
        a_rep2[b,2] <- 0
      }
      if(a_rep2[b,3] > b_fai[1,2]) {
        a_rep2[b,3] <- b_fai[1,2]
      }
    }
    write.table(a_rep2, file=paste(out_dir, "/", a_name, sep=""), sep="\t", quote=F, row.names=F, col.names=F)
  }	
}
