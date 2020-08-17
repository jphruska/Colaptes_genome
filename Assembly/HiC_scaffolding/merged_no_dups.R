x_files <- list.files(pattern="split*")

for(a in 1:length(x_files)) {
a_rep <- read.table(x_files[a], sep=" ", stringsAsFactors=F)

# remove NA column on end
a_rep <- a_rep[,1:16]

# add the last line from the previous file (if this is not the first file) to check for dups
if(a != 1) { a_rep <- rbind(a_continue, a_rep) }

# remove identical alignments
a_rep1 <- paste(a_rep[,2], a_rep[,3], a_rep[,6], a_rep[,7])
a_rep <- a_rep[match(unique(a_rep1), a_rep1),]

# map quality >= 1
a_rep <- a_rep[a_rep[,9] >= 10 & a_rep[,12] >= 10, ]
# remove lines that have less than a 500 bp alignment distance on same contig
a_rep1 <- a_rep[(a_rep[,2] == a_rep[,6] & abs(a_rep[,3] - a_rep[,7]) < 100) == FALSE, ]

# remove NA column on end
a_rep1 <- a_rep1[,1:16]

# keep last line to check for matches in next file (since it is sorted)
a_continue <- a_rep1[nrow(a_rep1), ]

# remove the first line that was kept from the previous file
if(a != 1) { a_rep1 <- a_rep1[2:nrow(a_rep1), ] }

# create output file
write.table(a_rep1, file=paste("no_dupes_", a, ".txt", sep=""), sep=" ", quote=F, row.names=F, col.names=F)
}
