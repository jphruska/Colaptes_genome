# head
#	SW	perc	perc	perc	query	position	in	query	matching	repeat	position	in	repeat
#score	div.	del.	ins.	sequence	begin	end	(left)	repeat	class/family	begin	end	(left)	ID
#

require(GenomicRanges)

genome_index <- read.table("c05_colaptes_arcs3.FINAL.reordered_2.fasta.fai", sep="\t", stringsAsFactors=F)
genome_index <- genome_index[,1:2]
options(scipen=999)

# read repeatmasker file
x <- read.table("c05_colaptes_arcs3.FINAL.reordered_2.fasta.out", sep="", stringsAsFactors=F, fill=T)

# start output
output_name <- "c05_colaptes_arcs3.FINAL.reordered_2.fasta.no_dups.out"
write(c("SW_score", "perc_div", "perc_delete", "perc_insert", "query", "q_start", "q_end", "q_leftover", "orientation", "repeat",
        "class/family", "s_start", "s_end", "s_leftover", "ID", "overlap"), file=output_name, ncolumns=16, sep="\t")

for(a in 1:nrow(genome_index)) {
  print(a)
  x_rep <- x[x[,5] == genome_index[a,1], ]
  
  if(nrow(x_rep) > 0) {
    # make genomic ranges object
    x_temp <- x_rep[,5:7]
    x_temp <- data.frame(chr=as.character(x_temp[,1]), start=as.numeric(x_temp[,2]), end=as.numeric(x_temp[,3]))
    refGR <- makeGRangesFromDataFrame(x_temp)	
    testGR <- makeGRangesFromDataFrame(x_temp)
    # overlap finding
    overlaps <- findOverlaps(refGR, testGR)
    # convert to matrix and remove self matches
    potential_overlaps <- cbind(overlaps@from, overlaps@to)
    potential_overlaps <- potential_overlaps[potential_overlaps[,1] != potential_overlaps[,2], ]
    
    # go through each row of xrep and look at any overlaps to determine if parts need to be removed or not
    for(b in 1:nrow(x_rep)) {
      b_rep <- x_rep[b,]
      # possible higher matches
      b_rep2 <- x_rep[potential_overlaps[potential_overlaps[,1] %in% b,2],]
      if(nrow(b_rep2) > 0) {
        # keep only possible alternatives with 
        b_rep2 <- b_rep2[b_rep2[,2] < b_rep[1,2],]
        # if no better matches write, otherwise carry on
        if(nrow(b_rep2) > 0) {
          # get bp of element we are testing
          b_rep_bp <- seq(from=b_rep[1,6], to=b_rep[1,7], by=1)
          # get bp of possible better matches
          b_alt_bp <- list()
          for(d in 1:nrow(b_rep2)) {
            b_alt_bp[[d]] <- seq(from=b_rep2[d,6], to=b_rep2[d,7], by=1)
          }
          b_alt_bp <- unique(unlist(b_alt_bp))
          # find the bp in this rep that are not in the possible alternatives
          b_rep_bp <- b_rep_bp[b_rep_bp %in% b_alt_bp == F]
          # determine if the element has been split or not and write one or more intervals
          if(length(unique(diff(b_rep_bp))) == 1) { # only one interval
            if(unique(diff(b_rep_bp)) == 1) {
              b_interval <- c(min(b_rep_bp), max(b_rep_bp))
              # input new interval into the other information
              b_output <- c(as.character(b_rep[1,1:5]), as.character(b_interval), as.character(b_rep[1,8:ncol(b_rep)]))
              write(b_output, file=output_name, sep="\t", ncolumns=16, append=T)
            }
          } else if(length(unique(diff(b_rep_bp))) > 1) { # more than one interval
            # if multiple intervals, write multiple lines split into the intervals
            b_diff <- c(1,diff(b_rep_bp))
            b_diff <- seq(from=1, to=length(b_rep_bp), by=1)[b_diff != 1]
            b_diff <- sort(c(b_diff, b_diff - 1, 1, length(b_rep_bp)))
            # write output for each interval
            for(d in 1:(length(b_diff) / 2)) {
              d1 <- d * 2 - 1
              d2 <- d * 2
              b_interval <- c(b_rep_bp[b_diff[d1]], b_rep_bp[b_diff[d2]])
              # input new interval into the other information
              b_output <- c(as.character(b_rep[1,1:5]), as.character(b_interval), as.character(b_rep[1,8:ncol(b_rep)]))
              write(b_output, file=output_name, sep="\t", ncolumns=16, append=T)
            }
          }
        } else { # no better alternatives
          write(as.character(b_rep), file=output_name, sep="\t", ncolumns=16, append=T)
        }
      } else { # no overlapping alternatives
        write(as.character(b_rep), file=output_name, sep="\t", ncolumns=16, append=T)
      }	
    }
  }
}



