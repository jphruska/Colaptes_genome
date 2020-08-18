output_name <- "codeml_results_pvals_uncorrected.txt"
output_name2 <- "codeml_results_pvals_corrected.txt"
write(c("gene_number", "colaptes_un_p", "apaloderma_un_p", "buceros_un_p", "merops_un_p"), file=output_name, ncolumns=5, sep="\t")

x <- list.files(pattern="*fasta")
x2 <- list.files(pattern="*txt")
x_numbers <- as.numeric(sapply(strsplit(x, "_"), "[[", 1))

require(stats)
require(Biostrings)
require(seqinr)
require(rphast)

for(a in 1:max(x_numbers)) {
  x_match <- match(paste(a, "_aligned_trimmed.fasta", sep=""), x) # see if fasta file exists
  if(!is.na(x_match)) { # if yes read it
    a_rep <- readDNAStringSet(paste(a, "_aligned_trimmed.fasta", sep=""))
    if(a_rep@ranges@width[1] >= 150) { # only use those that are at least 50 AAs (150 nucleotides)
      x_match <- match(paste(a, "_total_output.txt", sep=""), x2)
      if(!is.na(x_match)) { # read in codeml output
        a_results <- read.table(paste(a, "_total_output.txt", sep=""), fill = T, stringsAsFactors=F)
        a_null_lnl <- a_results[1,5] # get the null model lnL
        a_alt_lnl <- a_results[c(3,5,7,9), 5] # get the alt models' lnL
        a_LRT <- 2 * (a_alt_lnl - a_null_lnl) # calculate the LRT
        a_uncorrected_p <- pchisq(a_LRT, df=1, lower.tail=FALSE) # get p-values for the LRT values (chi-square two tail)
        a_output <- c(a, a_uncorrected_p)
        write(a_output, file=output_name, ncolumns=5, append=T, sep="\t")
      }
    }
  } 
}

# read in previous output with p-values
output <- read.table(output_name, sep="\t", stringsAsFactors=F, header=T)
# calculate number of tests = number of genes * four tests
number_comparisons <- nrow(output) * 4
# multiple testing correction of the p-values using Benjamini & Hochberg (1995) (fdr)
output[,2] <- p.adjust(output[,2], method="fdr", n=number_comparisons)
output[,3] <- p.adjust(output[,3], method="fdr", n=number_comparisons)
output[,4] <- p.adjust(output[,4], method="fdr", n=number_comparisons)
output[,5] <- p.adjust(output[,5], method="fdr", n=number_comparisons)


# find minimum p-value for each gene and append that column to the output
min_p <- apply(output[,2:5], 1, min)
output <- cbind(output, min_p)
plot(min_p, pch=19, cex=0.1)

write.table(output, file=output_name2, sep="\t", quote=F, col.names=T, row.names=F)



# make the 4 fold degenerate sites output directory
dir.create("_4d_output")

# remove all significant tests and any rows missing info
filtered_output <- na.omit(output)
filtered_output <- filtered_output[filtered_output$min_p > 0.05, ]

# loop to
# read in multiple sequence alignments that are not under selection so as to get the four-fold degenerate sites
# for a later phylogeny
for(a in 1:nrow(filtered_output)) {
  a_rep <- read.msa(paste(filtered_output[a,1], "_aligned_trimmed.fasta", sep=""))
  a_feat <- feat(seqname="colaptes", feature="CDS", start=1, end=ncol(a_rep))
  
  a_4d_rep <- get4d.msa(a_rep, a_feat)
  write.msa(a_4d_rep, file=paste("_4d_output/", filtered_output[a,1], "_4d.fasta", sep=""), format="FASTA")
}



# list all the 4d alignments output
x_files <- list.files("_4d_output", full.names=T)
# loop to read in all alignments and concatenate
colaptes <- list()
apaloderma <- list()
buceros <- list()
merops <- list()
for(a in 1:length(x_files)) {
  a_rep <- readDNAStringSet(x_files[a])
  colaptes[[a]] <- as.character(a_rep)[1]
  apaloderma[[a]] <- as.character(a_rep)[2]
  buceros[[a]] <- as.character(a_rep)[3]
  merops[[a]] <- as.character(a_rep)[4]
}
colaptes <- paste(unlist(colaptes), collapse="")
apaloderma <- paste(unlist(apaloderma), collapse="")
buceros <- paste(unlist(buceros), collapse="")
merops <- paste(unlist(merops), collapse="")

output_name <- "_total_4d_sites.fasta"
write(">colaptes", file=output_name, ncolumns=1)
write(colaptes, file=output_name, ncolumns=1, append=T)
write(">apaloderma", file=output_name, ncolumns=1, append=T)
write(apaloderma, file=output_name, ncolumns=1, append=T)
write(">buceros", file=output_name, ncolumns=1, append=T)
write(buceros, file=output_name, ncolumns=1, append=T)
write(">merops", file=output_name, ncolumns=1, append=T)
write(merops, file=output_name, ncolumns=1, append=T)
