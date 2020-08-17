output_directory <- "cds_fasta_files"
output_summary <- "cds_summary_and_mapping.txt"
dir.create(output_directory)

# read in fasta files
colaptes_fasta <- scan("colaptes_cds_renamed.fasta", what="character")
colaptes_fasta <- gsub(">", "", colaptes_fasta)
apaloderma_fasta <- scan("apaloderma_cds_renamed.fasta", what="character")
apaloderma_fasta <- gsub(">", "", apaloderma_fasta)
buceros_fasta <- scan("buceros_cds_renamed.fasta", what="character")
buceros_fasta <- gsub(">", "", buceros_fasta)
merops_fasta <- scan("merops_cds_renamed.fasta", what="character")
merops_fasta <- gsub(">", "", merops_fasta)

# put each fasta into a dataframe with column 1 = name and column 2 = sequence
colaptes_cds <- data.frame(name=as.character(colaptes_fasta[1:(length(colaptes_fasta) / 2) * 2 - 1]),
                          sequence=as.character(colaptes_fasta[1:(length(colaptes_fasta) / 2) * 2]))
apaloderma_cds <- data.frame(name=as.character(apaloderma_fasta[1:(length(apaloderma_fasta) / 2) * 2 - 1]),
                           sequence=as.character(apaloderma_fasta[1:(length(apaloderma_fasta) / 2) * 2]))
buceros_cds <- data.frame(name=as.character(buceros_fasta[1:(length(buceros_fasta) / 2) * 2 - 1]),
                        sequence=as.character(buceros_fasta[1:(length(buceros_fasta) / 2) * 2]))
merops_cds <- data.frame(name=as.character(merops_fasta[1:(length(merops_fasta) / 2) * 2 - 1]),
                              sequence=as.character(merops_fasta[1:(length(merops_fasta) / 2) * 2]))



# read in blast files
q_colaptes_s_apaloderma <- read.table("q_colaptes_s_apaloderma.blast", sep="\t", stringsAsFactors=F, quote="\"")
q_colaptes_s_buceros <- read.table("q_colaptes_s_buceros.blast", sep="\t", stringsAsFactors=F, quote="\"")
q_colaptes_s_merops <- read.table("q_colaptes_s_merops.blast", sep="\t", stringsAsFactors=F, quote="\"")
q_apaloderma_s_colaptes <- read.table("q_apaloderma_s_colaptes.blast", sep="\t", stringsAsFactors=F, quote="\"")
q_buceros_s_colaptes <- read.table("q_buceros_s_colaptes.blast", sep="\t", stringsAsFactors=F, quote="\"")
q_merops_s_colaptes <- read.table("q_merops_s_colaptes.blast", sep="\t", stringsAsFactors=F, quote="\"")


# write header of summary file
write(c("colaptes_gene_number", "apaloderma_gene_name", "buceros_gene_name", "merops_gene_name", "number"), file=output_summary, ncolumns=5)

# list each colaptes predicted gene
colaptes_genes <- as.character(unique(colaptes_cds[,1]))

# start counter
counter <- 1

# loop for each predicted gene
for(a in 1:length(colaptes_genes)) {
  
  # matches where the subject is the other species
  s_apaloderma <- q_colaptes_s_apaloderma[grep(colaptes_genes[a], q_colaptes_s_apaloderma[,1]),]
  s_buceros <- q_colaptes_s_buceros[grep(colaptes_genes[a], q_colaptes_s_buceros[,1]),]
  s_merops <- q_colaptes_s_merops[grep(colaptes_genes[a], q_colaptes_s_merops[,1]),]
  # keep going if there is a match for all three species
  if(nrow(s_apaloderma) > 0 & nrow(s_buceros) > 0 & nrow(s_merops) > 0) {
    # check how many unique gene matches for each search
    if(length(unique(s_apaloderma[,2])) == 1) {
      s_apaloderma <- as.character(s_apaloderma[1,2])
    } else { # choose best match by minimum e-value, if two have same e-value, choose higher bit score
      s_apaloderma <- s_apaloderma[s_apaloderma[,11] == min(s_apaloderma[,11]),]
      s_apaloderma <- s_apaloderma[s_apaloderma[,12] == max(s_apaloderma[,12]),]
      s_apaloderma <- as.character(s_apaloderma[1,2])
    }
    if(length(unique(s_buceros[,2])) == 1) {
      s_buceros <- as.character(s_buceros[1,2])
    } else { # choose best match by minimum e-value, if two have same e-value, choose higher bit score
      s_buceros <- s_buceros[s_buceros[,11] == min(s_buceros[,11]),]
      s_buceros <- s_buceros[s_buceros[,12] == max(s_buceros[,12]),]
      s_buceros <- as.character(s_buceros[1,2])
    }
    if(length(unique(s_merops[,2])) == 1) {
      s_merops <- as.character(s_merops[1,2])
    } else { # choose best match by minimum e-value, if two have same e-value, choose higher bit score
      s_merops <- s_merops[s_merops[,11] == min(s_merops[,11]),]
      s_merops <- s_merops[s_merops[,12] == max(s_merops[,12]),]
      s_merops <- as.character(s_merops[1,2])
    }
    
    # search for the matched genes in the searches where colaptes was the subject to check for reciprocal matches
    q_apaloderma <- q_apaloderma_s_colaptes[grep(s_apaloderma, q_apaloderma_s_colaptes[,1], fixed = TRUE),]
    q_buceros <- q_buceros_s_colaptes[grep(s_buceros, q_buceros_s_colaptes[,1], fixed = TRUE),]
    q_merops <- q_merops_s_colaptes[grep(s_merops, q_merops_s_colaptes[,1], fixed = TRUE),]
    # keep going if there is a match for all three species
    if(nrow(q_apaloderma) > 0 & nrow(q_buceros) > 0 & nrow(q_merops) > 0) {
      if(length(unique(q_apaloderma[,2])) == 1) {
        q_apaloderma <- as.character(q_apaloderma[1,2])
      } else { # choose best match by minimum e-value, if two have same e-value, choose higher bit score
        q_apaloderma <- q_apaloderma[q_apaloderma[,11] == min(q_apaloderma[,11]),]
        q_apaloderma <- q_apaloderma[q_apaloderma[,12] == max(q_apaloderma[,12]),]
        q_apaloderma <- as.character(q_apaloderma[1,2])
      }
      if(length(unique(q_buceros[,2])) == 1) {
        q_buceros <- as.character(q_buceros[1,2])
      } else { # choose best match by minimum e-value, if two have same e-value, choose higher bit score
        q_buceros <- q_buceros[q_buceros[,11] == min(q_buceros[,11]),]
        q_buceros <- q_buceros[q_buceros[,12] == max(q_buceros[,12]),]
        q_buceros <- as.character(q_buceros[1,2])
      }
      if(length(unique(q_merops[,2])) == 1) {
        q_merops <- as.character(q_merops[1,2])
      } else { # choose best match by minimum e-value, if two have same e-value, choose higher bit score
        q_merops <- q_merops[q_merops[,11] == min(q_merops[,11]),]
        q_merops <- q_merops[q_merops[,12] == max(q_merops[,12]),]
        q_merops <- as.character(q_merops[1,2])
      }
      
      # keep going if all three best matches 
      if(q_buceros == colaptes_genes[a] & q_merops == colaptes_genes[a] & q_apaloderma == colaptes_genes[a]) {
        # write info to summary file
        write(c(colaptes_genes[a], s_apaloderma, s_buceros, s_merops, counter), file=output_summary, ncolumns=5, append=T)
        
        # write the fasta files for each gene that got this far
        write(">colaptes", file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1)
        write(as.character(colaptes_cds[colaptes_cds[,1] == colaptes_genes[a],2]), 
              file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)
        write(">apaloderma", file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)
        write(as.character(apaloderma_cds[apaloderma_cds[,1] == s_apaloderma,2]), 
              file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)
        write(">buceros", file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)
        write(as.character(buceros_cds[buceros_cds[,1] == s_buceros,2]), 
              file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)
        write(">merops", file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)
        write(as.character(merops_cds[merops_cds[,1] == s_merops,2]), 
              file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)	
        
        # add to counter
        counter <- counter + 1
      }
    }	
  }
}

