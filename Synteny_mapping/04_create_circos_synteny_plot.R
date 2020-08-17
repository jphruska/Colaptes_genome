library(OmicCircos)
options(stringsAsFactors=FALSE)
library(pals)
library(dplyr)

x <- read.table("colaptes_reordered_filtered.txt", sep="\t", header=T, row.names = NULL)

names(x) <- c("X.S1.", "X.E1.", "X.S2.", "X.E2.", "X.LEN.1.", "X.LEN.2." , 
              "X...IDY.", "X.LEN.R.", "X.LEN.Q.", "X.TAGS.", "X")

# subset scaffolds to largest 45 (Hi_C scaffold 1-45)
# generate list of scaffolds from 1-30 with for loop

hi_c_scaffold_1_40 <- c()

for (a in 1:40) {
  hi_c_scaffold_1_40[a] <- paste("Hi_C_scaffold_", a , sep = "")
  print(a)
}

# subset data frame so only largest 30 colaptes scaffolds are shown

x.1 <- filter(x, X %in% c("Ca_s1__1", "Ca_s3__1A", "Ca_s11__1B", "Ca_s17__1C", "Ca_s30__1D", 
                          "Ca_s5__2", "Ca_s13__2A", "Ca_s15__2B", "Ca_s24__2C", 
                          "Ca_s4__3", "Ca_s10__3A", "Ca_s28__3B", 
                          "Ca_s9__4", "Ca_s20__4A", "Ca_s27__4B", 
                          "Ca_s6__5", "Ca_s26__5A", "Ca_s8__6", "Ca_s7__7", 
                          "Ca_s12__8", "Ca_s16__9", "Ca_s19__10", "Ca_s22__11", 
                          "Ca_s18__12", "Ca_s21__13", "Ca_s25__14", "Ca_s23__18", 
                          "Ca_s29__20", "Ca_s2__Z"))

# subset data frame so only chromosomes 1-30 (including z) of gallus are shown

x.1 <- filter(x.1, X.TAGS. %in% c("Ga_1", "Ga_2", "Ga_3", "Ga_4", "Ga_5", "Ga_6", "Ga_7", 
                                  "Ga_8", "Ga_9", 
                                  "Ga_10", "Ga_11", "Ga_12", "Ga_13", 
                                  "Ga_14", 
                                  "Ga_18", "Ga_20",  
                                  "Ga_Z"))

# rearrange the input file to the following order
# target name, start, end, query name, start, end, identity, length target, length query
x.1 <- x.1[,c(10,1,2,11,3,4,7, 9,10)]

# subset scaffold names and start sites
x.2 <- x.1[,c(1,2,4,5)]

# rename colaptes scaffolds
x.2[,3] <- paste("Co_a", sapply(strsplit(x.2[,3], "_"), "[[", 3), sep="")


# bin the segments of the circos plot
segf_lengths <- c()
x3 <- c()
for(a in 1:length(unique(x.2[,1]))) {
  a_rep <- x.2[x.2[,1] == unique(x.2[,1])[a], ]
  test <- unique(sort(a_rep[,2]))
  a_rep[,2] <- match(a_rep[,2], test)
  x3 <- rbind(x3, a_rep)
  
  segf_lengths <- rbind(segf_lengths, c(unique(x.2[,1])[a], length(test)))
}


x4 <- c()
for(a in 1:length(unique(x.2[,3]))) {
  a_rep <- x3[x3[,3] == unique(x.2[,3])[a], ]
  test <- unique(sort(a_rep[,4]))
  a_rep[,4] <- match(a_rep[,4], test)
  x4 <- rbind(x4, a_rep)
  segf_lengths <- rbind(segf_lengths, c(unique(x.2[,3])[a], length(test)))
}
segf_lengths <- data.frame(chrom=as.character(segf_lengths[,1]), bins=as.numeric(segf_lengths[,2]))


# create segments object
seg_f <- c()
for(a in 1:nrow(segf_lengths)) {
  a_name <- segf_lengths[a,1]
  a_lengths <- 0:(segf_lengths[a,2] - 1)
  a_lengths2 <- a_lengths + 1
  a_name <- rep(a_name, length(a_lengths))
  a_v <- rep("NA", length(a_lengths))
  a_note <- rep("NA", length(a_lengths))
  a_output <- cbind(a_name, a_lengths, a_lengths2, a_v, a_note)
  seg_f <- rbind(seg_f, a_output)
}
seg_f <- data.frame(seg.name=as.character(seg_f[,1]), seg.Start=as.numeric(seg_f[,2]), seg.End=as.numeric(seg_f[,3]), the.v=as.character(seg_f[,4]), Note=as.character(seg_f[,5]))


# make the links object in the right format
link_names <- paste("n", seq(from=1, to=nrow(x4), by=1), sep="")
link_v <- data.frame(seg1=as.character(x4[,1]), po1=as.numeric(x4[,2]), name1=as.character(link_names), seg2=as.character(x4[,3]), po2=as.numeric(x4[,4]), name2=as.character(link_names))
# subset for easier plotting
link_subset <- link_v[sample(1:nrow(link_v), floor(nrow(link_v) / 100)), ]


# make the angular database
seg_names <- sort(unique(segf_lengths[,1]))

# modify for particular plot if necessary:
seg_names <- c("Ca_s1__1", "Ca_s3__1A", "Ca_s11__1B", "Ca_s17__1C", "Ca_s30__1D", 
               "Ca_s5__2", "Ca_s13__2A", "Ca_s15__2B", "Ca_s24__2C", 
               "Ca_s4__3", "Ca_s10__3A", "Ca_s28__3B", 
               "Ca_s9__4", "Ca_s20__4A", "Ca_s27__4B", 
               "Ca_s6__5", "Ca_s26__5A", "Ca_s8__6", "Ca_s7__7", 
               "Ca_s12__8", "Ca_s16__9", "Ca_s19__10", "Ca_s22__11", 
               "Ca_s18__12", "Ca_s21__13", "Ca_s25__14", "Ca_s23__18", 
               "Ca_s29__20", "Ca_s2__Z", "Ga_Z", "Ga_20", 
               "Ga_18", "Ga_14", "Ga_13", "Ga_12", "Ga_11", 
               "Ga_10", "Ga_9", "Ga_8", "Ga_7", "Ga_6", "Ga_5", "Ga_4", "Ga_3", "Ga_2", "Ga_1") 

seg_colors <- c(rep("darkorange", 29), rep("navyblue", 17))

# maybe for plotting opposite?
#seg_names <- rev(seg_names)
#seg_colors <- rev(seg_colors)

db <- segAnglePo(seg_f, seg=seg_names)



# colors for plotting
cols <- polychrome(length(unique(link_subset[,4])))
cols_rgb <- col2rgb(cols)
new_cols <- c()
for(a in 1:ncol(cols_rgb)) {
  a_col <- rgb(cols_rgb[1,a], cols_rgb[2,a], cols_rgb[3,a], max=255, alpha=0.4*255)
  new_cols <- c(new_cols, a_col)
}
cols <- new_cols

plot_colors <- rep("un", nrow(link_subset))
for(a in 1:length(unique(link_subset[,4]))) {
  a_rep <- sort(unique(link_subset[,4]))[a]
  plot_colors[link_subset[,4] == a_rep] <- cols[a]
}



par(mar=c(2,2,2,2)) 
plot(c(1,800),  c(1,800), type="n", axes=FALSE, xlab="XLAB", ylab="YLAB", main="")
circos(R=360, cir=db, type="chr", col=seg_colors, print.chr.lab=TRUE, W=40, scale = F)
circos(R=350, cir=db, W=40, mapping=link_subset, type="link", lwd=0.75, col=plot_colors)
