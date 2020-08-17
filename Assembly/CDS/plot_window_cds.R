options(scipen=999)
require(stats)

x <- 	read.table("colaptes_cds_props.txt", sep="\t", header=T, stringsAsFactors=F)
total_windows <- nrow(x)
x <- na.omit(x)

window_size <- 10


# what are the unique chromosomes and their bounding areas for plotting?
chr <- unique(x[,1])
chr_polygons_cds <- list()
# make the plotting polygons
for(a in 1:length(chr)) {
	a1 <- rownames(x)[x[,1] == chr[a]]
	a2 <- a1[length(a1)]
	a1 <- a1[1]
	chr_polygons_cds[[a]] <- rbind(c(a1, 0), c(a2, 0), c(a2, 0.8), c(a1, 0.8), c(a1, 0))
}

# set up plotting dimensions
par(mfrow=c(1,1))
par(mar=c(0.5,5,1,0))

# plot cds
plot(c(-1,-1), ylim=c(0,0.2), xlim=c(1, total_windows), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="Prop. CDS")
odd <- 0
for(a in 1:length(chr_polygons_cds)) {
	if(odd == 1) {
		polygon(chr_polygons_cds[[a]], col="snow2", border="white")
		odd <- 0	
	} else {
		odd <- 1
	}
}
# plot
points(rownames(x), x$Prop_CDS, pch=19, cex=0.1, col="azure3")	

# sliding windows
total_rep <- c()
place_rep <- c()
sliding_windows <- ceiling(total_windows / window_size)
for(b in 0:sliding_windows) {
	b_rep <- seq(from=(b*window_size - (window_size/2 -1)), to=(b*window_size + (window_size/2)), by=1)
	b_rep <- b_rep[b_rep >= 1 & b_rep <= total_windows]
	b_rep <- b_rep[x[b_rep,1] %in% x[b*window_size,1]]
	b_rep <- na.omit(match(b_rep, rownames(x)))
	total_rep <- c(total_rep, mean(x$Prop_CDS[b_rep]))
	place_rep <- c(place_rep, b*window_size)
}
lines(place_rep, total_rep, lwd=1.2, col="navyblue")

