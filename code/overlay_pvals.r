# set root directory
AVLROOT <- "."

commandlineargs <- commandArgs(TRUE)

nfields <- 5
filenames <- vector(length=nfields)

# first nfields command line arguments are file names
for(i in 1:nfields) {
    filenames[i] <- commandlineargs[i]
}
# then the next is the parameter that was being varied
param <- commandlineargs[nfields+1]

#print(filenames)

# set x labels for different parameters
labels <- list(depth="Max depth of mobile BCBPs", rel_force_mean="Mean relative force")

# open pdf file
pdf(file=paste("pvals_", param, "_combined_fields.pdf", sep=""))

colours <- rainbow(6)

# get maximum and minimum x-values and y-values for plotting
n <- length(filenames)
tables <- list()  # initialise a list of tables for reuse later
minx <- NA
maxx <- NA
miny <- NA
maxy <- NA
for(j in 1:n) {
	file <- filenames[j]
        path <- sprintf("%s/%s", AVLROOT, file)
        #print(path)
	table <- read.table(path, header=TRUE)
        local_minx <- min(as.numeric(row.names(table)))
        if(is.na(minx) || local_minx < minx) {
            minx <- local_minx
        }
        local_maxx <- max(as.numeric(row.names(table)))
        if(is.na(maxx) || local_maxy > maxx) {
            maxx <- local_maxx
        }
        local_miny <- min(table$pBCBCBP)
        if(is.na(miny) || local_miny < miny) {
            miny <- local_miny
        }
        local_maxy <- max(table$pBCBCBP)
        if(is.na(maxy) || local_maxy > maxy) {
            maxy <- local_maxy
        }
        tables[[j]] <- table
}

#print(minx)
#print(maxx)
#print(miny)
#print(maxy)

first <- TRUE
for(k in 1:n) {
    table <- tables[[k]]
    x <- row.names(table)  # range of values iof 'param' for this field
    y <- table$pBCBCBP  # range of p values for this field
    if(!first) {
        par(new=TRUE)  # plot on same graph
        plot(x, y, type="l", lty=k%%6, lwd=3, col=colours[k%%6],
             xlab=labels[[param]], ylab="Cross p value",
             xlim=c(minx, maxx), ylim=c(miny, maxy))
    } else {
        plot(x, y, type="l", lty=k%%6, lwd=3, col=colours[k%%6],
             xlab=NULL, ylab=NULL, # stop r superimposing 'x' and 'y' on labels
             xlim=c(minx, maxx), ylim=c(miny, maxy),
             main=paste("Sensitivity to ", param, ", all fields", sep=""))
        first <- FALSE
    }
}

dev.off()
