## Verify the BCBP data files.

check.drp <- TRUE                       #Compute Density Recovery Profile?
if (check.drp) {
  library(sjedrp)
}



data.dir <- "~/mosaics/data/bcbp2"      #Where files are stored.


show.map <- function(field) {


  pch.cone <- 19; pch.bp <- 1           #plotting type of cone, bipolar
  
  filename  <- sprintf("%s/bc_bcbp_f%d.dat", data.dir, field)
  dat <- read.table(filename)

  cone.n <- length(which(dat[,3] == 1))
  bp.n    <- length(which(dat[,3] == 2))

  plot(dat[,1], dat[,2], asp=1, pch=ifelse(dat[,3]==1, pch.cone, pch.bp))


  bb.file <- sprintf("%s/bc_bcbp_f%d.w", data.dir, field)
  if (file.exists(bb.file)) {
    bb <- scan(bb.file)
    area <- (bb[2] - bb[1]) * (bb[4] - bb[3]) * 1e-6
    cat(sprintf("area in mm^2: %.3f\n", area))
    rect(bb[1], bb[3], bb[2], bb[4], lty=2)
    legend("topright",  c("BC", "BCBP"), pch=c(pch.cone, pch.bp))
    cat(sprintf("Number of BCs\t%d\n", cone.n))
    cat(sprintf("Number of BCBPs\t%d\n", bp.n))

    if (check.drp) {
      op <- par(mfrow=c(2,1))
      cones <- dat[which(dat[,3]==1),1:2]
      bp <-    dat[which(dat[,3]==2),1:2]

      nbins <- 10; r <- 10
      drp.c <- autodrp(cones[,1], cones[,2], nbins, r, bb)
      drp.b <- autodrp(bp[,1],       bp[,2], nbins, r, bb)
      plot(drp.c)
      plot(drp.b)
      par(op)
    }
    cat(sprintf("  BC effrad: %.1f um\n", drp.c$effrad))
    cat(sprintf("BCBP effrad: %.1f um\n", drp.b$effrad))
  }



}

show.map(1)
show.map(2)
show.map(3)
show.map(4)
show.map(5)
show.map(6)
