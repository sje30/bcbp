###################################################################################################
######### R script
######### Author: Andreas Lønborg, avl26@cam.ac.uk
######### Last updated: 2010.06.26
###################################################################################################

## SJE: Set root directory for where data files are.
date()
##AVLROOT <- Sys.getenv("AVLROOT")
AVLROOT <- "."


#### Importing command line arguments ###
#### The format for command line arguments is: 'fieldnumber' 'mean BC exclusion zone' 'st BC exclusion zone'

commandlineargs <- as.numeric(commandArgs(TRUE))
field <- commandlineargs[1]

depth.low <- commandlineargs[2]
depth.high <- commandlineargs[3]
depth.incr <- commandlineargs[4]
loopsize <- commandlineargs[5]                # Number of iterations in simulation study
ulimit <- commandlineargs[6]			# Upper limit for calculation of p scores (mu meter).

depths <- seq(depth.low,depth.high,depth.incr)

pBC <- matrix(nrow=length(depths),ncol=1)
rownames(pBC) <- depths
pBCBP <- matrix(nrow=length(depths),ncol=1)
rownames(pBCBP) <- depths
pBCBCBP <- matrix(nrow=length(depths),ncol=1)
rownames(pBCBCBP) <- depths

#################################
### Assumptions
#################################

# The neurons in this study are located in mid-peripheral retina of the macaque
# The pedicles of BCs can be considered vertical immobile extensions of the BC cell bodies
# The pedicles of BCs have the same vertical position and the BC pedicles are located in the OPL, 10 mu meter from the OPL-INL interface (Koyama 1992, figure)
# The perikarya of BCBPs are located at 33+-11% from OPL to IPL (Koyama 1992)
# The INL is 40 mu meter thick in mid peripheral retina (Koyama 1992, figure 2A, left), remember to adjust this according to what field.
# BCBPs located more than 12 mu meter from the INL-OPL interface (ie. d > 22) cannot move because the are restrained from the layer above (12 = 3*BCBP radius). This is in accordance with Figure 2 in Koyama 1992.
# The field simulated here is intended to mimic the third data field, which has 119 BCs and 268 BCBP in a 633x579 mu meter field. Located 5.5 mumeter from fovea, temporal.
# The cones form a mosaic in order to cover the visual field homogenously, the BCBP are not governed by such principle, the are only governed my a small exclusion zone and then dendritic tension.
# The BC exclusion zone is drawn from a Gaussian distribution, to "soften" the transition from forbidden to allowed. This improves the fit with data.

#################################
### Notes
#################################

#Field 3: the BCBP parms could be set back to 16, 4.
# Fields with low BCBP density have more flat cdrp. This is intuitive according to my model. If there are less BCBP, then there will be less connections to BC and thus less scope for movement that can generate the peak.
# my model fits well with the thougts from Kouyama.

#################################
### Sourcing functions
#################################
## SJE:
extra.functions <- sprintf("%s/functions_mimicsthesis.R", AVLROOT)
source(extra.functions)


#################################
### Loading packages
#################################
library(stats)
library(class)
library(graphics)
library(splancs)
library(sjedrp)

##################################
### Initialising all list variables
#################################
QBCBPPos.new <- list()
Qa <- list() 
Qbin <- list() 
Qcountcrossall <- list() 
Qdata1 <- list() 
Qdendritic.length <- list() 
Qdendr.length <- list() 
QdexcBC <- list() 
QdexcBCBP <- list() 
QdexcBC.sd <- list() 
QdrpBCall <- list() 
QdrpBCBCBPall <- list() 
QdrpBCBCBPold <- list() 
QdrpBCBCBPsjeafter.density <- list() 
QdrpBCBCBPsjeafter.effrad <- list() 
QdrpBCBCBPsjeafter.maxr <- list() 
QdrpBCBCBPsjebefore <- list() 
QdrpBCBCBPsjebefore.ds <- list() 
QdrpBCBCBPsjebefore.k <- list() 
QdrpBCBCBPsjebefore.p <- list() 
QdrpBCBPall <- list() 
QdrpBCBPsjeafter <- list() 
QdrpBCBPsjeafter.ds <- list() 
QdrpBCBPsjeafter.k <- list() 
QdrpBCBPsjeafter.p <- list() 
QdrpBCBPsjebefore.density <- list() 
QdrpBCBPsjebefore.effrad <- list() 
QdrpBCBPsjebefore.maxr <- list() 
QdrpBCdatasje <- list() 
QdrpBCsje.density <- list() 
QdrpBCsje.effrad <- list() 
QdrpBCsje.maxr <- list() 
Qeffrad_BCBP_ttest <- list() 
QGBCBPafter <- list() 
QGBCBPbefore <- list() 
QGBCBPdata <- list() 
QGBCenvelope <- list() 
QINLdepth <- list() 
Qk <- list() 
QKBCall <- list() 
QKBCBCBPall <- list() 
QKBCBPall <- list() 
QKBCBPmean <- list() 
QKBCmean <- list() 
QLBCBCBPafter <- list() 
QLBCBCBPbefore <- list() 
QLBCBCBP.data <- list() 
QLBCBPafterenvelope <- list() 
QLBCBPbeforeenvelope <- list() 
QLBC.data <- list() 
Qlinks <- list() 
Qloopsize <- list() 
Qly <- list() 
Qmaxdendr <- list() 
Qmaximumy <- list() 
Qmaxr_BCBP_ttest <- list() 
Qmean.dendritic.length <- list() 
Qmean.drpBCBCBPsjeafter.density <- list() 
Qmean.drpBCBCBPsjeafter.effrad <- list() 
Qmean.drpBCBCBPsjeafter.maxr <- list() 
Qmean.drpBCBCBPsjebefore.density <- list() 
Qmean.drpBCBCBPsjebefore.effrad <- list() 
Qmean.drpBCBCBPsjebefore.maxr <- list() 
Qmean.drpBCBPsjeafter.density <- list() 
Qmean.drpBCBPsjeafter.effrad <- list() 
Qmean.drpBCBPsjeafter.maxr <- list() 
Qmean.drpBCBPsjebefore.density <- list() 
Qmean.drpBCBPsjebefore.effrad <- list() 
Qmean.drpBCBPsjebefore.maxr <- list() 
Qmean.drpBCsje.density <- list() 
Qmean.drpBCsje.effrad <- list() 
Qmean.drpBCsje.maxr <- list() 
Qmeanfield.dendritic.length <- list() 
Qmeanfield.migration.dist <- list() 
Qmean.prop.movedBCBPs <- list() 
Qmean.synapses <- list() 
Qminimumy <- list() 
QmovedBCBPs <- list() 
QnbBCBP <- list() 
QnbBCBP.data <- list() 
QnBC <- list() 
QnbDistBC <- list() 
QnbDistBCBPbefore <- list() 
QnbDistBC.data <- list() 
QnC <- list() 
QnonconnectingBCBP <- list() 
QpBC <- list() 
Qp_BCBCBP_ttest <- list() 
Qp_BCBP_ttest <- list() 
Qpolygon <- list() 
Qpos1.cross.all.data <- list() 
Qpos2.cross.all.data <- list() 
Qr <- list() 
Qreac_store <- list() 
Qrel.p <- list() 
QRI_BC <- list() 
QRI_BCBPbefore <- list() 
QRI_BCBP_ttest <- list() 
QRI_BC_ttest <- list() 
Qsd.dendritic.length.post <- list() 
Qsd.drpBCBCBPsjebefore.ds <- list() 
Qsd.drpBCBPsjebefore.ds <- list() 
Qsdfield.dendritic.length <- list() 
Qsdfield.migration.dist <- list() 
Qsd.prop.movedBCBPs <- list() 
Qsd.synapses <- list() 
Qsums <- list() 
QTBC <- list() 
QTBCBP <- list() 
Qulimit <- list() 
Qvert.d.mean <- list() 
Qvert.d.thr <- list() 
QBCBPPos <- list() 
QBCindex <- list() 
QconePos <- list() 
Qdata <- list() 
Qdata2 <- list() 
Qdendritic.length.post <- list() 
Qdendr.length.diff <- list() 
Qdendr.length.post <- list() 
QdexcBC.trunc <- list() 
QdexcBCBP.trunc <- list() 
QdexcBCBP.sd <- list() 
QdexcCone <- list() 
Qdist.cross.all.data <- list() 
QdrpBC <- list() 
QdrpBCBCBP <- list() 
QdrpBCBCBPdatasje <- list() 
QdrpBCBCBPsjeafter <- list() 
QdrpBCBCBPsjeafter.ds <- list() 
QdrpBCBCBPsjeafter.k <- list() 
QdrpBCBCBPsjeafter.p <- list() 
QdrpBCBCBPsjebefore.density <- list() 
QdrpBCBCBPsjebefore.effrad <- list() 
QdrpBCBCBPsjebefore.maxr <- list() 
QdrpBCBP <- list() 
QdrpBCBPdatasje <- list() 
QdrpBCBPsjeafter.density <- list() 
QdrpBCBPsjeafter.effrad <- list() 
QdrpBCBPsjeafter.maxr <- list() 
QdrpBCBPsjebefore <- list() 
QdrpBCBPsjebefore.ds <- list() 
QdrpBCBPsjebefore.k <- list() 
QdrpBCBPsjebefore.p <- list() 
QdrpBCsje <- list() 
QdrpBCsje.ds <- list() 
QdrpBCsje.k <- list() 
QdrpBCsje.p <- list() 
Qeffrad_BC_ttest <- list() 
QGBC <- list() 
QGBCBPafterenvelope <- list() 
QGBCBPbeforeenvelope <- list() 
QGBCdata <- list() 
Qi <- list() 
Qj <- list() 
QK <- list() 
QKBCBCBPafter <- list() 
QKBCBCBPbefore <- list() 
QKBCBCBPmean <- list() 
QKBCBP.data <- list() 
QKBC.data <- list() 
QKbefore <- list() 
QLBC <- list() 
QLBCBCBPafterenvelope <- list() 
QLBCBCBPbeforeenvelope <- list() 
QLBCBPafter <- list() 
QLBCBPbefore <- list() 
QLBCBP.data <- list() 
QLBCenvelope <- list() 
Qlinks.new <- list() 
Qlx <- list() 
Qmaximumx <- list() 
Qmaxr_BCBCBP_ttest <- list() 
Qmaxr_BC_ttest <- list() 
Qmean.dendritic.length.post <- list() 
Qmean.drpBCBCBPsjeafter.ds <- list() 
Qmean.drpBCBCBPsjeafter.k <- list() 
Qmean.drpBCBCBPsjeafter.p <- list() 
Qmean.drpBCBCBPsjebefore.ds <- list() 
Qmean.drpBCBCBPsjebefore.k <- list() 
Qmean.drpBCBCBPsjebefore.p <- list() 
Qmean.drpBCBPsjeafter.ds <- list() 
Qmean.drpBCBPsjeafter.k <- list() 
Qmean.drpBCBPsjeafter.p <- list() 
Qmean.drpBCBPsjebefore.ds <- list() 
Qmean.drpBCBPsjebefore.k <- list() 
Qmean.drpBCBPsjebefore.p <- list() 
Qmean.drpBCsje.ds <- list() 
Qmean.drpBCsje.k <- list() 
Qmean.drpBCsje.p <- list() 
Qmeanfield.dendritic.length.post <- list() 
Qmean.migration.dist <- list() 
Qmean.rel.force <- list() 
Qminimumx <- list() 
QnbBC <- list() 
QnbBCBPbefore <- list() 
QnbBC.data <- list() 
QnBCBP <- list() 
QnbDistBCBPafter <- list() 
QnbDistBCBP.data <- list() 
Qnbins <- list() 
QpBCBCBP <- list() 
QpBCBP <- list() 
Qp_BC_ttest <- list() 
Qpolygon.data <- list() 
Qpos1.cross.all <- list() 
Qpos2.cross.all <- list() 
QPos.new <- list() 
Qreac <- list() 
QRI_BCBPafter <- list() 
QRI_BCBP.data <- list() 
QRI_BC.data <- list() 
Qsd.dendritic.length <- list() 
Qsd.drpBCBCBPsjeafter.ds <- list() 
Qsd.drpBCBPsjeafter.ds <- list() 
Qsd.drpBCsje.ds <- list() 
Qsdfield.dendritic.length.post <- list() 
Qsd.migration.dist <- list() 
Qsd.rel.force <- list() 
Qseed <- list() 
Qsynapses <- list() 
QTBCBCBP <- list() 
Qvert.d <- list() 
Qvert.d.sd <- list() 
Qx <- list() 


#for (m in 1:length(BCexc)) {


#################################
### Global parameters
#################################
a <- 0.10				# Fraction of cone which are blue cones (Roorda2001)
r <- 4					# BCBP radius (in mid peripheral retina) (Koyama 1992)
GdexcCone <- c(12,12,12,12,12)		# Radius of cone exclusion zone (mu metre)
GdexcBC <- c(30,38.3,28,33,36)		# Mean radius of blue cone exclusion zone (mu metre)
GdexcBC.sd <- c(6,3.9,3,7,7)		# Standard deviation of blue cone exclusion zone (mu metre)
GdexcBC.trunc <- c(17,32,23,18,23)	# Blue clone exclusion zone truncation, based on min observed nearest neighbour distance
GdexcBCBP <- c(18,22,19,18,22)		# Radius of BCBP exclusion zone (mu metre) (this parameter is important for the even distribution of links), used to be 8
GdexcBCBP.sd <- c(5,6,5,5,8)		# BCBP exclusion zone standard deviation
GdexcBCBP.trunc <- c(8.2,7.2,5.4,5.7,5.4)		# BCBP exclusion zone truncation, based on min. observed nearest neighbour distance
GINLdepth <- c(40,20,40,30,15)		# Depth of inner nuclear layer (Kouyama)
maxdendr <- 44				# Maxmimal horisontal dendrite length
Gvert.d.mean <- 10 + 0.33*GINLdepth	# Mean vertical distance from BC pedicles to BCPB perikaryons (Koyama 1992)
Gvert.d.sd <- GINLdepth*0.11		# Standard deviation of vertical distrance from BC pedicles to BCBP perikaryons
vert.d.thr <- 3*r+10			# Threshold for vertical distance between BCBP perikaryons and BC pedicles, if distance is greater there is no BCBP motility
Gmean.rel.force <- c(0.75,0.90,0.82,0.84,0.97)	# Mean value of relative force (friction/string)
sd.rel.force <- 0.1			# Standard deviation of relative force (friction/string)
# loopsize <- 49				# Number of iterations in simulation study
bin <- 10				# Bin size
nbins <- 10				# Number of bins
Gseed <- c(1,2,2,2,2)			# Seends for random number generation

#################################
### Looping over datafield, repeating the entire analysis for each field
#################################
x <- commandlineargs[1]
cat("Here starts the depth optimization for field",x,"\n",file=paste("LogField",x,".txt",sep=""),append=FALSE)

#for (x in 2:2) {		# This loop is very big

#seed <- c(1,2,2,2,2)
#set.seed(seed[x])


#################################
### Reading data
#################################
exptfile <- sprintf("%s/data/bc_bcbp_f%d.dat", AVLROOT, x)
##exptfile <- paste("data/bc_bcbp_f",x,".dat", sep="")
data <- read.table(exptfile)
#data <- data[sample(1:length(data[,1]),80),] ## This is temporary, only here to make the script run faster
data1 <- data[which(data[,3] == 1),c(1,2)]
data2 <- data[which(data[,3] == 2),c(1,2)]

#################################
### Setting field parameters
#################################
nBC <- sum(data[,3]==1)			# Number of blue cones
nBCBP <- sum(data[,3]==2)		# Number of blue cone bipolar cells
nC <- nBC*(1/a)				# Number of cones
maximumx <-  max(data[,1])		# Max x value (mu metre)
maximumy <-  max(data[,2])		# Max y value (mu metre)
minimumx <-  min(data[,1])		# Min x value (mu metre)
minimumy <-  min(data[,2])		# Min y value (mu metre)
lx <- maximumx - minimumx		# Horizontal dimension of retinal area (mu metre)
ly <- maximumy - minimumy		# Vertical dimension of retinal area (mu metre)

dexcCone <- GdexcCone[x]
dexcBC <- GdexcBC[x]
dexcBC.sd <- GdexcBC.sd[x]
dexcBC.trunc <- GdexcBC.trunc[x]
dexcBCBP <- GdexcBCBP[x]
dexcBCBP.sd <- GdexcBCBP.sd[x]
dexcBCBP.trunc <- GdexcBCBP.trunc[x]
INLdepth <- GINLdepth[x]
vert.d.mean <- Gvert.d.mean[x]
vert.d.sd <- Gvert.d.sd[x]
mean.rel.force <- Gmean.rel.force[x]
polygon <- matrix( c(0,0,0,ly,lx,ly,lx,0,0,0), nrow=5, ncol=2, byrow=T)
polygon.data <- matrix( c(minimumx,minimumy,minimumx,maximumy,maximumx,maximumy,maximumx,minimumy,minimumx,minimumy), nrow=5, ncol=2, byrow=T)

#################################
### Initializing variables to store output from simulation
#################################

drpBC <- matrix(NA, nrow=loopsize, ncol=nbins)				# Matrix with drp, for BCs, one entry for each simulation (nonedge)
drpBCall <-matrix(NA, nrow=loopsize, ncol=nbins)			# Matrix with drp, for BCs, one entry for each simulation (all)
drpBCsje <- list()							# List with drp, for BCs, one entry for each simulation (all, using sjedrp package)

drpBCBP <- matrix(NA, nrow=loopsize, ncol=nbins)			# Matrix with drp, for BCBPs, one entry for each simulation (nonedge)
drpBCBPall <- matrix(NA, nrow=loopsize, ncol=nbins)			# Matrix with drp, for BCBPs, one entry for each simulation (all)
drpBCBPsjeafter <- list() 						# List with drp, for BCBPs, one entry for each simulation (all, using sjedrp package)
drpBCBPsjebefore <- list()						# List with drp, for BCBPs, one entry for each simulation (all, using sjedrp package)

drpBCBCBP <- matrix(NA, nrow=loopsize, ncol=nbins)			# Matrix with drp, for BC to BCBPs, one entry for each simulation (nonedge BCs)
drpBCBCBPall <- matrix(NA, nrow=loopsize, ncol=nbins)			# Matrix with drp, for BC to BCBPs, one entry for each simulation (all)
drpBCBCBPold <- matrix(NA, nrow=loopsize, ncol=nbins)			# Matrix with drp, for BC to BCBPs, one entry for each simulation (old)
drpBCBCBPsjeafter <- list()						# List with drp, for BC to BCBPs, one entry for each simulation (all, using sjedrp package)
drpBCBCBPsjebefore <- list()						# List with drp, for BC to BCBPs, one entry for each simulation (all, using sjedrp package)

LBC <- matrix(NA, nrow=loopsize, ncol=(nbins*bin))			# Matrix with L function for BCs, one row for each simulation
LBCBPafter <- matrix(NA, nrow=loopsize, ncol=(nbins*bin))		# Matrix with L function for BCBPs, one row for each simulation
LBCBPbefore <- matrix(NA, nrow=loopsize, ncol=(nbins*bin))		# Matrix with L function for BCBPs, one row for each simulation
LBCBCBPafter <- matrix(NA, nrow=loopsize, ncol=(nbins*bin))		# Matrix with L function for cross between BCs and BCBPs, one row for each simulation
LBCBCBPbefore <- matrix(NA, nrow=loopsize, ncol=(nbins*bin))		# Matrix with L function for cross between BCs and BCBPs, one row for each simulation

GBC <- matrix(NA, nrow=loopsize, ncol=(nbins*bin))			# Matrix with G function for BCs, one row for each simulation
GBCBPbefore <- matrix(NA, nrow=loopsize, ncol=(nbins*bin))		# Matrix with G function for BCBPs, one row for each simulation
GBCBPafter <- matrix(NA, nrow=loopsize, ncol=(nbins*bin))		# Matrix with G function for BCBPs, one row for each simulation

synapses <- vector()							# Number of synapses in each simulation
movedBCBPs <- vector()							# Number of BCBPS which have moved in each simulation
mean.migration.dist <- vector()						# Average migration distance among BCBPs that move
sd.migration.dist <- vector()						# Standad deviation of migration distance among BCBPs that move

dendritic.length <- matrix(NA, ncol=nBCBP, nrow=loopsize)
dendritic.length.post <- matrix(NA, ncol=nBCBP, nrow=loopsize)
mean.dendritic.length <- vector()
mean.dendritic.length.post <- vector()
sd.dendritic.length <- vector()
sd.dendritic.length.post <- vector()

RI_BC <- vector()							# Regularity index for blue cones
RI_BCBPbefore <- vector()						# Regularity index for BCBP before migration
RI_BCBPafter <- vector()						# Regularity index for BCBP after migration


#################################
### Looping over parameter values
#################################
for (m in 1:length(depths)) {
#	for (o in 1:length(BCBPexc.sd)) {

#################################
### Setting seed
#################################
set.seed(Gseed[x])

#################################
### Generating cone mosaic
#################################

conePos <- position(nC,lx,ly,dexcCone,0,dexcCone)	# lx and ly are the field dimensions, nC is the number of cones and dexcCone is the cone exclusion zone (which is also the truncation in this case, where there is no standard deviation, the exclusion zone is constant


# mean.rel.force  <- forceratio[m]              # Mean radius of blue cone exclusion zone (mu metre)

vert.d.thr <- depths[m] 	# Vertical threshold depth, BCBP further away from PC pedicles than this will not move

#################################
### Simulations starts here, replications
#################################
for (k in 1:loopsize) {

##############################
### Differentiation of a sub set of cones into blue cones
##############################
BCindex <- bcDiff(nC,nBC,dexcBC,dexcBC.sd,dexcBC.trunc,conePos) # Index of cones which have now differentiated into blue cones
bcPos <- conePos[BCindex,]

################
### BC output ##
################

# DRP
drpBCsje[[k]] <- autodrp(bcPos[,1],bcPos[,2],nbins,bin)

# Calculating L function for BCs (this produces simulation output: LBC)
LBC[k,] <- lFunc(bcPos,polygon,nbins,bin)

# Calculating the distance to the nearest neighbour for each BC
nbDistBC <- nbDist(bcPos)

# Calculating G function for BCs (this produces simulation output: GBC)
GBC[k,] <- gFunc(nbDistBC,bin,nbins,nBC)

# Calculating regularity index
RI_BC[k] <- regInd(nbDistBC)


###################################
### Modeling the BCBPs
###################################
BCBPPos <- position(nBCBP,lx,ly,dexcBCBP,dexcBCBP.sd,dexcBCBP.trunc)

#####################################
## Dendrites between BCs and BCBPs
#####################################
# For each BCBP, the BC it connects to
reac_store <- connections(bcPos,BCBPPos,nBC,nBCBP,maxdendr)

# Matrix giving links between each BCBP to BCs
links <- cbind(BCBPPos,bcPos[reac_store,])

# even with perfect BCBP mosaic, too many cones have 3 or more connections (ass. nearest neighbour). This indicates that in real life, something more than nearest neighbour is at work.

####################################
## Modeling BCBPs movement
####################################

# Calculating length of dendrites
dendr.length <- apply(links,1,function(x) { sqrt((x[1]-x[3])^2+(x[2]-x[4])^2)})
dendr.length[which(is.na(dendr.length))] <- 0

# Post migration BC and BCBP positions
Pos.new <- migration(nBCBP,vert.d.mean,vert.d.sd,vert.d.thr,mean.rel.force,sd.rel.force,dexcBCBP.trunc,dendr.length,links,BCBPPos)

# Post migration BCBP positions
BCBPPos.new <- Pos.new[(nBC+1):(nBC+nBCBP),]

# New links between BC and BCBPs
links.new <- cbind(BCBPPos.new,bcPos[reac_store,])

# Calculating dendritic length after migration
dendr.length.post <- apply(links.new,1,function(x) { sqrt((x[1]-x[3])^2+(x[2]-x[4])^2)})
dendr.length.post[which(is.na(dendr.length.post))] <- 0

##########################
### BCBP output
##########################

### Before migration

# DRP
drpBCBPsjebefore[[k]] <- autodrp(BCBPPos[,1],BCBPPos[,2],nbins,bin)

# L function
LBCBPbefore[k,] <- lFunc(BCBPPos,polygon,nbins,bin)

# Finding nearest neighbour distances
nbDistBCBPbefore <- nbDistBCBP(BCBPPos)

# Calculate G function for BCBP, before
GBCBPbefore[k,] <- gFunc(nbDistBCBPbefore,bin,nbins,nBCBP)

# Calculating regularity index, before
RI_BCBPbefore[k] <- regInd(nbDistBCBPbefore)

### After migration

# Saving dendritic length
dendritic.length[k,] <- abs(dendr.length)
dendritic.length.post[k,] <- abs(dendr.length.post)

# Mean dendritic length before and after motility 
mean.dendritic.length[k] <- mean(abs(dendr.length), na.rm=TRUE)
mean.dendritic.length.post[k] <- mean(abs(dendr.length.post), na.rm=TRUE)

# Standard deviation of dendritic length, before and after
sd.dendritic.length[k] <- sd(abs(dendr.length), na.rm=TRUE)
sd.dendritic.length.post[k] <- sd(abs(dendr.length.post), na.rm=TRUE)

# Average migration distance
dendr.length.diff <- dendr.length-dendr.length.post
mean.migration.dist[k] <- mean(dendr.length.diff[which(dendr.length.diff != 0)])	# Average migration distance amon BCBPs that move
sd.migration.dist[k] <- sd(dendr.length.diff[which(dendr.length.diff != 0)])	# Average migration distance amon BCBPs that move

# Number of BCBPs which have moved
movedBCBPs[k] <- sum(dendr.length-dendr.length.post != 0)

# DRP
drpBCBPsjeafter[[k]] <- autodrp(BCBPPos.new[,1],BCBPPos.new[,2],nbins,bin)

# Calculating L function for BCBPs
LBCBPafter[k,] <- lFunc(BCBPPos.new,polygon,nbins,bin)

# Nearest neighbour distance for BCBPs after migration
nbDistBCBPafter <- nbDist(BCBPPos.new)

# Calculating G function for BCBPs
GBCBPafter[k,] <- gFunc(nbDistBCBPafter,bin,nbins,nBCBP)

# Calculating regularity index
RI_BCBPafter[k] <- regInd(nbDistBCBPafter)

######################
### Cross output
######################

# Counting the number of synaptic connections
synapses[k] <- nBCBP-sum(is.na(reac_store))

# DRP before
drpBCBCBPsjebefore[[k]] <- crossdrp(BCBPPos[,1],BCBPPos[,2],bcPos[,1],bcPos[,2],nbins,bin)

# DRP after
drpBCBCBPsjeafter[[k]] <- crossdrp(BCBPPos.new[,1],BCBPPos.new[,2],bcPos[,1],bcPos[,2],nbins,bin)

# Calculating the L function for BCBCBPs, before
LBCBCBPbefore[k,] <- lFuncCross(bcPos,BCBPPos,polygon,nbins,bin)

# Calculating the L function for BCBCBPs, after
LBCBCBPafter[k,] <- lFuncCross(bcPos,BCBPPos.new,polygon,nbins,bin)

print(paste("End of loop k =",k))

cat("Here ends loop",k,"\n",file=paste("LogField",x,".txt",sep=""),append=TRUE)


}		# Here ends the iteration loop


#####################
#####################
### Results from all iterations
#####################
#####################

#####################
### Dendritic lengths
#####################
meanfield.dendritic.length <- mean(mean.dendritic.length)		# The average dendritic length of a field and averaged over simulations
meanfield.dendritic.length.post <- mean(mean.dendritic.length.post)

sdfield.dendritic.length <- sd(mean.dendritic.length)
sdfield.dendritic.length.post <- sd(mean.dendritic.length.post)

sd.dendritic.length <- mean(sd.dendritic.length)
sd.dendritic.length.post <- mean(sd.dendritic.length.post)

#####################
## How many BCBPs that have moved
#####################
mean.prop.movedBCBPs <- mean(movedBCBPs/nBCBP)	# Mean proportion of migrated BCBPs
sd.prop.movedBCBPs <- sd(movedBCBPs/nBCBP)	# SD of proportion of migrated BCBPs

#####################
## How many synaptic connections
#####################
mean.synapses <- mean(synapses)	# Mean number of synapses (across iterations)
sd.synapses <- sd(synapses)	# SD of number of synapses (across iterations)

#####################
## Migration distance
#####################
meanfield.migration.dist <- mean(mean.migration.dist)
sdfield.migration.dist <- sd(mean.migration.dist)
sd.migration.dist <- mean(sd.migration.dist)

print("check 1")

#####################
## Simulated data L and G function, envelopes
#####################

# Calculating the envelope of the L function for the simulated BCs
LBCenvelope <- apply(LBC, 2, function(x) {quantile(x, probs=c(0.025,0.975))})

# Calculating the envelope of the L function for the simulated BCBPs, before
LBCBPbeforeenvelope <- apply(LBCBPbefore, 2, function(x) {quantile(x, probs=c(0.025,0.975))})

# Calculating the envelope of the L function for the simulated BCBPs, after
LBCBPafterenvelope <- apply(LBCBPafter, 2, function(x) {quantile(x, probs=c(0.025,0.975))})

F# Calculating the envelope of the L function for the simulated BC to BCBPs, before
LBCBCBPbeforeenvelope <- apply(LBCBCBPbefore, 2, function(x) {quantile(x, probs=c(0.025,0.975))})

# Calculating the envelope of the L function for the simulated BC to BCBPs, after
LBCBCBPafterenvelope <- apply(LBCBCBPafter, 2, function(x) {quantile(x, probs=c(0.025,0.975))})

# Calculating the envelope of the G function for the simulated BCs
GBCenvelope <- apply(GBC, 2, function(x) {quantile(x, probs=c(0.025,0.975))})

# Calculating the envelope of the G function for the simulated BCBPs, before
GBCBPbeforeenvelope <- apply(GBCBPbefore, 2, function(x) {quantile(x, probs=c(0.025,0.975))})

# Calculating the envelope of the G function for the simulated BCBPs, after
GBCBPafterenvelope <- apply(GBCBPafter, 2, function(x) {quantile(x, probs=c(0.025,0.975))})


#####################
## Simulated data drp
#####################

## BCs

drpBCsje.ds <- matrix(NA, nrow=loopsize, ncol=nbins)	#rep(0,nbins)
drpBCsje.effrad <- vector()
drpBCsje.p <- vector()
drpBCsje.maxr <- vector()
drpBCsje.k <- vector()
drpBCsje.density <- vector()

for (i in 1:loopsize) {
	drpBCsje.ds[i,] <- drpBCsje[[i]]$ds
	drpBCsje.effrad <- c(drpBCsje.effrad,drpBCsje[[i]]$effrad)
	drpBCsje.p <- c(drpBCsje.p,drpBCsje[[i]]$p)
	drpBCsje.maxr <- c(drpBCsje.maxr,drpBCsje[[i]]$maxr)
	drpBCsje.k <- c(drpBCsje.k,drpBCsje[[i]]$k)
	drpBCsje.density <- c(drpBCsje.density,drpBCsje[[i]]$density)
}

mean.drpBCsje.ds <- apply(drpBCsje.ds,2,mean)
sd.drpBCsje.ds <- apply(drpBCsje.ds,2,sd)

mean.drpBCsje.effrad <- mean(drpBCsje.effrad)
mean.drpBCsje.p <- mean(drpBCsje.p)
mean.drpBCsje.maxr <- mean(drpBCsje.maxr)
mean.drpBCsje.k <- mean(drpBCsje.k) 
mean.drpBCsje.density <- mean(drpBCsje.density)

print("check 2")
## BCBPs (after)

drpBCBPsjeafter.ds <- matrix(NA, nrow=loopsize, ncol=nbins)	#rep(0,nbins)
drpBCBPsjeafter.effrad <- vector()
drpBCBPsjeafter.p <- vector()
drpBCBPsjeafter.maxr <- vector()
drpBCBPsjeafter.k <- vector()
drpBCBPsjeafter.density <- vector()

for (i in 1:loopsize) {
	drpBCBPsjeafter.ds[i,] <- drpBCBPsjeafter[[i]]$ds
	drpBCBPsjeafter.effrad <- c(drpBCBPsjeafter.effrad,drpBCBPsjeafter[[i]]$effrad)
	drpBCBPsjeafter.p <- c(drpBCBPsjeafter.p,drpBCBPsjeafter[[i]]$p)
	drpBCBPsjeafter.maxr <- c(drpBCBPsjeafter.maxr,drpBCBPsjeafter[[i]]$maxr)
	drpBCBPsjeafter.k <- c(drpBCBPsjeafter.k,drpBCBPsjeafter[[i]]$k)
	drpBCBPsjeafter.density <- c(drpBCBPsjeafter.density,drpBCBPsjeafter[[i]]$density)
}

mean.drpBCBPsjeafter.ds <- apply(drpBCBPsjeafter.ds,2,mean)
sd.drpBCBPsjeafter.ds <- apply(drpBCBPsjeafter.ds,2,sd)
mean.drpBCBPsjeafter.effrad <- mean(drpBCBPsjeafter.effrad)
mean.drpBCBPsjeafter.p <- mean(drpBCBPsjeafter.p)
mean.drpBCBPsjeafter.maxr <- mean(drpBCBPsjeafter.maxr)
mean.drpBCBPsjeafter.k <- mean(drpBCBPsjeafter.k) 
mean.drpBCBPsjeafter.density <- mean(drpBCBPsjeafter.density)


## BCBPs (before)

drpBCBPsjebefore.ds <- matrix(NA, nrow=loopsize, ncol=nbins)
drpBCBPsjebefore.effrad <- vector()
drpBCBPsjebefore.p <- vector()
drpBCBPsjebefore.maxr <- vector()
drpBCBPsjebefore.k <- vector()
drpBCBPsjebefore.density <- vector()

for (i in 1:loopsize) {
	drpBCBPsjebefore.ds[i,] <- drpBCBPsjebefore[[i]]$ds
	drpBCBPsjebefore.effrad <- c(drpBCBPsjebefore.effrad,drpBCBPsjebefore[[i]]$effrad)
	drpBCBPsjebefore.p <- c(drpBCBPsjebefore.p,drpBCBPsjebefore[[i]]$p)
	drpBCBPsjebefore.maxr <- c(drpBCBPsjebefore.maxr,drpBCBPsjebefore[[i]]$maxr)
	drpBCBPsjebefore.k <- c(drpBCBPsjebefore.k,drpBCBPsjebefore[[i]]$k)
	drpBCBPsjebefore.density <- c(drpBCBPsjebefore.density,drpBCBPsjebefore[[i]]$density)
}

mean.drpBCBPsjebefore.ds <- apply(drpBCBPsjebefore.ds,2,mean)
sd.drpBCBPsjebefore.ds <- apply(drpBCBPsjebefore.ds,2,sd)
mean.drpBCBPsjebefore.effrad <- mean(drpBCBPsjebefore.effrad)
mean.drpBCBPsjebefore.p <- mean(drpBCBPsjebefore.p)
mean.drpBCBPsjebefore.maxr <- mean(drpBCBPsjebefore.maxr)
mean.drpBCBPsjebefore.k <- mean(drpBCBPsjebefore.k) 
mean.drpBCBPsjebefore.density <- mean(drpBCBPsjebefore.density)


## BCBCBPs (after)
print("check 3")

drpBCBCBPsjeafter.ds <- matrix(NA, nrow=loopsize, ncol=nbins)
drpBCBCBPsjeafter.effrad <- vector()
drpBCBCBPsjeafter.p <- vector()
drpBCBCBPsjeafter.maxr <- vector()
drpBCBCBPsjeafter.k <- vector()
drpBCBCBPsjeafter.density <- vector()

for (i in 1:loopsize) {
	drpBCBCBPsjeafter.ds[i,] <- drpBCBCBPsjeafter[[i]]$ds
	drpBCBCBPsjeafter.effrad <- c(drpBCBCBPsjeafter.effrad,drpBCBCBPsjeafter[[i]]$effrad)
	drpBCBCBPsjeafter.p <- c(drpBCBCBPsjeafter.p,drpBCBCBPsjeafter[[i]]$p)
	drpBCBCBPsjeafter.maxr <- c(drpBCBCBPsjeafter.maxr,drpBCBCBPsjeafter[[i]]$maxr)
	drpBCBCBPsjeafter.k <- c(drpBCBCBPsjeafter.k,drpBCBCBPsjeafter[[i]]$k)
	drpBCBCBPsjeafter.density <- c(drpBCBCBPsjeafter.density,drpBCBCBPsjeafter[[i]]$density)
}


mean.drpBCBCBPsjeafter.ds <- apply(drpBCBCBPsjeafter.ds,2,mean)
sd.drpBCBCBPsjeafter.ds <- apply(drpBCBCBPsjeafter.ds,2,sd)
mean.drpBCBCBPsjeafter.effrad <- mean(drpBCBCBPsjeafter.effrad)
mean.drpBCBCBPsjeafter.p <- mean(drpBCBCBPsjeafter.p)
mean.drpBCBCBPsjeafter.maxr <- mean(drpBCBCBPsjeafter.maxr)
mean.drpBCBCBPsjeafter.k <- mean(drpBCBCBPsjeafter.k) 
mean.drpBCBCBPsjeafter.density <- mean(drpBCBCBPsjeafter.density)


## BCBCBPs (before)

drpBCBCBPsjebefore.ds <- matrix(NA, nrow=loopsize, ncol=nbins)
drpBCBCBPsjebefore.effrad <- vector()
drpBCBCBPsjebefore.p <- vector()
drpBCBCBPsjebefore.maxr <- vector()
drpBCBCBPsjebefore.k <- vector()
drpBCBCBPsjebefore.density <- vector()

for (i in 1:loopsize) {
	drpBCBCBPsjebefore.ds[i,] <- drpBCBCBPsjebefore[[i]]$ds
	drpBCBCBPsjebefore.effrad <- c(drpBCBCBPsjebefore.effrad,drpBCBCBPsjebefore[[i]]$effrad)
	drpBCBCBPsjebefore.p <- c(drpBCBCBPsjebefore.p,drpBCBCBPsjebefore[[i]]$p)
	drpBCBCBPsjebefore.maxr <- c(drpBCBCBPsjebefore.maxr,drpBCBCBPsjebefore[[i]]$maxr)
	drpBCBCBPsjebefore.k <- c(drpBCBCBPsjebefore.k,drpBCBCBPsjebefore[[i]]$k)
	drpBCBCBPsjebefore.density <- c(drpBCBCBPsjebefore.density,drpBCBCBPsjebefore[[i]]$density)
}

mean.drpBCBCBPsjebefore.ds <- apply(drpBCBCBPsjebefore.ds,2,mean)
mean.drpBCBCBPsjebefore.effrad <- mean(drpBCBCBPsjebefore.effrad)
sd.drpBCBCBPsjebefore.ds <- apply(drpBCBCBPsjebefore.ds,2,sd)
mean.drpBCBCBPsjebefore.p <- mean(drpBCBCBPsjebefore.p)
mean.drpBCBCBPsjebefore.maxr <- mean(drpBCBCBPsjebefore.maxr)
mean.drpBCBCBPsjebefore.k <- mean(drpBCBCBPsjebefore.k) 
mean.drpBCBCBPsjebefore.density <- mean(drpBCBCBPsjebefore.density)


##########################################
##### real data drps
##########################################
print("check 4")

### BCs
drpBCdatasje <- autodrp(data1[,1],data1[,2],nbins,bin)

### BCBPs
drpBCBPdatasje <- autodrp(data2[,1],data2[,2],nbins,bin)

## BC to BCBPs
drpBCBCBPdatasje <- crossdrp(data1[,1],data1[,2],data2[,1],data2[,2],nbins,bin)


##########################################
### real data L and G function and RI
##########################################

# Calculating the L function for the real BCs
LBC.data <- lFunc(data1,polygon.data,nbins,bin)

# Calculating the L function for the real BCBPs
LBCBP.data <- lFunc(data2,polygon.data,nbins,bin)

# Calculating the L function for the real BCBCBPs
LBCBCBP.data <- lFuncCross(data1,data2,polygon.data,nbins,bin)

# Nearest neighbour distances for data BCs
nbDistBC.data <- nbDist(data1)

# G function for BC
GBCdata <- gFunc(nbDistBC.data,bin,nbins,nBC)

# Calculating regularity index for BC data
RI_BC.data <- regInd(nbDistBC.data)

# Nearest neighbour distances for data BCs
nbDistBCBP.data <- nbDist(data2)

# G function for BCBP
GBCBPdata <- gFunc(nbDistBCBP.data,bin,nbins,nBCBP)

# Calculating regularity index for BCBP data
RI_BCBP.data <- regInd(nbDistBCBP.data)


###########################
## Testing hypotheses
###########################


### K function hypothesis: Is model different from data

# ulimit <- 50	# Upper limit for region of interest (mu meter)

## BC

KBCall <- rbind(pi*(LBC.data)^2,pi*(LBC)^2)

KBCmean <- matrix(NA, nrow=(loopsize+1), ncol=ulimit)
for (i in 1:(loopsize+1)) {
	for (j in 1:ulimit) {
		KBCmean[i,j] <- mean(KBCall[-i,j])
	}
}

TBC <- vector()
for (i in 1:(loopsize+1)) {
	TBC[i] <- sum((sqrt(KBCall[i,1:ulimit]) - sqrt(KBCmean[i,]))^2)
}
pBC[m,1] <- ((loopsize+1-rank(TBC))/(loopsize+1))[1]


## BCBP

KBCBPall <- rbind(pi*(LBCBP.data)^2,pi*(LBCBPafter)^2)

KBCBPmean <- matrix(NA, nrow=(loopsize+1), ncol=ulimit)
for (i in 1:(loopsize+1)) {
	for (j in 1:ulimit) {
		KBCBPmean[i,j] <- mean(KBCBPall[-i,j])
	}
}

TBCBP <- vector()
for (i in 1:(loopsize+1)) {
	TBCBP[i] <- sum((sqrt(KBCBPall[i,1:ulimit]) - sqrt(KBCBPmean[i,]))^2)
}
pBCBP[m,1] <- ((loopsize+1-rank(TBCBP))/(loopsize+1))[1]


## BCBCBP


KBCBCBPall <- rbind(pi*(LBCBCBP.data)^2,pi*(LBCBCBPafter)^2)

KBCBCBPmean <- matrix(NA, nrow=(loopsize+1), ncol=ulimit)
for (i in 1:(loopsize+1)) {
	for (j in 1:ulimit) {
		KBCBCBPmean[i,j] <- mean(KBCBCBPall[-i,j])
	}
}

print("check 6")

TBCBCBP <- vector()
for (i in 1:(loopsize+1)) {
	TBCBCBP[i] <- sum(  (   sqrt(KBCBCBPall[i,1:ulimit]) - sqrt(KBCBCBPmean[i,])  )^2   )
}
pBCBCBP[m,1] <- ((loopsize+1-rank(TBCBCBP))/(loopsize+1))[1]

pBCBCBP

### Testing whether the RI is different from data
RI_BC_ttest <- t.test(RI_BC,mu=RI_BC.data)$p.value		# two-sided t test
RI_BCBP_ttest <- t.test(RI_BCBPafter,mu=RI_BCBP.data)$p.value	# two-sided t test

### Testing whether the eff rad is different from data
effrad_BC_ttest <- t.test(drpBCsje.effrad, mu=drpBCdatasje$effrad)$p.value
effrad_BCBP_ttest <- t.test(drpBCBPsjeafter.effrad, mu=drpBCBPdatasje$effrad)$p.value

### Testing whether the packing is different from data
p_BC_ttest <- t.test(drpBCsje.p, mu=drpBCdatasje$p)$p.value
p_BCBP_ttest <- t.test(drpBCBPsjeafter.p, mu=drpBCBPdatasje$p)$p.value
p_BCBCBP_ttest <- t.test(drpBCBCBPsjeafter.p, mu=drpBCBCBPdatasje$p)$p.value

### Testing whether the maxr is different from data
maxr_BC_ttest <- t.test(drpBCsje.maxr, mu=drpBCdatasje$maxr)$p.value
maxr_BCBP_ttest <- t.test(drpBCBPsjeafter.maxr, mu=drpBCBPdatasje$maxr)$p.value
maxr_BCBCBP_ttest <- t.test(drpBCBCBPsjeafter.maxr, mu=drpBCBCBPdatasje$maxr)$p.value

#Ending the loop over parameter values
cat("Here endeth: depth",depths[m],"\n",file=paste("LogField",x,".txt",sep=""),append=TRUE)

## Making plots to get a visual impression of the output


postscript(file=paste("field",x,"_", vert.d.thr,"_",loopsize,"_",ulimit,".ps",sep=""))	#   postscript(file=paste("field1,",y,",",z,".pdf",sep="")) #	 	
par(mfrow=c(3,3))

plot(LBC.data, type="l", lwd=1, col="red", xlab="mu meter", ylab="L", main="L function for BC")
lines(LBCenvelope[1,])
lines(LBCenvelope[2,])
text(80,20,label=paste("p =",round(pBC[m,1],digits=3)))

plot(LBCBP.data, type="l", lwd=1, col="red", xlab="mu meter", ylab="L", main="L function for BCBP")
lines(LBCBPafterenvelope[1,])
lines(LBCBPafterenvelope[2,])
text(80,20,label=paste("p =",round(pBCBP[m,1],digits=3)))

plot(LBCBCBP.data, type="l", lwd=1, col="red", xlab="mu meter", ylab="L", main="L function for BC to BCBP")
lines(LBCBCBPafterenvelope[1,])
lines(LBCBCBPafterenvelope[2,])
text(80,20,label=paste("p =",round(pBCBCBP[m,1],digits=3)))

barplot(mean.drpBCsje.ds, main="BC", sub=paste("Dens",round(mean.drpBCsje.density, digits=5),"Eff.rad",round(mean.drpBCsje.effrad,1),"Pack",round(mean.drpBCsje.p,2),"Maxr",round(mean.drpBCsje.maxr,1),"RI",round(mean.drpBCsje.k,2)), cex.main=0.8, cex.sub=0.8, xlab=)

plot(drpBCdatasje)

barplot(mean.drpBCBPsjeafter.ds, main="BCBP", sub=paste("Dens",round(mean.drpBCBPsjeafter.density, digits=5),"Eff.rad",round(mean.drpBCBPsjeafter.effrad,1),"Pack",round(mean.drpBCBPsjeafter.p,2),"Maxr",round(mean.drpBCBPsjeafter.maxr,1),"RI",round(mean.drpBCBPsjeafter.k,2)), cex.main=0.8, cex.sub=0.8)

plot(drpBCBPdatasje)

barplot(mean.drpBCBCBPsjeafter.ds, main="BCBCBP", sub=paste("Dens",round(mean.drpBCBCBPsjeafter.density, digits=5),"Eff.rad",round(mean.drpBCBCBPsjeafter.effrad,1),"Pack",round(mean.drpBCBCBPsjeafter.p,2),"Maxr",round(mean.drpBCBCBPsjeafter.maxr,1),"RI",round(mean.drpBCBCBPsjeafter.k,2)), cex.main=0.8, cex.sub=0.8)

plot(drpBCBCBPdatasje)

dev.off()


	}
#}

cat("Here ends the depth optimization for field",x,"\n",file=paste("LogField",x,".txt",sep=""),append=TRUE)

# Collecting all p values in one matrix
pvalues <- cbind(pBC,pBCBP,pBCBCBP)
colnames(pvalues) <- c("pBC","pBCBP","pBCBCBP")
rownames(pvalues) <- depths

##############
### Redirecting p values for BC L function to a file
#############
sink(paste("pvalues_depth_",commandlineargs[1],"_",commandlineargs[2],"_",commandlineargs[3],"_",commandlineargs[4],"_",commandlineargs[5],"_",commandlineargs[6],"_",dexcBC,"_",dexcBC.sd,"_",dexcBC.trunc,"_",dexcBCBP,"_",dexcBCBP.sd,"_",dexcBCBP.trunc,".txt",sep=""))

print("p values")
pvalues
#BC
#print("pBCBP")
#pBCBP
#print("pBCBCBP")
#pBCBCBP
sink()

cat("p values printed to file",file=paste("LogField",x,".txt",sep=""),append=TRUE)
###############
## Merging fields
###############

print("check 7")

QBCBPPos.new[[x]] <- BCBPPos.new
Qa[[x]] <- a
Qbin[[x]] <- bin
Qdata1[[x]] <- data1
Qdendritic.length[[x]] <- dendritic.length
Qdendr.length[[x]] <- dendr.length
QdexcBC[[x]] <- dexcBC
QdexcBCBP[[x]] <- dexcBCBP
QdexcBC.sd[[x]] <- dexcBC.sd
QdrpBCall[[x]] <- drpBCall
QdrpBCBCBPall[[x]] <- drpBCBCBPall
QdrpBCBCBPold[[x]] <- drpBCBCBPold
QdrpBCBCBPsjeafter.density[[x]] <- drpBCBCBPsjeafter.density
QdrpBCBCBPsjeafter.effrad[[x]] <- drpBCBCBPsjeafter.effrad
QdrpBCBCBPsjeafter.maxr[[x]] <- drpBCBCBPsjeafter.maxr
QdrpBCBCBPsjebefore[[x]] <- drpBCBCBPsjebefore
QdrpBCBCBPsjebefore.ds[[x]] <- drpBCBCBPsjebefore.ds
QdrpBCBCBPsjebefore.k[[x]] <- drpBCBCBPsjebefore.k
QdrpBCBCBPsjebefore.p[[x]] <- drpBCBCBPsjebefore.p
QdrpBCBPall[[x]] <- drpBCBPall
QdrpBCBPsjeafter[[x]] <- drpBCBPsjeafter
QdrpBCBPsjeafter.ds[[x]] <- drpBCBPsjeafter.ds
QdrpBCBPsjeafter.k[[x]] <- drpBCBPsjeafter.k
QdrpBCBPsjeafter.p[[x]] <- drpBCBPsjeafter.p
QdrpBCBPsjebefore.density[[x]] <- drpBCBPsjebefore.density
QdrpBCBPsjebefore.effrad[[x]] <- drpBCBPsjebefore.effrad
QdrpBCBPsjebefore.maxr[[x]] <- drpBCBPsjebefore.maxr
QdrpBCdatasje[[x]] <- drpBCdatasje
QdrpBCsje.density[[x]] <- drpBCsje.density
QdrpBCsje.effrad[[x]] <- drpBCsje.effrad
QdrpBCsje.maxr[[x]] <- drpBCsje.maxr
Qeffrad_BCBP_ttest[[x]] <- effrad_BCBP_ttest
QGBCBPafter[[x]] <- GBCBPafter
QGBCBPbefore[[x]] <- GBCBPbefore
QGBCBPdata[[x]] <- GBCBPdata
QGBCenvelope[[x]] <- GBCenvelope
QINLdepth[[x]] <- INLdepth
Qk[[x]] <- k
QKBCall[[x]] <- KBCall
QKBCBCBPall[[x]] <- KBCBCBPall
QKBCBPall[[x]] <- KBCBPall
QKBCBPmean[[x]] <- KBCBPmean
QKBCmean[[x]] <- KBCmean
QLBCBCBPafter[[x]] <- LBCBCBPafter
QLBCBCBPbefore[[x]] <- LBCBCBPbefore
QLBCBCBP.data[[x]] <- LBCBCBP.data
QLBCBPafterenvelope[[x]] <- LBCBPafterenvelope
QLBCBPbeforeenvelope[[x]] <- LBCBPbeforeenvelope
QLBC.data[[x]] <- LBC.data
Qlinks[[x]] <- links
Qloopsize[[x]] <- loopsize
Qly[[x]] <- ly
Qmaxdendr[[x]] <- maxdendr
Qmaximumy[[x]] <- maximumy
Qmaxr_BCBP_ttest[[x]] <- maxr_BCBP_ttest
Qmean.dendritic.length[[x]] <- mean.dendritic.length
Qmean.drpBCBCBPsjeafter.density[[x]] <- mean.drpBCBCBPsjeafter.density
Qmean.drpBCBCBPsjeafter.effrad[[x]] <- mean.drpBCBCBPsjeafter.effrad
Qmean.drpBCBCBPsjeafter.maxr[[x]] <- mean.drpBCBCBPsjeafter.maxr
Qmean.drpBCBCBPsjebefore.density[[x]] <- mean.drpBCBCBPsjebefore.density
Qmean.drpBCBCBPsjebefore.effrad[[x]] <- mean.drpBCBCBPsjebefore.effrad
Qmean.drpBCBCBPsjebefore.maxr[[x]] <- mean.drpBCBCBPsjebefore.maxr
Qmean.drpBCBPsjeafter.density[[x]] <- mean.drpBCBPsjeafter.density
Qmean.drpBCBPsjeafter.effrad[[x]] <- mean.drpBCBPsjeafter.effrad
Qmean.drpBCBPsjeafter.maxr[[x]] <- mean.drpBCBPsjeafter.maxr
Qmean.drpBCBPsjebefore.density[[x]] <- mean.drpBCBPsjebefore.density
Qmean.drpBCBPsjebefore.effrad[[x]] <- mean.drpBCBPsjebefore.effrad
Qmean.drpBCBPsjebefore.maxr[[x]] <- mean.drpBCBPsjebefore.maxr
Qmean.drpBCsje.density[[x]] <- mean.drpBCsje.density
Qmean.drpBCsje.effrad[[x]] <- mean.drpBCsje.effrad
Qmean.drpBCsje.maxr[[x]] <- mean.drpBCsje.maxr
Qmeanfield.dendritic.length[[x]] <- meanfield.dendritic.length
Qmeanfield.migration.dist[[x]] <- meanfield.migration.dist
Qmean.prop.movedBCBPs[[x]] <- mean.prop.movedBCBPs
Qmean.synapses[[x]] <- mean.synapses
Qminimumy[[x]] <- minimumy
QmovedBCBPs[[x]] <- movedBCBPs
QnBC[[x]] <- nBC
QnbDistBC[[x]] <- nbDistBC
QnbDistBCBPbefore[[x]] <- nbDistBCBPbefore
QnbDistBC.data[[x]] <- nbDistBC.data
QnC[[x]] <- nC
QpBC[[x]] <- pBC
Qp_BCBCBP_ttest[[x]] <- p_BCBCBP_ttest
Qp_BCBP_ttest[[x]] <- p_BCBP_ttest
Qpolygon[[x]] <- polygon
Qr[[x]] <- r
Qreac_store[[x]] <- reac_store
QRI_BC[[x]] <- RI_BC
QRI_BCBPbefore[[x]] <- RI_BCBPbefore
QRI_BCBP_ttest[[x]] <- RI_BCBP_ttest
QRI_BC_ttest[[x]] <- RI_BC_ttest
Qsd.dendritic.length.post[[x]] <- sd.dendritic.length.post
Qsd.drpBCBCBPsjebefore.ds[[x]] <- sd.drpBCBCBPsjebefore.ds
Qsd.drpBCBPsjebefore.ds[[x]] <- sd.drpBCBPsjebefore.ds
Qsdfield.dendritic.length[[x]] <- sdfield.dendritic.length
Qsdfield.migration.dist[[x]] <- sdfield.migration.dist
Qsd.prop.movedBCBPs[[x]] <- sd.prop.movedBCBPs
Qsd.synapses[[x]] <- sd.synapses
QTBC[[x]] <- TBC
QTBCBP[[x]] <- TBCBP
Qulimit[[x]] <- ulimit
Qvert.d.thr[[x]] <- vert.d.thr
QBCBPPos[[x]] <- BCBPPos
QBCindex[[x]] <- BCindex
QconePos[[x]] <- conePos
Qdata[[x]] <- data
Qdata2[[x]] <- data2
Qdendritic.length.post[[x]] <- dendritic.length.post
Qdendr.length.diff[[x]] <- dendr.length.diff
Qdendr.length.post[[x]] <- dendr.length.post
QdexcBC.trunc[[x]] <- dexcBC.trunc
QdexcBCBP.trunc[[x]] <- dexcBCBP.trunc
QdexcBCBP.sd[[x]] <- dexcBCBP.sd
QdexcCone[[x]] <- dexcCone
QdrpBC[[x]] <- drpBC
QdrpBCBCBP[[x]] <- drpBCBCBP
QdrpBCBCBPdatasje[[x]] <- drpBCBCBPdatasje
QdrpBCBCBPsjeafter[[x]] <- drpBCBCBPsjeafter
QdrpBCBCBPsjeafter.ds[[x]] <- drpBCBCBPsjeafter.ds
QdrpBCBCBPsjeafter.k[[x]] <- drpBCBCBPsjeafter.k
QdrpBCBCBPsjeafter.p[[x]] <- drpBCBCBPsjeafter.p
QdrpBCBCBPsjebefore.density[[x]] <- drpBCBCBPsjebefore.density
QdrpBCBCBPsjebefore.effrad[[x]] <- drpBCBCBPsjebefore.effrad
QdrpBCBCBPsjebefore.maxr[[x]] <- drpBCBCBPsjebefore.maxr
QdrpBCBP[[x]] <- drpBCBP
QdrpBCBPdatasje[[x]] <- drpBCBPdatasje
QdrpBCBPsjeafter.density[[x]] <- drpBCBPsjeafter.density
QdrpBCBPsjeafter.effrad[[x]] <- drpBCBPsjeafter.effrad
QdrpBCBPsjeafter.maxr[[x]] <- drpBCBPsjeafter.maxr
QdrpBCBPsjebefore[[x]] <- drpBCBPsjebefore
QdrpBCBPsjebefore.ds[[x]] <- drpBCBPsjebefore.ds
QdrpBCBPsjebefore.k[[x]] <- drpBCBPsjebefore.k
QdrpBCBPsjebefore.p[[x]] <- drpBCBPsjebefore.p
QdrpBCsje[[x]] <- drpBCsje
QdrpBCsje.ds[[x]] <- drpBCsje.ds
QdrpBCsje.k[[x]] <- drpBCsje.k
QdrpBCsje.p[[x]] <- drpBCsje.p
Qeffrad_BC_ttest[[x]] <- effrad_BC_ttest
QGBC[[x]] <- GBC
QGBCBPafterenvelope[[x]] <- GBCBPafterenvelope
QGBCBPbeforeenvelope[[x]] <- GBCBPbeforeenvelope
QGBCdata[[x]] <- GBCdata
Qi[[x]] <- i
Qj[[x]] <- j
QKBCBCBPmean[[x]] <- KBCBCBPmean
QLBC[[x]] <- LBC
QLBCBCBPafterenvelope[[x]] <- LBCBCBPafterenvelope
QLBCBCBPbeforeenvelope[[x]] <- LBCBCBPbeforeenvelope
QLBCBPafter[[x]] <- LBCBPafter
QLBCBPbefore[[x]] <- LBCBPbefore
QLBCBP.data[[x]] <- LBCBP.data
QLBCenvelope[[x]] <- LBCenvelope
Qlinks.new[[x]] <- links.new
Qlx[[x]] <- lx
Qmaximumx[[x]] <- maximumx
Qmaxr_BCBCBP_ttest[[x]] <- maxr_BCBCBP_ttest
Qmaxr_BC_ttest[[x]] <- maxr_BC_ttest
Qmean.dendritic.length.post[[x]] <- mean.dendritic.length.post
Qmean.drpBCBCBPsjeafter.ds[[x]] <- mean.drpBCBCBPsjeafter.ds
Qmean.drpBCBCBPsjeafter.k[[x]] <- mean.drpBCBCBPsjeafter.k
Qmean.drpBCBCBPsjeafter.p[[x]] <- mean.drpBCBCBPsjeafter.p
Qmean.drpBCBCBPsjebefore.ds[[x]] <- mean.drpBCBCBPsjebefore.ds
Qmean.drpBCBCBPsjebefore.k[[x]] <- mean.drpBCBCBPsjebefore.k
Qmean.drpBCBCBPsjebefore.p[[x]] <- mean.drpBCBCBPsjebefore.p
Qmean.drpBCBPsjeafter.ds[[x]] <- mean.drpBCBPsjeafter.ds
Qmean.drpBCBPsjeafter.k[[x]] <- mean.drpBCBPsjeafter.k
Qmean.drpBCBPsjeafter.p[[x]] <- mean.drpBCBPsjeafter.p
Qmean.drpBCBPsjebefore.ds[[x]] <- mean.drpBCBPsjebefore.ds
Qmean.drpBCBPsjebefore.k[[x]] <- mean.drpBCBPsjebefore.k
Qmean.drpBCBPsjebefore.p[[x]] <- mean.drpBCBPsjebefore.p
Qmean.drpBCsje.ds[[x]] <- mean.drpBCsje.ds
Qmean.drpBCsje.k[[x]] <- mean.drpBCsje.k
Qmean.drpBCsje.p[[x]] <- mean.drpBCsje.p
Qmeanfield.dendritic.length.post[[x]] <- meanfield.dendritic.length.post
Qmean.migration.dist[[x]] <- mean.migration.dist
Qmean.rel.force[[x]] <- mean.rel.force
Qminimumx[[x]] <- minimumx
QnBCBP[[x]] <- nBCBP
QnbDistBCBPafter[[x]] <- nbDistBCBPafter
QnbDistBCBP.data[[x]] <- nbDistBCBP.data
Qnbins[[x]] <- nbins
QpBCBCBP[[x]] <- pBCBCBP
QpBCBP[[x]] <- pBCBP
Qp_BC_ttest[[x]] <- p_BC_ttest
Qpolygon.data[[x]] <- polygon.data
QPos.new[[x]] <- Pos.new
QRI_BCBPafter[[x]] <- RI_BCBPafter
QRI_BCBP.data[[x]] <- RI_BCBP.data
QRI_BC.data[[x]] <- RI_BC.data
Qsd.dendritic.length[[x]] <- sd.dendritic.length
Qsd.drpBCBCBPsjeafter.ds[[x]] <- sd.drpBCBCBPsjeafter.ds
Qsd.drpBCBPsjeafter.ds[[x]] <- sd.drpBCBPsjeafter.ds
Qsd.drpBCsje.ds[[x]] <- sd.drpBCsje.ds
Qsdfield.dendritic.length.post[[x]] <- sdfield.dendritic.length.post
Qsd.migration.dist[[x]] <- sd.migration.dist
Qsd.rel.force[[x]] <- sd.rel.force
Qsynapses[[x]] <- synapses
QTBCBCBP[[x]] <- TBCBCBP
Qx[[x]] <- x
	

print(paste("Finished with field x =",x))

print("check 8")




#}

print("check 9")

save.image(file=paste("rdata_depth_",commandlineargs[1],"_",commandlineargs[2],"_",commandlineargs[3],"_",commandlineargs[4],"_",commandlineargs[5],"_",commandlineargs[6],"_",dexcBC,"_",dexcBC.sd,"_",dexcBC.trunc,"_",dexcBCBP,"_",dexcBCBP.sd,"_",dexcBCBP.trunc,".RData",sep=""))

##### Make plot of pBCBCBP (cross value) as function of depth
postscript(file=paste("plot_depth_",commandlineargs[1],"_",commandlineargs[2],"_",commandlineargs[3],"_",commandlineargs[4],"_",commandlineargs[5],"_",commandlineargs[6],"_",dexcBC,"_",dexcBC.sd,"_",dexcBC.trunc,"_",dexcBCBP,"_",dexcBCBP.sd,"_",dexcBCBP.trunc,".ps",sep=""))

plot(depths,pBCBCBP, type="b", lwd=1, xlab="Max depth of mobile BCBPs", ylab="Cross p value", main=paste("Sensitivity to depth, field ",field,sep=""))

text((commandlineargs[2]+commandlineargs[3])/2,(max(pBCBCBP)+min(pBCBCBP))/2,paste("Field = ",commandlineargs[1],"\n","Lower depth limit =",commandlineargs[2],"\n","Upper depth limit =",commandlineargs[3],"\n","Depth step size =",commandlineargs[4],"\n","Number of iterations =",commandlineargs[5],"\n","Upper integation limit =",commandlineargs[6]))

dev.off()


#### End of Script ####
date()
proc.time()
