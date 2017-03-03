# This file should take input stating which parameter to test for sensitivity
# and also other relevant parameters and do the simulations and the sensitivity
# test

##AVLROOT <- Sys.getenv("AVLROOT")
AVLROOT <- "."

#################################
### Sourcing functions
#################################
## SJE:
extra.functions <- sprintf("%s/functions_mimicsthesis2.R", AVLROOT)
source(extra.functions)

more.functions <- sprintf("%s/bcbp_sim.R", AVLROOT)
source(more.functions)

#################################
### Loading packages
#################################
library(stats)
library(class)
library(graphics)
library(splancs)
library(sjedrp)

#### Importing command line arguments ###
#### The format for command line arguments is:
#### field number, number of loops, upper limit for calculating p scores,
#### parameter to be tested for sensitivity, minimum (of that parameter), maximum, increment

# get commandline arguments as strings
commandlineargs <- commandArgs(TRUE)

analysis <- function(commandlineargs, AVLROOT) {
	sens.param <- commandlineargs[4]       # which parameter is to be tested for sensitivity (this will be a string)
	# can be: 'depth', 'cone_excl', 'bc_excl_mean', 'bc_excl_sd', 'bc_excl_trunc', 'bcbp_excl_mean',
	#         'bcbp_excl_sd', 'bcbp_excl_trunc', 'max_dendr_length', 'rel_force_mean'
	
	# convert all numeric commandline arguments to numeric values
	commandlineargs <- as.numeric(commandlineargs[-4])
	
	x <- commandlineargs[1]                # field number
	loopsize <- commandlineargs[2]         # number of iterations in simulation study
	ulimit <- commandlineargs[3]	       # upper limit for calculation of p scores (mu meter)
	low <- commandlineargs[4]              # minimum value of 'sens.param' to be tested
	high <- commandlineargs[5]             # maximum value of 'sens.param' to be tested
	incr <- commandlineargs[6]             # increment

	#################################
	### Global parameters
	#################################
	a <- 0.10				            # Fraction of cone which are blue cones (Roorda2001)
	r <- 4					            # BCBP radius (in mid peripheral retina) (Koyama 1992)
	GdexcCone <- c(12,12,12,12,12)		# Radius of cone exclusion zone (mu metre)
	GdexcBC <- c(30,38.3,28,33,36)		# Mean radius of blue cone exclusion zone (mu metre)
	GdexcBC.sd <- c(6,3.9,3,7,7)		# Standard deviation of blue cone exclusion zone (mu metre)
	GdexcBC.trunc <- c(17,32,23,18,23)	# Blue clone exclusion zone truncation, based on min observed nearest neighbour distance
	GdexcBCBP <- c(18,22,19,18,22)		# Radius of BCBP exclusion zone (mu metre) (this parameter is important for the even distribution of links), used to be 8
	GdexcBCBP.sd <- c(5,6,5,5,8)		# BCBP exclusion zone standard deviation
	GdexcBCBP.trunc <- c(8.2,7.2,5.4,5.7,5.4)		# BCBP exclusion zone truncation, based on min. observed nearest neighbour distance
	GINLdepth <- c(40,20,40,30,15)		# Depth of inner nuclear layer (Kouyama)
	maxdendr <- 44				        # Maxmimal horisontal dendrite length
	Gvert.d.mean <- 10 + 0.33*GINLdepth	# Mean vertical distance from BC pedicles to BCPB perikaryons (Koyama 1992)
	Gvert.d.sd <- GINLdepth*0.11		# Standard deviation of vertical distrance from BC pedicles to BCBP perikaryons
	vert.d.thr <- 3*r+10			    # Threshold for vertical distance between BCBP perikaryons and BC pedicles, if distance is greater there is no BCBP motility
	Gmean.rel.force <- c(0.75,0.90,0.82,0.84,0.97)	# Mean value of relative force (friction/string)
	sd.rel.force <- 0.1			        # Standard deviation of relative force (friction/string)
	bin <- 10				            # Bin size
	nbins <- 10				            # Number of bins
	Gseed <- c(1,2,2,2,2)			    # Seeds for random number generation

	#################################
	### Reading data
	#################################
	exptfile <- sprintf("%s/data/bc_bcbp_f%d.dat", AVLROOT, x)
	##exptfile <- paste("data/bc_bcbp_f",x,".dat", sep="")
	data <- read.table(exptfile)
	#data <- data[sample(1:length(data[,1]),80),] ## This is temporary, only here to make the script run faster
	data1 <- data[which(data[,3] == 1),c(1,2)]
	data2 <- data[which(data[,3] == 2),c(1,2)]
	# convert data1 and data2 to matrices so lFunc works on it
	data1 <- as.matrix(data1)
	data2 <- as.matrix(data2)

	#################################
	### Setting field parameters (these may be overwritten if they relate to the parameter to be tested)
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

	cat("Here starts the", sens.param, "optimization for field",x,"\n",file=paste("LogField",x,"_",sens.param,".txt",sep=""),append=FALSE)

	sens <- seq(low, high, incr) # vector of values of sens.param
	n <- length(sens)

	pBC = matrix(NA, nrow = n, ncol = 1)
	pBCBP = matrix(NA, nrow = n, ncol = 1)
	pBCBCBP = matrix(NA, nrow = n, ncol = 1)

	for(m in 1:n) {
		#################################
		### Overwrite default values of the parameters we are testing
		#################################
		if(sens.param == "depth") {
			vert.d.thr <- sens[m]
		} else if(sens.param == "cone_excl") {
			dexcCone <- sens[m]
		} else if(sens.param == "bc_excl_mean") {
			dexcBC <- sens[m]
		} else if(sens.param == "bc_excl_sd") {
			dexcBC.sd <- sens[m]
		} else if(sens.param == "bc_excl_trunc") {
			dexcBC.trunc <- sens[m]
		} else if(sens.param == "bcbp_excl_mean") {
			dexcBCBP <- sens[m]
		} else if(sens.param == "bcbp_excl_sd") {
			dexcBCBP.sd <- sens[m]
		} else if(sens.param == "bcbp_excl_trunc") {
			dexcBCBP.trunc <- sens[m]
		} else if(sens.param == "max_dendr_length") {
			maxdendr <- sens[m]
		} else if(sens.param == "rel_force_mean") {
			mean.rel.force <- sens[m]
		}

		#################################
		### Setting seed
		#################################
		set.seed(Gseed[x])

		#################################
		### Generating cone mosaic
		#################################
		conePos <- position(nC, lx, ly, dexcCone, 0, dexcCone)    # lx and ly are the field dimensions, nC is the number of cones and
		                                                          # dexcCone is the cone exclusion zone (which is also the truncation in
		                                                          # this case, where there is no standard deviation, the exclusion zone is constant

		################################
		### Packing variables for use in function 'sim'
		################################
		coneNumbers <- c(nC, nBC, nBCBP)
		dexcValues <- c(dexcBC, dexcBC.sd, dexcBC.trunc, dexcBCBP, dexcBCBP.sd, dexcBCBP.trunc)
		vert.d.values <- c(vert.d.mean, vert.d.sd, vert.d.thr)
		rel.force.values <- c(mean.rel.force, sd.rel.force)

		################################
		### call simulation function with 'loopsize' loops
		################################
		simOut <- sim(coneNumbers, dexcValues, conePos, bin, nbins, polygon, maxdendr, vert.d.values, rel.force.values, lx, ly, loopsize)

		################################
		### Unpacking output from 'sim'
		################################
		drpBCsje <- simOut[[1]]
		drpBCBPsjeafter <- simOut[[2]]
		drpBCBPsjebefore <- simOut[[3]]
		drpBCBCBPsjeafter <- simOut[[4]]
		drpBCBCBPsjebefore <- simOut[[5]]

		LBC <- simOut[[6]]
		LBCBPafter <- simOut[[7]]
		LBCBPbefore <- simOut[[8]]
		LBCBCBPafter <- simOut[[9]]
		LBCBCBPbefore <- simOut[[10]]

		GBC <- simOut[[11]]
		GBCBPbefore <- simOut[[12]]
		GBCBPafter <- simOut[[13]]

		synapses <- simOut[[14]]
	    movedBCBPs <- simOut[[15]]
	    mean.migration.dist <- simOut[[16]]
	    sd.migration.dist <- simOut[[17]]

	    dendritic.length <- simOut[[18]]
	    dendritic.length.post <- simOut[[19]]
	    mean.dendritic.length <- simOut[[20]]
	    mean.dendritic.length.post <- simOut[[21]]
	    sd.dendritic.length <- simOut[[22]]
	    sd.dendritic.length.post <- simOut[[23]]

	    RI_BC <- simOut[[24]]
	    RI_BCBPbefore <- simOut[[25]]
	    RI_BCBPafter <- simOut[[26]]

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

	    #print("check 1")

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

	    #print("check 2")
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
	    #print("check 3")

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
	    #print("check 4")

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

	    #print("check 6")

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
	    cat("Here endeth", sens.param, ":",sens[m],"\n",file=paste("LogField",x,"_",sens.param,".txt",sep=""),append=TRUE)

	    ## Making plots to get a visual impression of the output


	    postscript(file=paste("field",x,"_", sens.param, "_", sens[m], "_",loopsize,"_",ulimit,".ps",sep=""))	#   postscript(file=paste("field1,",y,",",z,".pdf",sep="")) #
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

	cat("Here ends the", sens.param, "optimization for field",x,"\n",file=paste("LogField",x,"_",sens.param,".txt",sep=""),append=TRUE)

	# Collecting all p values in one matrix
	pvalues <- cbind(pBC,pBCBP,pBCBCBP)
	colnames(pvalues) <- c("pBC","pBCBP","pBCBCBP")
	rownames(pvalues) <- sens

	##############
	### Redirecting p values for BC L function to a file
	#############
	sink(paste("pvalues_", sens.param, "_",x,"_",low,"_",high,"_",incr,"_",loopsize,"_",ulimit,"_",dexcBC,"_",dexcBC.sd,"_",dexcBC.trunc,"_",dexcBCBP,"_",dexcBCBP.sd,"_",dexcBCBP.trunc,".txt",sep=""))

	#print(paste("p values for", sens.param))
	print(pvalues)

	sink()

	cat("p values printed to file",file=paste("LogField",x,"_",sens.param,".txt",sep=""),append=TRUE)

	save.image(file=paste("rdata_", sens.param, "_",x,"_",low,"_",high,"_",incr,"_",loopsize,"_",ulimit,"_",dexcBC,"_",dexcBC.sd,"_",dexcBC.trunc,"_",dexcBCBP,"_",dexcBCBP.sd,"_",dexcBCBP.trunc,".RData",sep=""))

	##### Make plot of pBCBCBP (cross value) as function of sens.param
	postscript(file=paste("plot_", sens.param, "_",x,"_",low,"_",high,"_",incr,"_",loopsize,"_",ulimit,"_",dexcBC,"_",dexcBC.sd,"_",dexcBC.trunc,"_",dexcBCBP,"_",dexcBCBP.sd,"_",dexcBCBP.trunc,".ps",sep=""))

	# set labels according to parameter tested
	if(sens.param == "depth") {
		xlabel <- "Max depth of mobile BCBPs"
		title <- "depth"
	} else if(sens.param == "cone_excl") {
		xlabel <- "Cone exclusion zone"
		title <- paste("Sensitivity to cone exclusion zone, field ",x,sep="")
	} else if(sens.param == "bc_excl_mean") {
		xlabel <- "Blue Cone exclusion zone mean"
		title <- paste("Sensitivity to blue cone exclusion zone mean, field ",x,sep="")
	} else if(sens.param == "bc_excl_sd") {
		xlabel <- "Blue Cone exclusion zone sd"
		title <- paste("Sensitivity to blue cone exclusion zone sd, field ",x,sep="")
	} else if(sens.param == "bc_excl_trunc") {
		xlabel <- "Blue Cone exclusion zone truncation"
		title <- paste("Sensitivity to blue cone exclusion zone truncation, field ",x,sep="")
	} else if(sens.param == "bcbp_excl_mean") {
		xlabel <- "Blue Cone Bipoler exclusion zone mean"
		title <- paste("Sensitivity to blue cone bipoler exclusion zone mean, field ",x,sep="")
	} else if(sens.param == "bcbp_excl_sd") {
		xlabel <- "Blue Cone Bipoler exclusion zone sd"
		title <- paste("Sensitivity to blue cone bipoler exclusion zone sd, field ",x,sep="")
	} else if(sens.param == "bcbp_excl_trunc") {
		xlabel <- "Blue Cone Bipoler exclusion zone truncation"
		title <- paste("Sensitivity to blue cone bipoler exclusion zone truncation, field ",x,sep="")
	} else if(sens.param == "max_dendr_length") {
		xlabel <- "Maximum dendrite length"
		title <- paste("Sensitivity to maximum dendrite length, field ",x,sep="")
	} else if(sens.param == "rel_force_mean") {
		xlabel <- "Mean relative force"
		title <- paste("Sensitivity to mean relative force, field ",x,sep="")
	}

	plot(sens, pBCBCBP, type="b", lwd=1, xlab=xlabel, ylab="Cross p value", main=title)

	text((low+high)/2,(max(pBCBCBP)+min(pBCBCBP))/2,paste("Field = ",x,"\n","Lower limit =",low,"\n","Upper limit =",high,"\n","Step size =",incr,"\n","Number of iterations =",loopsize,"\n","Upper integation limit =",ulimit))

	dev.off()
}

###########################
## Run simulations
###########################
analysis(commandlineargs, AVLROOT)












