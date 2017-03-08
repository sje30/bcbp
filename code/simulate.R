# This file should contain a function to run m simulations in a loop (calling 'sim' from 'bcbp_sim')

########################################################################
### simulate						                                 ### 
###------------------------------------------------------------------###
### Does 'loopsize' simulations			                             ###
###------------------------------------------------------------------###
### Args:	commandlineargs = commandline arguments	to               ###
###                           'sensitivity_analysis' without         ###
###                           text based args                        ###
###		param1 = first parameter to be tested			             ###
###		param2 = second parameter to be tested			             ###
###		param2.value = current value of param2                       ###
###		AVLROOT = directory which files are read from                ###
########################################################################
simulate <- function(commandlineargs, param1, param2, param2.value, AVLROOT) {
	x <- commandlineargs[1]                # field number
	loopsize <- commandlineargs[2]         # number of iterations in simulation study
	ulimit <- commandlineargs[3]	       # upper limit for calculation of p scores (mu meter)
	low <- commandlineargs[4]              # minimum value of 'param1' to be tested
	high <- commandlineargs[5]             # maximum value of 'param1' to be tested
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
	# if the 2nd parameter to be tested is max_dendr_length, set this to the correct value
	if(!is.na(param2) && param2 == "max_dendr_length") {
		maxdendr <- param2.value
	} else {
		maxdendr <- 44				        # Maxmimal horisontal dendrite length
	}
	Gvert.d.mean <- 10 + 0.33*GINLdepth	# Mean vertical distance from BC pedicles to BCPB perikaryons (Koyama 1992)
	Gvert.d.sd <- GINLdepth*0.11		# Standard deviation of vertical distrance from BC pedicles to BCBP perikaryons
	# if the 2nd parameter to be tested is depth, set this to the correct value
	if(!is.na(param2) && param2 == "depth") {
		vert.d.thr <- param2.value
	} else {
		vert.d.thr <- 3*r+10			    # Threshold for vertical distance between BCBP perikaryons and BC pedicles, if distance is greater there is no BCBP motility
	}
	Gmean.rel.force <- c(0.75,0.90,0.82,0.84,0.97)	# Mean value of relative force (friction/string)
	sd.rel.force <- 0.1			        # Standard deviation of relative force (friction/string)
	bin <- 10				            # Bin size
	nbins <- 10				            # Number of bins
	Gseed <- c(1,2,2,2,2)			    # Seeds for random number generation

	#################################
	### Reading data
	#################################
	exptfile <- sprintf("%s/data/bc_bcbp_f%s.dat", AVLROOT, x)
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

	# set correct values depending on whether the parameter in question is the 2nd
	# parameter to be tested
	if(!is.na(param2) && param2 == "cone_excl") {
		dexcCone <- param2.value
	} else {
		dexcCone <- GdexcCone[x]
	}
	if(!is.na(param2) && param2 == "bc_excl_mean") {
		dexcBC <- param2.value
	} else {
		dexcBC <- GdexcBC[x]
	}
	if(!is.na(param2) && param2 == "bc_excl_sd") {
		dexcBC <- param2.value
	} else {
		dexcBC.sd <- GdexcBC.sd[x]
	}
	if(!is.na(param2) && param2 == "bc_excl_trunc") {
		dexcBC.trunc <- param2.value
	} else {
		dexcBC.trunc <- GdexcBC.trunc[x]
	}
	if(!is.na(param2) && param2 == "bcbp_exl_mean") {
		dexcBCBP <- param2.value
	} else {
		dexcBCBP <- GdexcBCBP[x]
	}
	if(!is.na(param2) && param2 == "bcbp_excl_sd") {
		dexcBCBP.sd <- param2.value
	} else {
		dexcBCBP.sd <- GdexcBCBP.sd[x]
	}
	if(!is.na(param2) && param2 == "bcbp_excl_trunc") {
		dexcBCBP.trunc <- param2.value
	} else {
		dexcBCBP.trunc <- GdexcBCBP.trunc[x]
	}
	INLdepth <- GINLdepth[x]
	vert.d.mean <- Gvert.d.mean[x]
	vert.d.sd <- Gvert.d.sd[x]
	if(!is.na(param2) && param2 == "rel_force_mean") {
		mean.rel.force <- param2.value
	} else {
		mean.rel.force <- Gmean.rel.force[x]
	}
	polygon <- matrix( c(0,0,0,ly,lx,ly,lx,0,0,0), nrow=5, ncol=2, byrow=T)
	polygon.data <- matrix( c(minimumx,minimumy,minimumx,maximumy,maximumx,maximumy,maximumx,minimumy,minimumx,minimumy), nrow=5, ncol=2, byrow=T)

	sens <- seq(low, high, incr) # vector of values of param1
	n <- length(sens)
	
	simOut <- list()

	for(m in 1:n) {
		#################################
		### Overwrite default values of the parameters we are testing
		#################################
		if(param1 == "depth") {
			vert.d.thr <- sens[m]
		} else if(param1 == "cone_excl") {
			dexcCone <- sens[m]
		} else if(param1 == "bc_excl_mean") {
			dexcBC <- sens[m]
		} else if(param1 == "bc_excl_sd") {
			dexcBC.sd <- sens[m]
		} else if(param1 == "bc_excl_trunc") {
			dexcBC.trunc <- sens[m]
		} else if(param1 == "bcbp_excl_mean") {
			dexcBCBP <- sens[m]
		} else if(param1 == "bcbp_excl_sd") {
			dexcBCBP.sd <- sens[m]
		} else if(param1 == "bcbp_excl_trunc") {
			dexcBCBP.trunc <- sens[m]
		} else if(param1 == "max_dendr_length") {
			maxdendr <- sens[m]
		} else if(param1 == "rel_force_mean") {
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
		simOut[[m]] <- sim(coneNumbers, dexcValues, conePos, bin, nbins, polygon, maxdendr, vert.d.values, rel.force.values, lx, ly, loopsize)
		
		print(paste("Outer loop m =", m, "completed"))
	}
	
	# add these to output because they're needed later
	simOut[[(n+1)]] <- data1
	simOut[[(n+2)]] <- data2
	simOut[[(n+3)]] <- bin
	simOut[[(n+4)]] <- nbins
	simOut[[(n+5)]] <- polygon.data
	simOut[[(n+6)]] <- nBC
	simOut[[(n+7)]] <- nBCBP
	
	return(simOut)
}
