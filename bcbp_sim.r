# This file should contain a function which calculates simulations of movement
# of blue cone bipolar cells in retina based on the model in the paper.

##sim <- function(nC, nBC, nBCBP, dexcBC, dexcBC.sd, dexcBC.trunc, dexcBCBP, dexcBCBP.sd, dexcBCBP.trunc,
##               conePos, nbins, bin, polygon, maxdendr,
##                vert.d.mean, vert.d.sd, vert.d.thr, mean.rel.force, sd.rel.force) {

sim <- function(coneNumbers, dexcValues, conePos, bin, nbins, polygon, maxdendr, vert.d.values, rel.force.values, lx, ly, loopsize) {
	#############################
	### Unpack parameters
	#############################
	nC <- coneNumbers[1]
	nBC <- coneNumbers[2]
	nBCBP <- coneNumbers[3]
	
	dexcBC <- dexcValues[1]
	dexcBC.sd <- dexcValues[2]
	dexcBC.trunc <- dexcValues[3]
	dexcBCBP <- dexcValues[4]
	dexcBCBP.sd <- dexcValues[5]
	dexcBCBP.trunc <- dexcValues[6]
	
	vert.d.mean <- vert.d.values[1]
	vert.d.sd <- vert.d.values[2]
	vert.d.thr <- vert.d.values[3]
	
	mean.rel.force <- rel.force.values[1]
	sd.rel.force <- rel.force.values[2]
	
	#############################
	### Initialise output variables
	#############################
	drpBCsje <- list()							                    # List with drp, for BCs, one entry for each simulation (all, using sjedrp package)
	
	drpBCBPsjeafter <- list() 						                # List with drp, for BCBPs, one entry for each simulation (all, using sjedrp package)
	drpBCBPsjebefore <- list()						                # List with drp, for BCBPs, one entry for each simulation (all, using sjedrp package)
	
	drpBCBCBPsjeafter <- list()						                # List with drp, for BC to BCBPs, one entry for each simulation (all, using sjedrp package)
	drpBCBCBPsjebefore <- list()						            # List with drp, for BC to BCBPs, one entry for each simulation (all, using sjedrp package)
	
	LBC <- matrix(NA, nrow=loopsize, ncol=(nbins*bin))			    # Matrix with L function for BCs, one row for each simulation
	LBCBPafter <- matrix(NA, nrow=loopsize, ncol=(nbins*bin))		# Matrix with L function for BCBPs, one row for each simulation
	LBCBPbefore <- matrix(NA, nrow=loopsize, ncol=(nbins*bin))		# Matrix with L function for BCBPs, one row for each simulation
	LBCBCBPafter <- matrix(NA, nrow=loopsize, ncol=(nbins*bin))		# Matrix with L function for cross between BCs and BCBPs, one row for each simulation
	LBCBCBPbefore <- matrix(NA, nrow=loopsize, ncol=(nbins*bin))	# Matrix with L function for cross between BCs and BCBPs, one row for each simulation
	
	GBC <- matrix(NA, nrow=loopsize, ncol=(nbins*bin))			    # Matrix with G function for BCs, one row for each simulation
	GBCBPbefore <- matrix(NA, nrow=loopsize, ncol=(nbins*bin))		# Matrix with G function for BCBPs, one row for each simulation
	GBCBPafter <- matrix(NA, nrow=loopsize, ncol=(nbins*bin))		# Matrix with G function for BCBPs, one row for each simulation
	
	synapses <- vector()							                # Number of synapses in each simulation
	movedBCBPs <- vector()							                # Number of BCBPS which have moved in each simulation
	mean.migration.dist <- vector()						            # Average migration distance among BCBPs that move
	sd.migration.dist <- vector()						            # Standad deviation of migration distance among BCBPs that move
	
	dendritic.length <- matrix(NA, ncol=nBCBP, nrow=loopsize)
	dendritic.length.post <- matrix(NA, ncol=nBCBP, nrow=loopsize)
	mean.dendritic.length <- vector()
	mean.dendritic.length.post <- vector()
	sd.dendritic.length <- vector()
	sd.dendritic.length.post <- vector()
	
	RI_BC <- vector()							    # Regularity index for blue cones
	RI_BCBPbefore <- vector()						# Regularity index for BCBP before migration
	RI_BCBPafter <- vector()						# Regularity index for BCBP after migration

	for(k in 1:loopsize) {
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
		Pos.new <- migration(nBCBP,vert.d.mean,vert.d.sd,vert.d.thr,mean.rel.force,sd.rel.force,dexcBCBP.trunc,dendr.length,links,BCBPPos, bcPos)
		
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
	}
	
	simOut <- list(drpBCsje, drpBCBPsjeafter, drpBCBPsjebefore, drpBCBCBPsjeafter, drpBCBCBPsjebefore, LBC,
	            LBCBPafter, LBCBPbefore, LBCBCBPafter, LBCBCBPbefore, GBC, GBCBPbefore, GBCBPafter, synapses,
	            movedBCBPs, mean.migration.dist, sd.migration.dist, dendritic.length, dendritic.length.post,
	            mean.dendritic.length, mean.dendritic.length.post, sd.dendritic.length, sd.dendritic.length.post,
	            RI_BC, RI_BCBPbefore, RI_BCBPafter)
	
	return(simOut)
}
