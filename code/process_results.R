# This file should a function which processes the results of simulations
# and 2 functions which call it - one for 1 parameter, one for 2

########################################################################
### process_results						                             ### 
###------------------------------------------------------------------###
### Gets statistics from the data and does lFunc and DRP plots	     ###
###------------------------------------------------------------------###
### Args:    simOut = results of simulations                         ###
###    commandlineargs = commandline arguments to                    ###
###                       'sensitivity_analysis' without             ###
###                       text based args                            ###
###		param1 = first parameter to be tested			             ###
###		param2 = second parameter to be tested			             ###
###		param2.value = current value of param2                       ###
###		plot_l_drp = whether or not to do the lFunc and DRP plots    ###
########################################################################
process_results <- function(simOut, commandlineargs, param1, param2, param2.value, plot_l_drp) {
	################################
	### Get command line arguments
	################################
	x <- commandlineargs[1]                # field number
	loopsize <- commandlineargs[2]         # number of iterations in simulation study
	ulimit <- commandlineargs[3]	       # upper limit for calculation of p scores (mu meter)
	low <- commandlineargs[4]              # minimum value of 'param1' to be tested
	high <- commandlineargs[5]             # maximum value of 'param1' to be tested
	incr <- commandlineargs[6]             # increment
	
	range1 <- seq(low, high, incr)
	
	n <- length(simOut)
	
	# for some reason these need to be converted back into numeric form
	# or, in the case of data1 and 2 and polygon.data, to data frames then matrices
	data1 <- simOut[n-6]
	data1 <- as.matrix(as.data.frame(data1))
	data2 <- simOut[n-5]
	data2 <- as.matrix(as.data.frame(data2))
	bin <- as.numeric(simOut[n-4])
	nbins <- as.numeric(simOut[n-3])
	polygon.data <- as.matrix(as.data.frame(simOut[n-2]))
	nBC <- as.numeric(simOut[n-1])
	nBCBP <- as.numeric(simOut[n])
	
	pBC = matrix(NA, nrow = (n-7), ncol = 1)
	pBCBP = matrix(NA, nrow = (n-7), ncol = 1)
	pBCBCBP = matrix(NA, nrow = (n-7), ncol = 1)

	for(m in 1:(n-7)) {
		################################
		### Unpacking output from 'sim'
		################################
		sim_current <- simOut[[m]]
		drpBCsje <- sim_current[[1]]
		drpBCBPsjeafter <- sim_current[[2]]
		drpBCBPsjebefore <- sim_current[[3]]
		drpBCBCBPsjeafter <- sim_current[[4]]
		drpBCBCBPsjebefore <- sim_current[[5]]

		LBC <- sim_current[[6]]
		LBCBPafter <- sim_current[[7]]
		LBCBPbefore <- sim_current[[8]]
		LBCBCBPafter <- sim_current[[9]]
		LBCBCBPbefore <- sim_current[[10]]

		GBC <- sim_current[[11]]
		GBCBPbefore <- sim_current[[12]]
		GBCBPafter <- sim_current[[13]]

		synapses <- sim_current[[14]]
	    movedBCBPs <- sim_current[[15]]
	    mean.migration.dist <- sim_current[[16]]
	    sd.migration.dist <- sim_current[[17]]

	    dendritic.length <- sim_current[[18]]
	    dendritic.length.post <- sim_current[[19]]
	    mean.dendritic.length <- sim_current[[20]]
	    mean.dendritic.length.post <- sim_current[[21]]
	    sd.dendritic.length <- sim_current[[22]]
	    sd.dendritic.length.post <- sim_current[[23]]

	    RI_BC <- sim_current[[24]]
	    RI_BCBPbefore <- sim_current[[25]]
	    RI_BCBPafter <- sim_current[[26]]
	    

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
		
		# Calculating the envelope of the L function for the simulated BC to BCBPs, before
		LBCBCBPbeforeenvelope <- apply(LBCBCBPbefore, 2, function(x) {quantile(x, probs=c(0.025,0.975))})
		
		# Calculating the envelope of the L function for the simulated BC to BCBPs, after
		LBCBCBPafterenvelope <- apply(LBCBCBPafter, 2, function(x) {quantile(x, probs=c(0.025,0.975))})
		
		### Calculating the envelope of the G function for the simulated BCs
		##GBCenvelope <- apply(GBC, 2, function(x) {quantile(x, probs=c(0.025,0.975))})
		##
		### Calculating the envelope of the G function for the simulated BCBPs, before
		##GBCBPbeforeenvelope <- apply(GBCBPbefore, 2, function(x) {quantile(x, probs=c(0.025,0.975))})
		##
		### Calculating the envelope of the G function for the simulated BCBPs, after
		##GBCBPafterenvelope <- apply(GBCBPafter, 2, function(x) {quantile(x, probs=c(0.025,0.975))})
		
		
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
		
		##pBCBCBP
		
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
		
	
		#####################
		### L-function and DRP plots (if wanted)
		#####################
	
		if(plot_l_drp == TRUE) {
			# name files depending on whether there is a second parameter
			if(is.na(param2)) {
				pdf(file=paste("field",x,"_", param1,"_", range1[m],"_",loopsize,"_",ulimit,".pdf",sep=""))	#   pdf(file=paste("field1,",y,",",z,".pdf",sep="")) #
			} else {
				pdf(file=paste("field",x,"_", param1,"_", range1[m],"_",param2,"_",param2.value,"_",loopsize,"_",ulimit,".pdf",sep=""))	#   pdf(file=paste("field1,",y,",",z,".pdf",sep="")) #
			}
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
	}
	
	# Collecting all p values in one matrix
	pvalues <- cbind(pBC,pBCBP,pBCBCBP)
	
	return(pvalues)
}

########################################################################
### process1param						                             ### 
###------------------------------------------------------------------###
### Processes data for 1 parameter incl. plotting pvalues			 ###
###------------------------------------------------------------------###
### Args:    simOut = results of simulations                         ###
###    commandlineargs = commandline arguments to                    ###
###                       'sensitivity_analysis' without             ###
###                       text based args                            ###
###		param = parameter to be tested			                     ###
###		plot_l_drp = whether or not to do the lFunc and DRP plots    ###
########################################################################
process1param <- function(simOut, commandlineagrs, param, plot_l_drp) {
	################################
	### Get command line arguments
	################################
	x <- commandlineargs[1]                # field number
	loopsize <- commandlineargs[2]         # number of iterations in simulation study
	ulimit <- commandlineargs[3]	       # upper limit for calculation of p scores (mu meter)
	low <- commandlineargs[4]              # minimum value of 'param' to be tested
	high <- commandlineargs[5]             # maximum value of 'param' to be tested
	incr <- commandlineargs[6]             # increment
	
	range1 <- seq(low, high, incr)
	
	# call process results and assign it to pvalues
	# NOTE: not a pure function - also does plots (L-function and DRP) if 'plot_l_drp' = TRUE
	pvalues <- process_results(simOut, commandlineargs, param, NA, NA, plot_l_drp)
	colnames(pvalues) <- c("pBC","pBCBP","pBCBCBP")
	rownames(pvalues) <- range1

	#####################
	### Plots of the cross p-value and record of all p-values
	#####################

	##############
	### Redirecting p values for BC L function to a file
	#############
	sink(paste("pvalues_", param, "_",x,"_",low,"_",high,"_",incr,"_",loopsize,"_",ulimit,".txt",sep=""))

	#print(paste("p values for", param))
	print(pvalues)

	sink()
	
	# change pvalues to data frame so its columns can be accessed by '$'
	pvalues <- as.data.frame(pvalues)

	#####################
	### Plot the cross p-values
	#####################
	plot_lines(pvalues, commandlineargs, param)
}

########################################################################
### process2params						                             ### 
###------------------------------------------------------------------###
### Processes data for 2 parameters incl. levelplot for pvalues	     ###
###------------------------------------------------------------------###
### Args:    simOut = results of simulations                         ###
###    commandlineargs = commandline arguments to                    ###
###                       'sensitivity_analysis' without             ###
###                       text based args                            ###
###		param1 = first parameter to be tested			             ###
###     param2 = second parameter to be tested                       ###
###		plot_l_drp = whether or not to do the lFunc and DRP plots    ###
########################################################################
process2params <- function(simOut, commandlineargs, param1, param2, plot_l_drp) {
	################################
	### Get command line arguments
	################################
	x <- commandlineargs[1]                # field number
	loopsize <- commandlineargs[2]         # number of iterations in simulation study
	ulimit <- commandlineargs[3]	       # upper limit for calculation of p scores (mu meter)
	low1 <- commandlineargs[4]              # minimum value of 'param1' to be tested
	high1 <- commandlineargs[5]             # maximum value of 'param1' to be tested
	incr1 <- commandlineargs[6]             # increment
	
	low2 <- commandlineargs[7]
	high2 <- commandlineargs[8]
	incr2 <- commandlineargs[9]
	
	range1 <- seq(low1, high1, incr1)

	range2 <- seq(low2, high2, incr2)
	n <- length(range2)
	
	pvalues <- matrix()
	
	################################
	### Get p-values for all variations of the 2 parameters
	################################
	for(k in 1:n) {
		# call process results and assign it to pvalues
		# NOTE: not a pure function - also does plots (L-function and DRP) if 'plot_l_drp' = TRUE
		pvalues_current <- process_results(simOut[[k]], commandlineargs, param1, param2, range2[k], plot_l_drp)
		m <- length(range1)
		# add columns in 'pvalues_current' to specify the values of 'param1' and 'param2'
		pvalues_current <- cbind(range1, rep(range2[k], m), pvalues_current)

		# add new pvalues to 'pvalues'
		if(all(is.na(pvalues))) {
			pvalues <- pvalues_current
		} else {
			pvalues <- rbind(pvalues, pvalues_current)
		}
	}
	
	# change pvalues to data frame so its columns can be accessed by '$'
	pvalues <- as.data.frame(pvalues)
	
	colnames(pvalues) <- c("param1_value", "param2_value", "pBC", "PBCBP", "pBCBCBP")
	
	##############
	### Redirecting p values for BC L function to a file
	#############
	sink(paste("pvalues_",param1,"_",param2,"_",x,"_",low1,"_",high1,"_",incr1,"_",low2,"_",high2,"_",incr2,"_",loopsize,"_",ulimit,".txt",sep=""))

	#print(paste("p values for", param))
	print(pvalues)

	sink()

	
	################################
	### Make a level plot
	################################
	plot_level(pvalues, commandlineargs, param1, param2)
}
		


