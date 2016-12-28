###################################################################
### Function definitions for friction/string force model	###
### Andreas Vejen Lønborg					###
### Updated: 2010.02.16						###
###								###
###################################################################

###################################################################
### position							### 
###-------------------------------------------------------------###
### Generation of mosaic of cells				###
###-------------------------------------------------------------###
### Args:	n = number of cells				###
###		lx = x dimension of field			###
###		ly = y dimension of field			###
###		dexc.mean = mean exclusion zone			###
###		dexc.sd = standard deviation of exclusion zone	###
###		dexc.trunc = truncation of exclusion zone	###
###################################################################

position <- function(n,lx,ly,dexc.mean,dexc.sd,dexc.trunc) {

# Positioning cells randomly in area
pos <- matrix(c(runif(n,0,lx),runif(n,0,ly)),ncol=2)

# Finding exclusion zones
if (dexc.sd == 0 ) {
	dexc.random <- rep(dexc.mean,n)
}
else {
	dexc.random <- rnorm(n,dexc.mean,dexc.sd)
}


# Repositioning cells such that there is no overlap
for (i in 1:n) {
	nbDistance <- 0
	cl <- (1:n)[-i]
	while (nbDistance < max(dexc.random[i],dexc.trunc)) {
		pos[i,] <- c(runif(1,0,lx),runif(1,0,ly))
		nb <- as.integer(as.character(knn(pos[-i,],pos[i,], cl, k=1)))
		nbDistance <- sqrt((pos[i,1]-pos[nb,1])^2+(pos[i,2]-pos[nb,2])^2)
	}
}
return(pos)
}




###################################################################
### bcDiff							### 
###-------------------------------------------------------------###
### Step 2: Differentiating the cones into blue cones and	###
###	    non-blue cones					###
###-------------------------------------------------------------###
### Args:	nC = number of cones				###
###		nBC = number of blue cones			###
###		dexcBC = mean blue cone exlusion zone		###
###		dexcBC.sd = standard deviation of BC excl. zone ###
###################################################################

bcDiff <- function(nC,nBC,dexcBC,dexcBC.sd,dexcBC.trunc,conePos) {

# Selecting the blue cones from the population of cones
BCindex <- sample(1:nC,nBC)

# Assigning an individual BC exclusion zone to each cone
random <- rnorm(nC,dexcBC,dexcBC.sd)

# Repositioning the blue cones such that blue cone specific exclusion zone is satisfied
for (i in 1:nBC) {
	nbDist <- 0
	cl <- (1:nBC)[-i]
	while (nbDist < max(random[BCindex[i]],dexcBC.trunc)) {
		BCindex[i] <- sample((1:nC)[-BCindex],1)
		nb <- BCindex[as.integer(as.character(knn(conePos[BCindex[-i],],conePos[BCindex[i],], cl, k=1)))]
		nbDist <- sqrt((conePos[BCindex[i],1]-conePos[nb,1])^2+(conePos[BCindex[i],2]-conePos[nb,2])^2)
	}
}
return(BCindex)
}


###################################################################
### lFunc							### 
###-------------------------------------------------------------###
### Calculate L function given point coordinates		###
###-------------------------------------------------------------###
### Args:	pos = x y coordinates				###
###		polygon = coordinates of corners of the field	###
###		nbins = number of bins in DRP			###
###		bin = binsize in DRP				###
###################################################################

lFunc <- function(pos,polygon,nbins,bin) {
K <- khat(pos, polygon, seq(1:(nbins*bin)) )
L <- sqrt(K/pi)
return(L)
}

###################################################################
### lFuncCross							### 
###-------------------------------------------------------------###
### Calculate cross L function given point coordinates		###
###-------------------------------------------------------------###
### Args:	pos1 = x y coordinates of cell group 1		###
###		pos2 = x y coordinates of cell group 2		###
###		polygon = coordinates of corners of the field	###
###		nbins = number of bins in DRP			###
###		bin = binsize in DRP				###
###################################################################

lFuncCross <- function(BCpos,BCBPpos,polygon,nbins,bin) {
K <- k12hat(BCpos,BCBPpos,polygon,seq(1:(nbins*bin)) )
L <- sqrt(K/pi)
return(L)
}


###################################################################
### nbDist							### 
###-------------------------------------------------------------###
### Calculate distance to nearest neighbour			###
###-------------------------------------------------------------###
### Args:	Pos = coordinates				###
###################################################################

#nbDist <- function(Pos) {

#test <- as.matrix(dist(Pos))
#nbDist <- apply(test,1,function(x) { x[order(x)[2]] } )
#return(nbDist)
#}

nbDist <- function(Pos) {

nbDist <- vector()
nb <- vector()
n <- dim(Pos)[1]
for (i in 1:n) {
	cl <- (1:n)[-i]
	nb[i] <- as.integer(as.character(knn(Pos[-i,],Pos[i,], cl, k=1)))
	nbDist[i] <- sqrt((Pos[i,1]-Pos[nb[i],1])^2+(Pos[i,2]-Pos[nb[i],2])^2)
}
return(nbDist)
}


## This function shall be deleted
nbDistBCBP <- function(Pos) {

test <- as.matrix(dist(BCBPPos))
diag(test) <- 10000
nbBCBPbefore <- apply(test,1,function(x) {as.integer(names(sort(x)[1]))})
nbDistBCBPbefore <- sqrt((BCBPPos[,1]-BCBPPos[nbBCBPbefore,1])^2+(BCBPPos[,2]-BCBPPos[nbBCBPbefore,2])^2)
return(nbDistBCBPbefore)
}
###




###################################################################
### gFunc							### 
###-------------------------------------------------------------###
### Calculate G function given point coordinates		###
###-------------------------------------------------------------###
### Args:	nbDist = distance to neighbour			###
###		bin = binsize in DRP				###
###		nbins = number of bins in DRP			###
###		n = number of cells				###
###################################################################

gFunc <- function(nbDist,bin,nbins,n) {

# Calculating G function
G <- vector()
for (i in 1:(bin*nbins)) {
	G[i] <- sum(nbDist <= i)
}
G <- G/n

return(G)
}

###################################################################
### RI								### 
###-------------------------------------------------------------###
### Calculate regularity index					###
###-------------------------------------------------------------###
### Args:	pos = x y coordinates				###
###		nbDistBC = distance to neighbour for BCs	###
###		bin = binsize in DRP				###
###		nbins = number of bins in DRP			###
###		nBC = number of blue cones			###
###################################################################


regInd <- function(nbDistBC) {
RI <- mean(nbDistBC)/sd(nbDistBC)
return(RI)
}







###################################################################
### connections							### 
###-------------------------------------------------------------###
### For each BCBP, the connection BC is determined		###
###-------------------------------------------------------------###
### Args:	pos = x y coordinates				###
###		nbDistBC = distance to neighbour for BCs	###
###		bin = binsize in DRP				###
###		nbins = number of bins in DRP			###
###		nBC = number of blue cones			###
###################################################################

connections <- function(bcPos,BCBPPos,nBC,nBCBP,maxdendr)   {

# Distance from BCBP to BCs
allPos <- rbind(bcPos,BCBPPos)
all.dist.matrix <- as.matrix(dist(allPos))
BCBPtoBCdist <- all.dist.matrix[1:(nBC),(nBC+1):(nBC+nBCBP)]

# Defining characteristic length as mean min distance between BCBP and BC
#r0 <- mean(apply(BCBPtoBCdist,2,min))

# Assuming exponential decay of connection probability, this essentially is the same as nearest neighbour.
exp.matrix <- exp(-BCBPtoBCdist)

# Calculate total probability
sums <- apply(exp.matrix,2,sum)

# Normalised probabilites
rel.p <- apply(exp.matrix,1, function(x) {x/sums})

# Making dendritic connection based on the probabilities (Gillespie algorithm)
reac_store <- vector()
for (i in 1:nBCBP) {
	rates <- rel.p[i,]
	lam <- sum(rates)
	rates <- rates / lam
	reac <- 1+sum(runif(1)>cumsum(rates))
	reac_store[i] <- reac
}

# BCBP more than 50 mu meter from BC don't make connections, ie. 44 mu meter in horizontal plane (ass. vertical dist of 23.2 mu meter), Kouyama 1992.
nonconnectingBCBP <- which(apply(BCBPtoBCdist,2,min) > maxdendr) # 44
reac_store[nonconnectingBCBP] <- NA

# Counting how many BCBP each BC connects to
nBCtoBCBP <- vector()
for (i in 1:nBC) {
	nBCtoBCBP[i] <- sum(reac_store == i, na.rm=TRUE )
}

# Removing dendrites to BCs which have more than four connections
nDendrRem <- nBCtoBCBP - 4
nDendrRem[which(nDendrRem < 0)] <- 0	 # Number of dendrites to remove from the various BCs

if(sum(nDendrRem)!=0) {

	mat <- matrix(ncol=2)
	for (i in 1:max(nDendrRem)) {
		mat <- rbind(mat,matrix(c(which(nDendrRem == i),rep(i,length(which(nDendrRem == i)))), ncol=2, nrow=length(which(nDendrRem == i))))
	}

	mat <- mat[-1,]

	if(is.matrix(mat)==FALSE) {mat <- matrix(mat, nrow=1, ncol=2)}	# If mat only has one row, we have to force back into matrix format.

	remove <- vector()
	for (i in 1:nrow(mat)) {
	remove <- c(remove,sample(which(reac_store == mat[i,1]),mat[i,2]))
	}

	reac_store[remove] <- NA
}

return(reac_store)
}








###################################################################
### migrations							### 
###-------------------------------------------------------------###
### Calculates post migrational positions of BCBPs (and BCs)	###
###-------------------------------------------------------------###
### Args:							###
###################################################################

migration <- function(nBCBP,vert.d.mean,vert.d.sd,vert.d.thr,mean.rel.force,sd.rel.force,dexcBCBP2,dendr.length,links,BCBPPos) {

# Assigning vertical distances between BCBP perikaryons and BC pedicles
vert.d <- rnorm(nBCBP,vert.d.mean,vert.d.sd)

# Selecting BCBPs which can move.
mobileBCBP <- which(vert.d < vert.d.thr)

# Assigning a friction force and a string force for each BC-BCBP connection
rel.force <- rnorm(nBCBP,mean.rel.force,sd.rel.force)

# Calculating new dendritic length
dendr.length.new <- (2*rel.force*sqrt(dendr.length^2+vert.d^2)-(rel.force^2+1)*dendr.length)/(1-rel.force^2)
dendr.length.new[which(rel.force > dendr.length/(sqrt(dendr.length^2+vert.d^2)))] <- dendr.length[which(rel.force > dendr.length/(sqrt(dendr.length^2+vert.d^2)))]
dendr.length.new[-mobileBCBP] <- dendr.length[-mobileBCBP]

# Finding new positions for BCBP
vector.diff <- cbind(links[,3]-links[,1],links[,4]-links[,2])
unit.vectors <- apply(vector.diff, 2, function(x) {x/dendr.length})
unit.vectors[which(is.na(unit.vectors))] <- 0	# Setting unit vectors of nonmoving BCBP to 0
change <- apply(unit.vectors, 2, function(x) {x*(dendr.length-dendr.length.new)})
change[which(is.na(change))] <- 0
BCBPPos.new <- BCBPPos + change
BCBPPos.new.old <- BCBPPos.new	# Saving the BCBP positions before fine tuning.

# Adjusting positions to make sure that BCBPs are not too close.
too.close <- 1

while (length(too.close) > 0) {

too.close <- which(apply(as.matrix(dist(BCBPPos.new)) < dexcBCBP2,1,sum) > 1)								# mobileBCBPs have to have a certain distance to all other BCBPs
too.close <- intersect(too.close, mobileBCBP)											# Only select BCBP among those in the moving layer.
too.close <- intersect(too.close, which(dendr.length.new!=dendr.length))							# Only select BCBP which have moved

# Looking for pairs in "too.close", only moving the member of a pair, which is the furthest away from the BC.
if (length(too.close)>1) {

temp <- as.matrix(dist(BCBPPos.new[too.close,]))
diag(temp) <- NA
pairs <- which(temp<8,arr.ind=T)

move.not <- vector()
if (nrow(pairs)!=0) {
	for (i in 1:nrow(pairs)  ) {
		if (dendr.length[too.close[pairs[i,1]]] > dendr.length[too.close[pairs[i,2]]]) {move.not <- c(move.not,too.close[pairs[i,2]])} else {move.not <- c(move.not,too.close[pairs[i,1]])}
	}
	move.not <- unique(move.not)
	}
else
	{move.not <- NULL}


too.close <- setdiff(too.close,move.not)

}	

BCBPPos.new[too.close,] <- BCBPPos.new[too.close,] - 0.5*unit.vectors[too.close,]		# Towards the end of movement, the speed is almost identical, hence same change for all.


}	# while loop ends

# Final positions of all blue cones and bipolars
newPos <- rbind(bcPos,BCBPPos.new)

return(newPos)

}
