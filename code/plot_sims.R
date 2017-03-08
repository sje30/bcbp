# This file should contain 2 functions:
# 1) a function to plot the cross pvalues when there is one parameter
#    whose values are being looped over
# 2) a function to do a levelplot when there are 2 parameters whose values
#    are being looped over

require(lattice)

########################################################################
### plot_lines					                                     ### 
###------------------------------------------------------------------###
### Plots the line graph for cross pvalues                           ###
###------------------------------------------------------------------###
### Args:    pvalues = data frame of pvalues                         ###
###    commandlineargs = commandline arguments to                    ###
###                       'sensitivity_analysis' without             ###
###                       text based args                            ###
###		param = parameter to be tested		          	             ###
########################################################################
plot_lines <- function(pvalues, commandlineargs, param) {
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
	
	colnames(pvalues) <- c("pBC","pBCBP","pBCBCBP")
	rownames(pvalues) <- range1
	
	pBCBCBP <- pvalues$pBCBCBP

	##### Make plot of pBCBCBP (cross value) as function of param
	pdf(file=paste("plot_", param, "_field",x,"_",low,"_",high,"_",incr,"_",loopsize,"_",ulimit,".pdf",sep=""))

	# set labels according to parameter tested
	if(param == "depth") {
		xlabel <- "Max depth of mobile BCBPs"
		title <- "depth"
	} else if(param == "cone_excl") {
		xlabel <- "Cone exclusion zone"
		title <- paste("Sensitivity to cone exclusion zone, field ",x,sep="")
	} else if(param == "bc_excl_mean") {
		xlabel <- "Blue Cone exclusion zone mean"
		title <- paste("Sensitivity to blue cone exclusion zone mean, field ",x,sep="")
	} else if(param == "bc_excl_sd") {
		xlabel <- "Blue Cone exclusion zone sd"
		title <- paste("Sensitivity to blue cone exclusion zone sd, field ",x,sep="")
	} else if(param == "bc_excl_trunc") {
		xlabel <- "Blue Cone exclusion zone truncation"
		title <- paste("Sensitivity to blue cone exclusion zone truncation, field ",x,sep="")
	} else if(param == "bcbp_excl_mean") {
		xlabel <- "Blue Cone Bipoler exclusion zone mean"
		title <- paste("Sensitivity to blue cone bipoler exclusion zone mean, field ",x,sep="")
	} else if(param == "bcbp_excl_sd") {
		xlabel <- "Blue Cone Bipoler exclusion zone sd"
		title <- paste("Sensitivity to blue cone bipoler exclusion zone sd, field ",x,sep="")
	} else if(param == "bcbp_excl_trunc") {
		xlabel <- "Blue Cone Bipoler exclusion zone truncation"
		title <- paste("Sensitivity to blue cone bipoler exclusion zone truncation, field ",x,sep="")
	} else if(param == "max_dendr_length") {
		xlabel <- "Maximum dendrite length"
		title <- paste("Sensitivity to maximum dendrite length, field ",x,sep="")
	} else if(param == "rel_force_mean") {
		xlabel <- "Mean relative force"
		title <- paste("Sensitivity to mean relative force, field ",x,sep="")
	}

	plot(range1, pBCBCBP, type="b", lwd=1, xlab=xlabel, ylab="Cross p value", main=title)

	text((low+high)/2,(max(pBCBCBP)+min(pBCBCBP))/2,paste("Field = ",x,"\n","Lower limit =",low,"\n","Upper limit =",high,"\n","Step size =",incr,"\n","Number of iterations =",loopsize,"\n","Upper integation limit =",ulimit))

	dev.off()
}

########################################################################
### plot_level					                                     ### 
###------------------------------------------------------------------###
### Plots a level plot for cross pvalues with 2 parameters varying   ###
###------------------------------------------------------------------###
### Args:    pvalues = data frame of pvalues                         ###
###    commandlineargs = commandline arguments to                    ###
###                       'sensitivity_analysis' without             ###
###                       text based args                            ###
###		param1 = first parameter to be tested		                 ###
###     param2 = second parameter to be tested                       ###
########################################################################
plot_level <- function(pvalues, commandlineargs, param1, param2) {
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

	################################
	### Make level plot
	################################
	p1 <- pvalues$param1_value
	p2 <- pvalues$param2_value
	g1 <- cbind(p1, p2)
	g1$pval <- pvalues$pBCBCBP
	fieldnum <- paste("field", x)
	
	g1$field <- rep(fieldnum, length(p1))

	df <- as.data.frame(g1)
	
	##### Make level plot of pBCBCBP (cross value) as function of param1 and param2
	pdf(file=paste("plot_",param1,"_",low1,"_",high1,"_",incr1,"_",param2,"_",low2,"_",
	               high2,"_",incr2,"_field",x,"_",loopsize,"_",ulimit,".pdf",sep=""))

	# print because levelplot alone leaves the pdf blank when used in conjunction with dev.off()
	print(levelplot(pval~p1*p2|field, data=df, xlab=param1, ylab=param2,
	          cuts=50, main="P values", colorkey=TRUE, region=TRUE))
	dev.off()
}
	




















