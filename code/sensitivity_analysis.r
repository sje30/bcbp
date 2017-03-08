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

simulation.func <- sprintf("%s/simulate.R", AVLROOT)
source(simulation.func)

process.funcs <- sprintf("%s/process_results.R", AVLROOT)
source(process.funcs)

plot.funcs <- sprintf("%s/plot_sims.R", AVLROOT)
source(plot.funcs)

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
#### parameter to be tested for sensitivity, minimum (of that parameter), maximum, increment,
#### 2nd parameter (optional),minimum2, maximum2, increment2

# get commandline arguments as strings
commandlineargs <- commandArgs(TRUE)

###########################
## Run simulations
###########################
param1 <- commandlineargs[4]       # which parameter is to be tested for sensitivity (this will be a string)
# can be: 'depth', 'cone_excl', 'bc_excl_mean', 'bc_excl_sd', 'bc_excl_trunc', 'bcbp_excl_mean',
#         'bcbp_excl_sd', 'bcbp_excl_trunc', 'max_dendr_length', 'rel_force_mean'

# fieldnumber
x <- commandlineargs[1]

# toggle this depending on whether you want the L-function and DRP plots
plot_l_drp = TRUE

if(!is.na(commandlineargs[8])) {
    param2 <- commandlineargs[8]    # second parameter to be tested
    # convert all the numeric arguments to numeric
    commandlineargs <- as.numeric(commandlineargs[-c(4, 8)])
    low2 <- commandlineargs[7]    # minimum value of 'param2' to be tested
    print(low2)
    high2 <- commandlineargs[8]   # maximum value of 'param2' to be tested
    print(high2)
    incr2 <- commandlineargs[9]  # increment for 'param2'
    print(incr2)
    param2.values <- seq(low2, high2, incr2)
    n <- length(param2.values)
    sim_out <- list()
    for(i in 1:n) {
        sim_out[[i]] <- simulate(commandlineargs, param1, param2, param2.values[i], AVLROOT)
    }
    process2params(sim_out, commandlineargs, param1, param2, plot_l_drp)
} else {
    commandlineargs <- as.numeric(commandlineargs[-4])
    sim_out <- simulate(commandlineargs, param1, NA, NA, AVLROOT)
    process1param(sim_out, commandlineagrs, param1, plot_l_drp)
}

# save R image
#save.image(file=paste("rdata_", param, "_",x,"_",low,"_",high,"_",incr,"_",loopsize,"_",ulimit,"_",dexcBC,"_",dexcBC.sd,"_",dexcBC.trunc,"_",dexcBCBP,"_",dexcBCBP.sd,"_",dexcBCBP.trunc,".RData",sep=""))

        











