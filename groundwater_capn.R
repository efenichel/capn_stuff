######################################################################################
# {capn} data: Kansas groundwater (gW.data) data
# Date: 12/30/2017
# Reference: Fenichel et al. (2016) PNAS
######################################################################################

#load required packages
if (!require("pacman")) install.packages("pacman") #pacman allows you to use the p_load function
#the p_load function checks is a library is installed, if not it installs it, then it attaches the 
#called library
p_load(capn, R.oo, repmis, ggplot2, devtools, RCurl)

#capn documentation: https://cran.r-project.org/web/packages/capn/capn.pdf 
#repmis documentation: https://cran.r-project.org/web/packages/repmis/index.html 

rm(list=ls()) #clear workspace
## troubleshooting branch
gitbranch <- "bygmd"
#get data set for the problem set from Github
source_data("https://github.com/efenichel/capn_stuff/raw/master/my_gw_data.RData")

ksdata <-readRDS(gzcon(url(paste0("https://github.com/efenichel/capn_stuff/raw/",gitbranch,"/KSwater_data.RDS")))) #RDS upload
#ksdata <- source_data("https://github.com/eaddicott/capn_stuff/raw/master/KSwater_data.RData") #Rdata file upload
#str(ksdata) # this line will show you the structure of the data

## STRUCTURE INFO
#The object will be a list of 7 lists. Each of the 7 lists corresponds to a groundwater management district (1-5), 
# the outgroup(6), or the entire state of Kansas(7).
# Each list will have 11 named elements:
# [1] $gmdnum: int   1:7
# [2] $mlogitcoeff: df with dim = (#cropcodes x 23) containing the coefficients and intercept terms from the multinomial logit model
# [3] $mlogitmeans: df containing the mean for each of the variables run in the mlogit model
# [4] $cropamts: df containing the mean acres planted to each of the 5 crop cover types for each cropcode
# [5] $watercoeff: df containing the water withdrawal regression coefficients and intercept
# [6] $wwdmeans: df containing the means for each of the variables in the water withdrawal regression
# [7] $costcropacre: df dim=(1x5) containing the cost per acre of planting each of the 5 crops
# [8] $cropprices: df dim=(5x1) containing the per unit prices of each of the crops
# [9] $meanwater: num mean AF water in the region of interest
# [10]$recharge: num recharge rate for the region of interest
# [11]$watermax: num upper bound of domain for node space, max water observed in region.

#load datasetup function
#source a separate script from github that contains the functions for gw system
if (!exists("datasetup", mode = "function")) { #only source if not exist
source(paste0("https://github.com/efenichel/capn_stuff/raw/",gitbranch,"/data_setup.R"))
}
#The elements of gw.data are
#The parameters from the multinomial logit for crop shares, 
#to know which is which see the labels on the additional csp.means data using the View(csp.means) command.
#The average amount planted in each field type defined by a wellhead 
#A water withdraw constant, wwd.alpha, which is element 3
#A water withdraw coefficient, wwd.beta, which is element 4
#Water withdraw crop coefficients, wwd.gamma, which is element 5
#Water withdraw crop coefficients main terms only, wwd.gamma.p1, which is element 6
#Water withdraw crop coefficients squared terms only, wwd.gamma.p2, which is element 7
#Per acre prices (gross revenues) used in the profit function, element 8
#Costs used in the profit function, element 9
#Also included are the raw coefficients for the crop shares models (csp.raw) and their mean values csp.means.  
##################################################################################################
##################################################################################################
##
##  Modeling parameters to control   #############################################################
##
##################################################################################################
###
region <- 6 # Select region,  1:5 are GMD, 6 is outgroup, 7 is state
###
# After setting the region, create capn data structure

region_data <- ksdata[[region]] #double-brackets here important. Load in region specific data
gw.data <- datasetup(region) # if running line by line, be sure to run this function at the bottom of this script
#Economic parameters
dr <- 0.03 #discount rate

#System parameters
recharge <- region_data[[10]] #units are inches per year constant rate 
#recharge <- 1.25 uncomment to input your own for sensitivity analysis

#capN parameters
order <- 10 # approximaton order
NumNodes <- 100 #number of nodes
wmax <- region_data[[11]]

#####################################################################
# Get the system model associated with Fenichel et al. 2016
# First load the functions 
#source('system_fns.R')  
#source a separate script from github that contains the functions for gw system
if (!exists("cropFwater", mode = "function")) { #only sources if doesn't exist
source(paste0("https://raw.githubusercontent.com/efenichel/capn_stuff/",gitbranch,"/system_fns.R")) #will need to fix in the end
}
#notice in the explore tab there are now a bunch of Functions.
#you can click on any function and see the functions that support the system model.
#These are the functions that you need to build and calibrate for your system.

######################################################################
######################################################################
######################################################################

#     TEST capn
######################################################################

#Prepare {capn}

#prepare capN
Aspace <- aproxdef(order,0,wmax,dr) #defines the approximation space
nodes <- chebnodegen(NumNodes,0,wmax) #define the nodes

#prepare for simulation
simuData <- matrix(0,nrow = NumNodes, ncol = 5)

#simulate at nodes
for(j in 1:NumNodes){
  simuData[j,1]<-nodes[j] #water depth nodes
  simuData[j,2]<-sdot(nodes[j],recharge,gw.data) #change in stock over change in time
  simuData[j,3]<- 0-WwdDs1(nodes[j],gw.data) # d(sdot)/ds, of the rate of change in the change of stock
  simuData[j,4]<-ProfDs1(nodes[j],gw.data) #Change in profit with the change in stock
  simuData[j,5]<-profit(nodes[j],gw.data) #profit
}


#recover approximating coefficents
pC<-paprox(Aspace, simuData[,1],simuData[,2],simuData[,3],simuData[,4])  #the approximated coefficent vector for prices

#project shadow prices, value function, and inclusive wealth
waterSim <- psim(pcoeff = pC,
                 stock = simuData[,1],
                 wval = simuData[,5],
                 sdot = simuData[,2])
#convert to data frame 
waterSim <- as.data.frame(waterSim)

#use ggplot to plot the water shadow price function  
ggplot() + 
  geom_line(data = waterSim, aes(x = stock, y = shadowp),
            color = 'blue') +
            labs( 
            x= "Stored groundwater",
            y = "Shadow price")  +
            theme(  #http://ggplot2.tidyverse.org/reference/theme.html
              axis.line = element_line(color = "black"), 
              panel.background = element_rect(fill = "transparent",colour = NA),
              plot.background = element_rect(fill = "transparent",colour = NA)
            )
cat("if everything runs well the next line should say 16.82325", "\n")
cat("At 21.5 acre feet of water, the shadow price is" , psim(pC,21.5)$shadowp, "\n")

#use ggplot plot the value function
ggplot() + 
  geom_line(data = waterSim[5:100,], aes(x = stock[5:100], y = vfun[5:100]),
            # four rows are removed because they are two far outside the data
            color = 'blue') +
  labs( 
    x= "Stored groundwater",
    y = "Intertemporial Welfare")  +
  theme(  #http://ggplot2.tidyverse.org/reference/theme.html
    axis.line = element_line(color = "black"), 
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )

cat("you may get an warning that some number or rows were removed from the plot. This is not important.", "\n")

testme<-psim(pcoeff = pC, 
             stock = c(18.5,21.5), 
             wval = c(profit(18.5,gw.data),profit(21.5, gw.data)),
             sdot = c(sdot(18.5, recharge, gw.data),sdot(21.5, recharge, gw.data))
             )
testme

rm(j, testme)

