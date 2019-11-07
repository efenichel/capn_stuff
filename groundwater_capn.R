######################################################################################
# {capn} data: Kansas groundwater (gW.data) data
# Date: 12/30/2017
# Reference: Fenichel et al. (2016) PNAS
######################################################################################

#load required packages
if (!require("pacman")) install.packages("pacman") #pacman allows you to use the p_load function
#the p_load function checks is a library is installed, if not it installs it, then it attaches the 
#called library
p_load(capn, R.oo, repmis, ggplot2, RCurl)

#capn documentation: https://cran.r-project.org/web/packages/capn/capn.pdf 
#repmis documentation: https://cran.r-project.org/web/packages/repmis/index.html 

rm(list=ls()) #clear workspace

#get data set for the problem set from Github
source_data("https://github.com/efenichel/capn_stuff/raw/master/my_gw_data.RData")


ksdata <-readRDS(gzcon(url("https://github.com/efenichel/capn_stuff/raw/forspi/KSwater_data.RDS")))

datasetup(2)

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
#additional modeling parameters
dr <- 0.03 #discount rate
recharge <- 1.25 #inches per year constant rate

#capN parameters
order <- 10 # approximaton order
NumNodes <- 100 #number of nodes

#####################################################################
# Get the system model associated with Fenichel et al. 2016
# First load the functions 
#source('system_fns.R')  

source("https://raw.githubusercontent.com/efenichel/capn_stuff/forspi/system_fns.R") #will need to fix in the end

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
Aspace <- aproxdef(order,0,116,dr) #defines the approximation space
nodes <- chebnodegen(NumNodes,0,116) #define the nodes

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

