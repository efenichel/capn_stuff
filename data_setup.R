# -----------------------------------------------------------------------------
# File Name: datasetup function
# Purpose: Process data from regressions for the groundwater problem
# Author:  Ethan Addicott
# Date:  11/06/19
# Notes: 

#------------------------------------------------------------------------------

# Convert regression output and summary statistics into capn parameters -------------------------------------
datasetup <- function(gmdnum, dataset = ksdata) {
  region_data <- dataset[[gmdnum]]
  # Region Specific Parameters ---------------------------------------------------
  #CASE 1: GMD 1-5
  if (gmdnum < 6) {
    # Crop Choice Parameters Setup ------------------------------------------------
    `mlogit.gmd.coeff` <- region_data[[2]]
    `mlogit.gmd.means` <- region_data[[3]]
    
    mlogit.gmd.coeff   <- data.frame(mlogit.gmd.coeff[-c(1,24,25),], 
                                     stringsAsFactors = FALSE)
    mlogit.gmd.means   <- data.frame(mlogit.gmd.means[-c(22,23), gmdnum + 1], 
                                     stringsAsFactors = FALSE)
    
    # Transpose betas, chop off variable names to create crop.coeff column
    crop.gmd.betas     <- t(mlogit.gmd.coeff[1,-c(1,ncol(mlogit.gmd.coeff))])
    crop.gmd.alphas    <- rep(NA, ncol(mlogit.gmd.coeff) - 2)
    crop.mean.vector <- append(as.numeric(unlist(mlogit.gmd.means[-c(1),])), 1)
    max.j <- ncol(mlogit.gmd.coeff) - 1
    for (j in 2:max.j) {
      crop.gmd.alphas[j - 1] <- as.numeric(mlogit.gmd.coeff[-c(1), j]) %*% 
        as.numeric(crop.mean.vector)
    }
    crop.coeff.gmd <- data.frame(crop.gmd.alphas,crop.gmd.betas, 
                                 stringsAsFactors = FALSE)
    colnames(crop.coeff.gmd) <- c("alpha","beta")
    
    # Crop Amounts Parameter Setup ------------------------------------------------
    `crop.amts.gmd` <- region_data[[4]]
    crop.amts.gmd <- data.matrix(crop.amts.gmd[-c(1,8),-c(1)])
    crop.amts.gmd <- t(crop.amts.gmd)
    crop.amts.gmd <- transform(crop.amts.gmd, numeric)
    #Check that crop amounts are one more than 
    if (nrow(crop.amts.gmd) != max.j) {
      stop("BAD LENGTH")
    }
    # Water Withdrawal Parameter Setup --------------------------------------------
    `water.coeff` <- region_data[[5]]
    rmse <- as.numeric(water.coeff[36,2])
    water.coeff <- data.frame(water.coeff[-c(1,34:37),] ,stringsAsFactors = FALSE)
    
    # Beta ------------------------------------------------------------------------
    beta <- as.numeric(water.coeff[21,2])
    
    # Gamma -----------------------------------------------------------------------
    gamma <- as.numeric(water.coeff[c(1:10),2])
    gamma1 <- gamma[c(1,3,5,7,9)]
    gamma2 <- gamma[c(2,4,6,8,10)]
    
    # Alpha -----------------------------------------------------------------------
    alpha.water.coeff <- water.coeff[-c(1:11,21),2]
    `wwd.means` <- region_data[[6]]
    wwd.means <- data.frame(wwd.means[-c(1,11,22,23),gmdnum + 2], stringsAsFactors = FALSE)
    wwd.mean.vector <- append(as.numeric(unlist(wwd.means)),1)
    alpha <- as.numeric(alpha.water.coeff) %*% as.numeric(wwd.mean.vector)
    alpha <- alpha + (rmse**2)/2 #correction for log estimation
    alpha <- as.vector(alpha, mode = "numeric")
  }
  #CASE 2: GMD 6 Outgroup
  else if (gmdnum == 6){
    # Crop Choice Parameters Setup ------------------------------------------------
    `mlogit.state.coeff` <- dataset[[7]][[2]]
    `mlogit.state.means` <- dataset[[7]][[3]]
    mlogit.state.coeff   <- data.frame(mlogit.state.coeff[-c(1,24,25),], 
                                       stringsAsFactors = FALSE)
    mlogit.state.means   <- data.frame(mlogit.state.means[-c(22,23),], 
                                       stringsAsFactors = FALSE)
    # Transpose betas, chop off variable names to create crop.coeff column
    crop.state.betas     <- t(mlogit.state.coeff[1,-c(1,ncol(mlogit.state.coeff))])
    crop.state.alphas    <- rep(NA, ncol(mlogit.state.coeff) - 2)
    crop.mean.vector <- append(as.numeric(mlogit.state.means$mean[-c(1)]),1)
    max.j <- ncol(mlogit.state.coeff) - 1
    for (j in 2:max.j) {
      crop.state.alphas[j - 1] <- as.numeric(mlogit.state.coeff[-c(1), j]) %*% 
        as.numeric(crop.mean.vector)
    }
    crop.coeff.state <- data.frame(crop.state.alphas,crop.state.betas, 
                                   stringsAsFactors = FALSE)
    colnames(crop.coeff.state) <- c("alpha","beta")
    
    # Crop Amounts Parameter Setup ------------------------------------------------
    `crop.amts.state` <- dataset[[7]][[4]]
    crop.amts.state <- data.matrix(crop.amts.state[-c(1,8),-c(1)])
    crop.amts.state <- t(crop.amts.state)
    
    #Check that crop amounts are one more than 
    if (nrow(crop.amts.state) != max.j) {
      stop("BAD LENGTH")
    }
    
    # Water Withdrawal Parameter Setup --------------------------------------------
    `water.coeff` <- region_data[[5]]
    rmse <- as.numeric(water.coeff[36,2])
    water.coeff <- data.frame(water.coeff[-c(1,34:37),] ,stringsAsFactors = FALSE)
    
    # Beta ------------------------------------------------------------------------
    beta <- as.numeric(water.coeff[21,2])
    
    # Gamma -----------------------------------------------------------------------
    gamma <- as.numeric(water.coeff[c(1:10),2])
    gamma1 <- gamma[c(1,3,5,7,9)]
    gamma2 <- gamma[c(2,4,6,8,10)]
    
    # Alpha -----------------------------------------------------------------------
    alpha.water.coeff <- water.coeff[-c(1:11,21),2]
    `wwd.means` <- region_data[[6]]
    wwd.means <- data.frame(wwd.means[-c(1,11,22,23),gmdnum + 2], stringsAsFactors = FALSE)
    wwd.mean.vector <- append(as.numeric(unlist(wwd.means)),1)
    alpha <- as.numeric(alpha.water.coeff) %*% as.numeric(wwd.mean.vector)
    alpha <- alpha + (rmse**2)/2 #correction for log estimation
    alpha <- as.vector(alpha, mode = "numeric")
  }
  #CASE 3: STATE
  else if (gmdnum == 7){
    # Crop Choice Parameters Setup ------------------------------------------------
    `mlogit.state.coeff` <- region_data[[2]]
    `mlogit.state.means` <- region_data[[3]]
    mlogit.state.coeff   <- data.frame(mlogit.state.coeff[-c(1,24,25),], 
                                       stringsAsFactors = FALSE)
    mlogit.state.means   <- data.frame(mlogit.state.means[-c(22,23),], 
                                       stringsAsFactors = FALSE)
    # Transpose betas, chop off variable names to create crop.coeff column
    crop.state.betas     <- t(mlogit.state.coeff[1,-c(1,ncol(mlogit.state.coeff))])
    crop.state.alphas    <- rep(NA, ncol(mlogit.state.coeff) - 2)
    crop.mean.vector <- append(as.numeric(mlogit.state.means$mean[-c(1)]),1)
    max.j <- ncol(mlogit.state.coeff) - 1
    for (j in 2:max.j) {
      crop.state.alphas[j - 1] <- as.numeric(mlogit.state.coeff[-c(1), j]) %*% 
        as.numeric(crop.mean.vector)
    }
    crop.coeff.state <- data.frame(crop.state.alphas,crop.state.betas, 
                                   stringsAsFactors = FALSE)
    colnames(crop.coeff.state) <- c("alpha","beta")
    
    # Crop Amounts Parameter Setup ------------------------------------------------
    `crop.amts.state` <- region_data[[4]]
    crop.amts.state <- data.matrix(crop.amts.state[-c(1,8),-c(1)])
    crop.amts.state <- t(crop.amts.state)
    
    #Check that crop amounts are one more than 
    if (nrow(crop.amts.state) != max.j) {
      stop("BAD LENGTH")
    }
    
    # Water Withdrawal Parameter Setup --------------------------------------------
    `water.coeff` <- region_data[[5]]
    rmse <- as.numeric(water.coeff[36,2])
    water.coeff <- data.frame(water.coeff[-c(1,34:37),] ,stringsAsFactors = FALSE)
    
    # Beta ------------------------------------------------------------------------
    beta <- as.numeric(water.coeff[21,2])
    
    # Gamma -----------------------------------------------------------------------
    gamma <- as.numeric(water.coeff[c(1:10),2])
    gamma1 <- gamma[c(1,3,5,7,9)]
    gamma2 <- gamma[c(2,4,6,8,10)]
    
    # Alpha -----------------------------------------------------------------------
    alpha.water.coeff <- water.coeff[-c(1:11,21),2]
    `wwd.means` <- region_data[[6]]
    wwd.means <- data.frame(wwd.means[-c(1,11,22,23),], stringsAsFactors = FALSE)
    wwd.mean.vector <- append(as.numeric(wwd.means$all),1)
    alpha <- as.numeric(alpha.water.coeff) %*% as.numeric(wwd.mean.vector)
    alpha <- alpha + (rmse**2)/2 #correction for log estimation
    alpha <- as.vector(alpha, mode = "numeric")
  }
  
  # CONSISTENT PARAMETERS ------------------------------------------------------
  # Cost Crop Acre 
  `cost.crop.acre` <- region_data[[7]]
  cost.crop.acre <- as.vector(unlist(cost.crop.acre), mode = "numeric")
  
  # Crop Prices 
  `crop.prices` <- region_data[[8]]
  crop.prices <- as.vector(unlist(crop.prices), mode = "numeric")
  
  
  # Store Parameters in Data Structure ------------------------------------------
  gw.data <- rep(NA,9) #preallocate
  # GMD 1-5
  if (gmdnum < 6) {
    gw.data <- list(crop.coeff.gmd, crop.amts.gmd, alpha, beta, gamma, 
                    gamma1, gamma2, crop.prices, cost.crop.acre)
  }
  # GMD 6 (Outgroup) or State
  if (gmdnum >= 6) {
    gw.data <- list(crop.coeff.state, crop.amts.state, alpha, beta, gamma, 
                    gamma1, gamma2, crop.prices, cost.crop.acre)
  }
  
  return(gw.data)
} 
