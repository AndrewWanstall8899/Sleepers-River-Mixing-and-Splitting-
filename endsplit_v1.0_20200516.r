################################################################################
#
# EndSplit -- Script for end-member mixing and end-member splitting
#
# Version 1.0  build 2020.05.16
# Author: James Kirchner, ETH Zurich
#
#
#
# Copyright (C) 2020 ETH Zurich and James Kirchner
# Public use of this script is permitted under GNU General Public License 3 (GPL3); for details see <https://www.gnu.org/licenses/>
#
# READ THIS CAREFULLY:
# ETH Zurich and James Kirchner make ABSOLUTELY NO WARRANTIES OF ANY KIND, including NO WARRANTIES, expressed or implied, that this software is
#    free of errors or is suitable for any particular purpose.  Users are solely responsible for determining the suitability and
#    reliability of this software for their own purposes.
#
#
# ALSO READ THIS:
# This script implements end-member mixing and end-member splitting, as described in Kirchner and Allen, "Seasonal partitioning of 
# precipitation between streamflow and evapotranspiration, inferred from end-member splitting analysis", Hydrology and Earth System Sciences.
# Users publishing results based on this script should cite that paper.
#
#
# The equation numbers here refer to the corresponding equations in that paper (and its supplement, for equation numbers beginning with "s").
# The equations here may differ from the text in certain details due to (for example) the array indexing conventions of the R language.
#
# For an illustration of how to apply this script, see the accompanying demonstration script EndSplit_demo_v1.0_20191024
#



##############################################
# calculate weighted mean and standard error
##############################################
wtd_mean <- function(x, wt=rep(1, length(x))) {
  # Calculates weighted mean and its standard error, using the correct standard error formula
  # for cases where weights are measures of importance or volume.  Note that this is different
  # from the common weighted standard error formula, where the weights are assumed to be 
  # inversely proportional to the variance of each input value.  For highly uneven weights,
  # the conventional weighted standard error formula will converge toward zero (which is wrong for our problem),
  # whereas the correct formula will converge toward infinity (which is correct).

  # receives as input:
  # x is a numeric vector of values to be averaged
  # wt is a numeric vector of weights.  Weights must be >=0
  # x and wt must be of equal length
  
  # returns a list with two elements:
  # mean is the weighted mean
  # SE is the standard error

  # first we null out all values of x that don't have weights, and vice versa.  
  x[is.na(wt)] <- NA
  wt[is.na(x)] <- NA
  
  #check for different lengths of vectors
  if (length(x) != length(wt)) stop("error in wtd_mean: x and wt have different length")
  
  #check for negative weights
  if (sum((wt<0), na.rm=TRUE)>0) stop("error in wtd_mean: negative weights")
  
  #sum of weights
  sumwt <- sum(wt, na.rm=TRUE)
  
  #sum of squared weights
  sumsq <- sum(wt*wt, na.rm=TRUE)
  
  #effective degrees of freedom (note that this equals n for even weights, and trends toward 1 for very uneven weights)
  n_eff <- sumwt*sumwt/sumsq                                              #Equation (S6)
    
  #weighted mean
  xbar <- sum(x*wt, na.rm=TRUE)/sumwt                                     #Equation (S4)
  
  #weighted variance
  varx <- (sum(wt*((x-xbar)^2), na.rm=TRUE)/sumwt) * n_eff/(n_eff-1)      #Equation (S8)
  
  return(
    list(mean = xbar,         #weighted mean
         se = sqrt(varx/n_eff)  #standard error of weighted mean          Equations (S8) and (S9)
         )
  )
 
}








##############################################
# perform end-member splitting calculations
##############################################
EndSplit <- function(Pdel, 
                     Qdel, 
                     Pwt=rep(1, length(Pdel)), 
                     Qwt=rep(1, length(Qdel)), 
                     Pdelcat, 
                     Qdelcat=rep(1, length(Qdel)), 
                     P, Q, 
                     Pcat, Qcat=rep(1, length(Qdel)), n.yr=1){
  
  # receives as input:
  #
  #    Pdel    - vector of precipitation isotope del values (or tracer concentrations)
  #    Pwt     - vector of weights for precipitation dels (usually total precip volume over sampling interval)
  #    Pdelcat - vector of categories into which precipitation is divided (seasons, precipitation types, etc.)
  #                There must be exactly two such categories.  Categories can be labeled by numbers or alphanumerics,
  #                    since we will be handling them as factors
  #    Pdel, Pwt, and Pdelcat must all be the same length (this is checked).
  #
  #    Qdel    - vector of streamflow isotope del values (or tracer concentration)
  #    Qwt     - vector of weights for streamflow dels (usually total streamflow volume over sampling interval, or instantaneous discharge)
  #    Qdelcat - vector of categories into which streamflow is divided (seasons, precipitation types, etc.)
  #                There must be at least one category.  Categories can be labeled by numbers or alphanumerics,
  #                    since we will be handling them as factors
  #    Qdel, Qwt, and Qdelcat must all be the same length (this is checked).
  #
  #    P    - a vector of precipitation flux (in volume/area/time or depth/time)
  #    Pcat - categories into which precipitation is divided (seasons, precipitation types, etc.)
  #                There must be exactly two such categories.  Categories can be labeled by numbers or alphanumerics,
  #                    since we will be handling them as factors.  The labels used in Pcat must be IDENTICAL to 
  #                    the labels used in Pdelcat; for example, you can't use "summer" and "winter" in Pcat, and
  #                    "s" and "W" in Pdelcat.
  #
  #    Q    - streamflow flux (in volume/area/time or depth/time)
  #    Qcat - categories into which streamflow is divided (seasons, precipitation types, etc.)
  #                There must be at least one category.  Categories can be named by numbers or alphanumerics,
  #                    since we will be handling them as factors.  The labels used in Qcat must be IDENTICAL to 
  #                    the labels used in Qdelcat; for example, you can't use "growing" and "dormant" in Qcat, and
  #                    "grow" and "form" in Qdelcat.
  #
  #    P, Q, Pcat, and Qcat should be of the same length.  They should also be regularly spaced (e.g., hourly, daily, etc.) and should
  #    cover the same span of time (ideally, complete water years).  P and Q must be measured in the same units.
  #
  #    n.yr - number of water years covered by P and Q data (to convert total fluxes to yearly fluxes).  Specifying the 
  #              wrong value of n.yr will yield the wrong yearly fluxes, but the end-member mixing and splitting fractions
  #              will be unaffected.
  #
  
  
  # returns a list consisting of the following:
  #
  #    table  - a table with the following four data frames
  #        f      - data frame of ("backward") end-member mixing fractions
  #        f.se   - data frame of standard errors of end-member mixing fractions
  #        eta    - data frame of ("forward") end-member splitting fractions
  #        eta.se - data frame of standard errors of end-member splitting fractions
  #        These data frames each have two columns (corresponding to the two categories (or seasons) of inputs) and
  #            one row for each discharge category (or discharge season), plus two more rows.  The next-to-last row,
  #            "allQ", is the sum of all discharge categories, and the last row, "nonQ", is the sum of water leaving
  #            the catchment by all other pathways (in practice this will usually be evapotranspiration, but it could
  #            include any other unmeasured fluxes).
  #    Pdel.bar  - average precipitation tracer signature (by precipitation category or season)
  #    Pdel.se   - standard error of Pdel.bar
  #    Qdel.bar  - average streamflow tracer signature (by discharge category or season)
  #    Qdel.se   - standard error of Qdel.bar
  #    Ptot      - total precipitation flux (by precipitation category or season)
  #    Ptot.se   - standard error of Ptot
  #    Qtot      - total discharge flux (by discharge category or season)
  #    Qtot.se   - standard error of Qtot

  
  
  
  
  # coerce category vectors to factors
  if (!is.factor(Pcat)) Pcat <- factor(Pcat, ordered=TRUE)
  if (!is.factor(Pdelcat)) Pdelcat <- factor(Pdelcat, ordered=TRUE)
  if (!is.factor(Qcat)) Qcat <- factor(Qcat, ordered=TRUE)
  if (!is.factor(Qdelcat)) Qdelcat <- factor(Qdelcat, ordered=TRUE)
  
  # check that category labels match
  if (sum(levels(Pcat)!=levels(Pdelcat))>0) stop ("fatal error: Pcat and Pdelcat use different labels")
  if (sum(levels(Qcat)!=levels(Qdelcat))>0) stop ("fatal error: Qcat and Qdelcat use different labels")
  if (nlevels(Pcat)!=2) stop ("fatal error: need exactly two precipitation categories")

  # check that vector lengths match  
  n <- length(Pdel)
  if ((length(Pwt) != n) | (length(Pdelcat) != n)) stop ("fatal error: Pdel, Pwt, and Pdelcat must have the same length")
  n <- length(Qdel)
  if ((length(Qwt) != n) | (length(Qdelcat) != n)) stop ("fatal error: Qdel, Qwt, and Qdelcat must have the same length")
  if (length(P) != length(Pcat)) stop ("fatal error: P and Pcat must have the same length")
  if (length(Q) != length(Qcat)) stop ("fatal error: Q and Qcat must have the same length")
  

  # calculate means and standard errors

  # firstcreate vectors and fill them with NA  
  Pdel.bar <- array(NA, dim = nlevels(Pcat), dimnames = list(levels(Pcat)))                      #Vector of average precipitation del values
  Pdel.se <- array(NA, dim = nlevels(Pcat), dimnames = list(paste(levels(Pcat),".se", sep="")))  #Vector of precipitation del standard errors
  Ptot <- Pdel.bar;                                                                              #Vector of total precipitation fluxes
  Ptot.se <- Pdel.se;                                                                             #Vector of precipitation flux standard errors
  Qdel.bar <- array(NA, dim = nlevels(Qcat)+2, dimnames = list(c(levels(Qcat),"allQ", "nonQ")))  #Vector of average streamflow del values
  Qdel.se <- array(NA, dim = nlevels(Qcat)+2, dimnames = list(c(paste(levels(Qcat),".se",sep=""),"allQ.se", "nonQ.se"))) #Vector of streamflow del standard errors
  Qtot <- Qdel.bar;                                                                              #Vector of total streamflow fluxes
  Qtot.se <- Qdel.se;                                                                            #Vector of streamflow flux standard errors
  
  
  # we have introduced two extra discharge categories: "AllQ" (total discharge), and "nonQ" (whatever is not sampled as discharge)
  # "nonQ" will normally be interpreted as ET, but fundamentally it is all unmeasured fluxes
  
  for(i in 1:nlevels(Pcat)){ # loop through precipitation categories
    x <- wtd_mean(Pdel[as.integer(Pdelcat)==i], Pwt[as.integer(Pdelcat)==i]) #calculate weighted mean
    Pdel.bar[i] <- x$mean   # copy wtd mean to correct element of vector
    Pdel.se[i] <- x$se     # copy wtd std. error to correct element of vector
    
    x <- wtd_mean(P[as.integer(Pcat)==i]) # calculate unweighted mean (which works b/c we are not passing a weight vector)
    Ptot[i] <- x$mean*sum(as.integer(Pcat)==i)/n.yr   # copy wtd mean to correct element of vector, multiply by time steps in that category, and divide by number of years to get category average per year
    Ptot.se[i] <- x$se*sum(as.integer(Pcat)==i)/n.yr     # copy wtd std. error to correct element of vector, multiply by time steps in that category, and divide by number of years to get category average per year
  }
  
  for(j in 1:nlevels(Qcat)){ # loop through streamflow categories
    x <- wtd_mean(Qdel[as.integer(Qdelcat)==j], Qwt[as.integer(Qdelcat)==j]) #calculate weighted mean
    Qdel.bar[j] <- x$mean   # copy wtd mean to correct element of vector
    Qdel.se[j] <- x$se     # copy wtd std. error to correct element of vector
    
    x <- wtd_mean(Q[as.integer(Qcat)==j]) # calculate unweighted mean (which works b/c we are not passing a weight vector)
    Qtot[j] <- x$mean*sum(as.integer(Qcat)==j)/n.yr   # copy wtd mean to correct element of vector, multiply by time steps in that category, and divide by number of years to get category average per year
    Qtot.se[j] <- x$se*sum(as.integer(Qcat)==j)/n.yr     # copy wtd std. error to correct element of vector, multiply by time steps in that category, and divide by number of years to get category average per year
  }

  # Now we need to calculate the Qdel and Qtot for the total discharge, and for the ET (or "nonQ") by mass balance.
  # Note that we should not do this from the raw data, because if the sampling frequency varies between categories that will introduce a bias.
  
  # For clarity, we will do these steps outside of the arrays Qdel.bar (etc.) and then copy the values over. 
  
  AllQ <- sum(Qtot, na.rm=TRUE) # sum all Q
  AllQ.se <- sqrt(sum(Qtot.se*Qtot.se , na.rm=TRUE))  # Gaussian error propagation
  AllQ.del <- sum( Qdel.bar*Qtot, na.rm=TRUE )/ AllQ  # weighted average del value
  AllQ.del.se <- sqrt( sum((Qdel.se*Qtot/AllQ)^2, na.rm=TRUE) + sum((Qdel.bar*Qtot.se*(AllQ-Qtot)/AllQ^2)^2, na.rm=TRUE)) # Gaussian error propagation
                             # note that the second sum recognizes that each Qtot is part of AllQ so they are handled jointly
  
  AllP <- sum(Ptot, na.rm=TRUE) # sum of all P
  AllP.se <- sqrt(sum(Ptot.se*Ptot.se , na.rm=TRUE))  # Gaussian error propagation
  AllP.del <- sum( Pdel.bar*Ptot, na.rm=TRUE )/ AllP  # weighted average del value
  AllP.del.se <- sqrt( sum((Pdel.se*Ptot/AllP)^2, na.rm=TRUE) + sum((Pdel.bar*Ptot.se*(AllP-Ptot)/AllP^2)^2, na.rm=TRUE)) # Gaussian error propagation
                            # note that the second sum recognizes that each Ptot is part of AllP so they are handled jointly
  
  nonQ <- AllP - AllQ  # mass balance
  nonQ.se <- sqrt( AllP.se^2 + AllQ.se^2 )  # Gaussian error propagation (note this overstates the uncertainty b/c it assumes that P and Q are uncorrelated)
  nonQ.del <- (AllP.del*AllP - AllQ.del*AllQ)/nonQ # isotope mass balance
  
  d_d.Pdel1 <- Ptot[1]/nonQ       #derivatives of nonQ.del for error propagation
  d_d.Pdel2 <- Ptot[2]/nonQ
  d_d.AllQ.del <- -AllQ/nonQ
  d_d.Ptot1 <- (Pdel.bar[1]-nonQ.del)/nonQ
  d_d.Ptot2 <- (Pdel.bar[2]-nonQ.del)/nonQ
  d_d.AllQ  <- (nonQ.del-AllQ.del)/nonQ

  nonQ.del.se <- sqrt( (Pdel.se[1]*d_d.Pdel1)^2      #Gaussian error propagation    Equation (S22)
                       + (Pdel.se[2]*d_d.Pdel2)^2 
                       + (AllQ.del.se*d_d.AllQ.del)^2 
                       + (Ptot.se[1]*d_d.Ptot1)^2 
                       + (Ptot.se[2]*d_d.Ptot2)^2
                       + (AllQ.se*d_d.AllQ)^2 )
  
  
  Qdel.bar[nlevels(Qcat)+1] <- AllQ.del  # here we copy the results into the relevant parts of the Qdel and Qtot vectors
  Qdel.se[nlevels(Qcat)+1] <- AllQ.del.se
  Qtot[nlevels(Qcat)+1] <- AllQ
  Qtot.se[nlevels(Qcat)+1] <- AllQ.se

  Qdel.bar[nlevels(Qcat)+2] <- nonQ.del
  Qdel.se[nlevels(Qcat)+2] <- nonQ.del.se
  Qtot[nlevels(Qcat)+2] <- nonQ
  Qtot.se[nlevels(Qcat)+2] <- nonQ.se
  
  
  
  
  
  # do end-member mixing
  
  f <- array(NA, dim=c(nlevels(Qcat)+2, nlevels(Pcat)), dimnames = list(c(levels(Qcat), "AllQ", "nonQ"), paste("f.",levels(Pcat), sep=""))) # array of end-member mixing fractions
  f.se <- array(NA, dim=c(nlevels(Qcat)+2, nlevels(Pcat)), dimnames = list(c(levels(Qcat), "AllQ", "nonQ"), paste("f.",levels(Pcat),".se",sep=""))) # array of standard errors
  
  denom <- Pdel.bar[1]-Pdel.bar[2] # we will be re-using this a lot
  for(j in 1:(nlevels(Qcat)+1)){
    f[j,1] <- (Qdel.bar[j]-Pdel.bar[2])/denom # end-member mixing: Equation (3)
    f[j,2] <- 1-f[j,1]                        # because the two fractions need to add to 1
    f.se[j,1] <- sqrt( (Qdel.se[j]/denom)^2 + (Pdel.se[1]*(-f[j,1]/denom))^2 + (Pdel.se[2]*(Pdel.bar[1]-Qdel.bar[j])/denom^2)^2 )  # Equation (S12)
    f.se[j,2] <- f.se[j,1] # Gaussian error propagation
  }
  
  # fraction of ET from each precipitation source
  nj <- nlevels(Qcat)+2
  f[nj,1] <- (Ptot[1] - Qtot[nj-1] * (Qdel.bar[nj-1]-Pdel.bar[2])/denom)/(Ptot[1]+Ptot[2]-Qtot[nj-1])  # Equation (8) or (18)

  f.se[nj,1] <- sqrt( (Ptot.se[1]*(1-f[nj,1])/nonQ)^2                                                  # Equation (S20)
                      + (Ptot.se[2]*f[nj,1]/nonQ)^2
                      + (Qtot.se[nj-1]*(f[nj,1]-f[nj-1,1])/nonQ)^2
                      + (Qdel.se[nj-1]*AllQ/(Pdel.bar[1]-Pdel.bar[2])/nonQ)^2
                      + (Pdel.se[1]*AllQ*f[nj-1,1]/(Pdel.bar[1]-Pdel.bar[2])/nonQ)^2
                      + (Pdel.se[2]*AllQ*f[nj-1,2]/(Pdel.bar[2]-Pdel.bar[1])/nonQ)^2 )
  
  f[nj,2] <- 1-f[nj,1]
  f.se[nj,2] <- f.se[nj,1]
  
  

  # do end-member splitting
  
  eta <- array(NA, dim=c(nlevels(Qcat)+2, nlevels(Pcat)), dimnames = list(c(levels(Qcat), "AllQ", "nonQ"), paste("eta.",levels(Pcat), sep=""))) # array of end-member splitting fractions
  eta.se <- array(NA, dim=c(nlevels(Qcat)+2, nlevels(Pcat)), dimnames = list(c(levels(Qcat), "AllQ", "nonQ"), paste("eta.",levels(Pcat),".se", sep=""))) # array of standard errors
  
  for(j in 1:(nlevels(Qcat)+1)){
    for(i in 1:nlevels(Pcat)){
      eta[j,i] <- f[j,i]*Qtot[j]/Ptot[i]                                                                            # end-member splitting: Equation (4) or (22)-(24)
      eta.se[j,i] <- abs(eta[j,i]) * sqrt( (f.se[j,i]/f[j,i])^2 + (Qtot.se[j]/Qtot[j])^2 + (Ptot.se[i]/Ptot[i])^2)  # Gaussian error propagation: Equation (S17)
    }
  }
  
  # and now do the "nonQ" category (i.e., ET or other unmeasured fluxes)
  for(i in 1:nlevels(Pcat)){
    eta[nj,i] <- 1 - eta[nj-1,i]      # end-member splitting: Equation (24)
    eta.se[nj,i] <- eta.se[nj-1,i]     # Gaussian error propagation
  }
  
  
  table <- cbind(f, f.se, eta, eta.se)
  
  
  
  return(list(table = table,              #return table of f's and eta's
              Pdel.bar = Pdel.bar,        #return fluxes and mean tracer signatures of P and Q, with standard errors
              Pdel.se = Pdel.se,
              Qdel.bar = Qdel.bar,
              Qdel.se = Qdel.se,
              Ptot = Ptot,
              Ptot.se = Ptot.se,
              Qtot = Qtot,
              Qtot.se = Qtot.se))
  
} # end EndSplit




