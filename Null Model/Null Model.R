library(tidyverse)
library(vegan)
library(FD)

################################################################################
# Null model which simulates randomized communities based on observed plots.   #
#                                                                              #
# Our null model fixes observed plot-level species richness and cover values,  # 
# but allows for random selection of species from the observed species pool.   #
# This approach simulates communities devoid of any species filtering or       #
# assembly rules.                                                              #
#                                                                              #          
# The model is composed of two functions. The first function simulates a       #
# single randomized community and provides a null taxonomic and functional     #
# beta diversity value. Taxonomic beta diversity is calculated as the mean     #
# Bray-Curtis dissimilarity. Functional beta diversity is calculated using     #
# Rao's quadratic entropy, scaled from 0 to 1, in an additive partition (see   #          
# main text for details). The second function iterates the first function,     #
# providing  n (999) null beta diversity values.                               #
#                                                                              #
# PLEASE NOTE!!! This code has sections designed to handle bare plots from the #
# first year after fire. These sections should be used when applied to the     #
# first year of data, but should be commented out for all other years.         #
# Comments pertaining to these sections will begin with "!!!" to help identify #
# them.                                                                        #
################################################################################

################################################################################
## FUNCTION ONE -- NULL BETA DIVERSITY GENERATOR ###############################
################################################################################

# x = a community matrix in long form, such that there are three columns: species, plot number, and abundance (cover). 
# trait = a trait table
# plots = plots that will be simulated (plot numbers provide species richness and cover constraints and correspond to fire severity class: 1-20 are unburned, 21-40 are low-severity, 41-60 are high-severity)

null_beta <- function(x, traits, plots){
  
  # Generate species list from observed species pool, which is represented in the trait table
  Gamma <- as.vector(row.names(traits))
  
  # Start building data frame with plot numbers and cover values. Missing and empty plots are absent.
  nullModel <- x %>%
    select(plot, cover)
  
  # Add column for species and randomly sample from species list according to plot richness. The repeat loop ensures all species in species pool are included, such that gamma diversity remains constant.
  repeat{
    nullModel. <- nullModel %>% # nullModel becomes nullModel. so that nullModel is preserved if Gamma condition is not met.
      group_by(plot) %>%
      mutate(Species=NA) %>%
      mutate(across(Species, ~sample(Gamma, size=length(plot), replace=F))) %>%
      ungroup()

    if(length(unique(nullModel.$Species)) == length(Gamma)){break}
    else (rm(nullModel.))
  }
  
  # !!! Attach bare plots for 2020 data. Comment this out for all other years
  # bare <- data.frame("plot"=c(33,42,45,46,48,51,53,56,57), cover=NA, Species=NA)
  # nullModel. <- nullModel. %>% rbind(bare)
  # nullModel. <- nullModel.[order(nullModel.$plot),]
  
  # Convert null community data frame into a matrix. Columns are alphabetized by species.
  nullModel. <- nullModel. %>%
    pivot_wider(names_from=Species, values_from=cover, values_fill=0)
  
  nullMatrix <- nullModel. %>% 
    # select(-"NA") %>% # !!! Comment this out for all other years after 2020
    column_to_rownames(var="plot")
  
  nullMatrix <- nullMatrix[,order(colnames(nullMatrix))]
  
  # Extract plots by severity. Missing species (abundance=0 in all selected rows) are removed.
  nullMatrix <- nullMatrix %>%
    filter(rownames(.) %in% plots)
  nullMatrix <- nullMatrix[, colSums(nullMatrix != 0) > 0]
  
  # Trim trait table to match community matrix
  traitMatrix <- traits %>%
    filter(rownames(.) %in% colnames(nullMatrix))
  
  # Compute taxonomic and functional beta diversity
  bray <- nullMatrix
    # mutate(z=0.001) # !!! Add dummy species for 2020. Remove for all other years.
  betaTAX <- mean(vegdist(bray))
  
  # !!! Use this version for 2020
  # rao <- dbFD(traitMatrix, nullMatrix[!rownames(nullMatrix) %in% bare$plot,], scale.RaoQ=T, message=F)$RaoQ
  # gammaFUN <- dbFD(traitMatrix, colMeans(nullMatrix), scale.RaoQ=T, message=F)$RaoQ
  # alphaFUN <- mean(c(rao, rep(0, 20-length(rao))))
  # betaFUN <- gammaFUN-alphaFUN
  
  # !!! Use this version for all other years
  gammaFUN <- dbFD(traitMatrix, colMeans(nullMatrix), scale.RaoQ=T, message=F)$RaoQ
  alphaFUN <- mean(dbFD(traitMatrix, nullMatrix, scale.RaoQ=T, message=F)$RaoQ)
  betaFUN <- gammaFUN-alphaFUN

  beta <- list("tax"=betaTAX, "fun"=betaFUN)
  
  return(beta)
  
}

################################################################################
# FUNCTION TWO -- ITERATE ######################################################
################################################################################

# x, traits, and plots are as above
# Nsim = number of desired iterations of the null model (default=999)

nulling_me_softly <- function(x, traits, plots, Nsim=999){
  
  tax=NULL
  fun=NULL
  
  for(i in 1:Nsim){
    y <- null_beta(x, traits, plots)
    tax[i]=y$tax
    fun[i]=y$fun
  }
  nullBeta <- cbind(tax,fun)
  return(nullBeta)
}

################################################################################
################################ END NULL MODEL ################################
################################################################################