library(tidyverse)
library(vegan)
library(FD)

################################################################################
# Null model which simulates randomized communities based on observed plots.   #
#                                                                              #
# Our null model fixes observed plot-level species richness, but allows for    # 
# random selection of species from the observed species pool at any relative   #
# abundance. This approach simulates communities devoid of any species         #
# filtering or assembly rules.                                                 #
#                                                                              #          
# The model is composed of two functions. The first function simulates a       #
# single randomized community and provides a null taxonomic and functional     #
# beta diversity value. Taxonomic beta diversity is calculated as the mean     #
# Bray-Curtis dissimilarity. Functional beta diversity is calculated using     #
# Rao's quadratic entropy, scaled from 0 to 1, in an additive partition (see   #          
# main text for details). The second function iterates the first function,     #
# providing  n (999) null beta diversity values.                               #
################################################################################

################################################################################
## FUNCTION ONE -- NULL BETA DIVERSITY GENERATOR ###############################
################################################################################

# x = a community matrix in long form, such that there are three columns: species, plot number, and abundance (cover). 
# trait = a trait table
# plots = plots that will be simulated (plot numbers provide species richness constraints and correspond to fire severity class: 1-20 are unburned, 21-40 are low-severity, 41-60 are high-severity)

null_beta <- function(x, traits, plots){
  
  # Generate species list from observed species pool, which is represented in the trait table
  Gamma <- as.vector(row.names(traits))
  
  # Start building data frame with plot numbers. Missing and empty plots are absent.
  nullModel <- data.frame("Plot" = unique(x$plot))
  
  # Generate species richness in each plot
  Richness <- as.vector(x %>%
    filter(cover > 0) %>%
    group_by(plot) %>%
    summarize(Richness=n()) %>%
    select(Richness))
  
  nullModel <- nullModel %>%
    cbind(Richness)
  
  # Add column for species and randomly sample from species list according to plot richness. The repeat loop ensures all species in species pool are included, such that gamma diversity remains constant.
  repeat{
    nullModel. <- nullModel %>% # nullModel becomes nullModel. so that nullModel is preserved if Gamma condition is not met.
      group_by(Plot) %>%
      mutate(Species=NA, across(Species, ~list(sample(Gamma, size=Richness, replace=F)))) %>%
      unnest(Species) %>%
      ungroup()
    
    if(length(unique(nullModel.$Species)) == length(Gamma)){break} 
    else (rm(nullModel.))
  }
  
  # Create vector of random species covers from uniform distribution (0.2 to 100) and bind to data frame.
  Cover <- as.vector(runif(n=length(nullModel.$Species), 0.2, 100))
  
  nullModel. <- nullModel. %>%
    cbind("Cover"=Cover)
  
  # Convert null community data frame into a matrix. Columns are alphabetized by species.
  nullModel. <- nullModel. %>%
    pivot_wider(names_from=Species, values_from=Cover, values_fill=0)
  
  nullMatrix <- nullModel. %>% 
    select(-Richness) %>%
    column_to_rownames(var="Plot")
  
  nullMatrix <- nullMatrix[,order(colnames(nullMatrix))]
  
  # Relativize using Wisconsin double standardization and extract plots by severity. Missing species (abundance=0 in all selected rows) are removed.
  nullMatrix <- wisconsin(nullMatrix)
  nullMatrix <- nullMatrix %>%
    filter(rownames(.) %in% plots)
  nullMatrix <- nullMatrix[, colSums(nullMatrix != 0) > 0]
  
  # Trim trait table to match community matrix
  traitMatrix <- traits %>%
    filter(rownames(.) %in% colnames(nullMatrix))
  
  # Compute taxonomic and functional beta diversity
  betaTAX <- mean(vegdist(nullMatrix, method="bray"))
  
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