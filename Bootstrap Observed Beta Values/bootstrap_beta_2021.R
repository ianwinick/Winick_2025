library(tidyverse)
library(vegan)
library(FD)

################################################################################
# R code for bootstrapping the observed taxonomic and functional beta          #
# diversity values for the 2021 community data (year 2 after fire)             #
################################################################################

# See bootstrap_beta_2020 for code comments

data <- read_csv("CommunityMatrix_2021.csv")

traits <- read_csv("traitMatrix.csv") %>%
  column_to_rownames(var="spp")

data_U <- species %>%
  filter(rownames(.) %in% 1:20)
data_U <- data_U[, colSums(data_U != 0) > 0]
data_U <- data_U[rowSums(data_U != 0) > 0, ]

data_L <- species %>%
  filter(rownames(.) %in% 21:40)
data_L <- data_L[, colSums(data_L != 0) > 0]
data_L <- data_L[rowSums(data_L != 0) > 0, ]

data_H <- species %>%
  filter(rownames(.) %in% 41:60)
data_H <- data_H[, colSums(data_H != 0) > 0]
data_H <- data_H[rowSums(data_H != 0) > 0, ]

################################################################################
# TAXONOMIC ####################################################################
################################################################################

# Unburned
dist_U <- as.data.frame(as.matrix(vegdist(data_U))) %>%
  rownames_to_column(var="y") %>%
  pivot_longer(-y, names_to="x", values_to="dist") %>%
  mutate(y=as.integer(y), x=as.integer(x)) %>%
  filter(x<y)

betaBoot = NULL
for(i in 1:999){
  betaBoot[i]=mean(sample(dist_U$dist, size=length(dist_U$dist), replace=TRUE))  
}

bootU <- data.frame("year"=2, "boot"=1:999, "type"="tax", "severity"="U", "beta"=betaBoot)

# Low Severity
dist_L <- as.data.frame(as.matrix(vegdist(data_L))) %>%
  rownames_to_column(var="y") %>%
  pivot_longer(-y, names_to="x", values_to="dist") %>%
  mutate(y=as.integer(y), x=as.integer(x)) %>%
  filter(x<y)

betaBoot = NULL
for(i in 1:999){
  betaBoot[i]=mean(sample(dist_L$dist, size=length(dist_L$dist), replace=TRUE))  
}

bootL <- data.frame("year"=2, "boot"=1:999, "type"="tax", "severity"="L", "beta"=betaBoot)

# High Severity
dist_H <- as.data.frame(as.matrix(vegdist(data_H))) %>%
  rownames_to_column(var="y") %>%
  pivot_longer(-y, names_to="x", values_to="dist") %>%
  mutate(y=as.integer(y), x=as.integer(x)) %>%
  filter(x<y)

betaBoot = NULL
for(i in 1:999){
  betaBoot[i]=mean(sample(dist_H$dist, size=length(dist_H$dist), replace=TRUE))  
}

bootH <- data.frame("year"=2, "boot"=1:999, "type"="tax", "severity"="H", "beta"=betaBoot)

################################################################################
# FUNCTIONAL ###################################################################
################################################################################

# Unburned
traitMatrix <- traits %>%
  filter(rownames(.) %in% colnames(data_U))

Q <- dbFD(traitMatrix, data_U, scale.RaoQ=T, messages=FALSE)$RaoQ
Qgamma <- dbFD(traitMatrix, colMeans(data_U), scale.RaoQ=T, messages=FALSE)$RaoQ

Qboot = NULL
for(i in 1:999){
  Qalpha=sample(Q, size=length(Q), replace=TRUE)
  Qboot[i]=mean(Qgamma-Qalpha)
}

Q_bootU <- data.frame("year"=2, "boot"=1:999, "type"="fun", "severity"="U", "beta"=Qboot)

# Low Severity
traitMatrix <- traits %>%
  filter(rownames(.) %in% colnames(data_L))

Q <- dbFD(traitMatrix, data_L, scale.RaoQ=T, messages=FALSE)$RaoQ
Qgamma <- dbFD(traitMatrix, colMeans(data_L), scale.RaoQ=T, messages=FALSE)$RaoQ

Qboot = NULL
for(i in 1:999){
  Qalpha=sample(Q, size=length(Q), replace=TRUE)
  Qboot[i]=mean(Qgamma-Qalpha)
}

Q_bootL <- data.frame("year"=2, "boot"=1:999, "type"="fun", "severity"="L", "beta"=Qboot)

# High Severity
traitMatrix <- traits %>%
  filter(rownames(.) %in% colnames(data_H))

Q <- dbFD(traitMatrix, data_H, scale.RaoQ=T, messages=FALSE)$RaoQ
Qgamma <- dbFD(traitMatrix, colMeans(data_H), scale.RaoQ=T, messages=FALSE)$RaoQ

Qboot = NULL
for(i in 1:999){
  Qalpha=sample(Q, size=length(Q), replace=TRUE)
  Qboot[i]=mean(Qgamma-Qalpha)
}

Q_bootH <- data.frame("year"=2, "boot"=1:999, "type"="fun", "severity"="H", "beta"=Qboot)

################################################################################
# Bind all of the bootstrap data frames together ###############################
################################################################################

bootstrap <- rbind(bootU, bootL, bootH, Q_bootU, Q_bootL, Q_bootH)
write.csv(bootstrap, "Bootstrap Observed Beta Values/bootstrap_2021.csv", row.names=FALSE)
