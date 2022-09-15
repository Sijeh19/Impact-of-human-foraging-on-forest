#-------------------------------------------------------------------------------
# Impact of human foraging on tree diversity, composition and abundance in a 
# tropical rainforest
# Sijeh A. Asuk, Thomas J. Matthews, Jonathan Sadler, Thomas A. M. Pugh, 
# Vincent T. Ebu, Nzube M. Ifebueme and Nicholas Kettridge 

#-------------------------------------------------------------------------------
# Clear the objects from R to start afresh 
rm(list=ls())

# clear plots output area
dev.off()

#-------------------------------------------------------------------------------
### install Required Packages

#-------------------------------------------------------------------------------
# Load required packages
library(stringi)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggmap)
library(gridExtra)
require(tidyverse)
library(plotly)
library(reshape2)
library(fields)
library(betapart)
library(gambin)
library(sads)
library(r2symbols)
library(hrbrthemes)
library(viridis)
library(ggpubr)

#-------------------------------------------------------------------------------
# set working directory
setwd("--")

#-------------------------------------------------------------------------------
####---NOTE----####
## "tree_data" is not provided.
### only summarized version of the data has been made available

## Import the csv file with into R
# tree_data <- read.csv("Oban_data.csv", header = TRUE, stringsAsFactors = FALSE)

## sort based on edible from community record, edible from secondary data
## and, non-human categories

# edible <- filter(tree_data, Community == "Edible")
# inedible <- filter(tree_data, Community != "Edible") #### or use this 


#------------------------------------------------------------------------------
### 3.2	SPECIES ABUNDANCE DISTRIBUTIONS AND EVENNESS
#-------------------------------------------------------------------------------

###############################################################
###SAD plots
library(gambin)
library(sads)

################
## generate values from a gambin distribution from the frequency of species name.

###################
### Total dataset
all_species <- read.table("all_species_gambin.csv", header = TRUE, sep = ",", row.names = 1)
ab <- table(all_species$x) %>% as.vector()
f <- fit_abundances(ab)

##Edible species
edible_species <- read.table("edible_species_gambin.csv", header = TRUE, sep = ",", row.names = 1)
ab2 <- table(edible_species$x) %>% as.vector()
f2 <- fit_abundances(ab2)

##Not edible
inedible_species <- read.table("inedible_species_gambin.csv", header = TRUE, sep = ",", row.names = 1)
ab3 <- table(inedible_species$x) %>% as.vector()
f3 <- fit_abundances(ab3)

## generate alpha value and confidence intervals for the plots
summary(f3)
summary(f2)
summary(f)



###############################################################
###Pielou's evenness index
##########
library(tidyverse)
library(vegan)

# Total
# ----------------------------------------------------------------
# load species abunance data for total forest
tot_data <- read.csv("tot_abundance.csv")
# compute species evenness per plot
tot_data_2 <- tot_data[,-1]
# calculate Pielou's evenness
H_tot <- diversity(tot_data_2)
J_tot <- H_tot/log(specnumber(tot_data_2))

# Edible
# ----------------------------------------------------------------
# load species abunance data for edible species category
ed_data <- read.csv("ed_abundance.csv")
# compute species evenness per plot
ed_data_2 <- ed_data[,-1]
# calculate Pielou's evenness
H_ed <- diversity(ed_data_2)
J_ed <- H_ed/log(specnumber(ed_data_2))

# Inedible 
# load species abunance data for inedible species category
ined_data <- read.csv("ined_abundance.csv")
# compute species evenness per plot
ined_data_2 <- ined_data[,-1]
# calculate Pielou's evenness
H_ined <- diversity(ined_data_2)
J_ined <- H_ed/log(specnumber(ined_data_2))

# import Combined evenness data from total forest, edible and inedible species into
total <- read.csv("combined_evenness_data.csv")

# compute student's t-test to compare edible and inedible
library(car)
test <- t.test(J_ed, J_ined,
               var.equal = T, alternative = "greater"
)
test


##########################
# Compute pooled evenness
pooled_data <- read.csv("pooled_total.csv")
pooled_data_2 <- pooled_data[,-1]

# calculate Pielou's evenness
H_pool <- diversity(pooled_data_2)
J_pool <- H_pool/log(specnumber(pooled_data_2))


###############################################################
###plot SADs and Pielou's evenness index

## vitualize the SADS and Pielou's evenness
windowsFonts(A = windowsFont("Times New Roman"))

palette.colors(palette = "Okabe-Ito")


par(mfrow = c(1, 2))

# box plot
{boxplot(evenness~category,data=total, main=NA,
        xlab="Species Category", ylab="Pielou's Evenness Index", family = "A", cex.lab=1.5,
        cex.main = 1.3, cex.axis = 1.4, col=c("#009E73", "#D55E00", "#0072B2"), notch=TRUE)
mtext('(a)', side=3, at=0.4, line=0.8, col = "black", cex=1.5)
text(2,0.99, "Edible~Inedible (",cex=1.3,family = "A") 
text(3,0.99, substitute(paste(italic('p value'), ' < 0.000)')),cex=1.3,family = "A") 
}

{barplot(f, ylim=c(0, 45), xlab = "Abundance Octaves", family = "A", cex.lab=1.5,
         cex.main = 1.3, cex.axis = 1.4, col="#0072B2")
  barplot(f3, ylim=c(0, 32), family = "A", cex.lab=1.5,
          cex.main = 1.3, cex.axis = 1.4, col="#D55E00", add=T)
  barplot(f2, ylim=c(0, 31), family = "A", cex.lab=1.5,
          cex.main = 1.3, cex.axis = 1.4, col = "#009E73", add=T)
  
  points(f, pch = 16, col="#000000")
  points(f2, pch = 16, col="#000000")
  points(f3, pch = 16, col="#000000")
  lines(f, pch = 10, col="#000000")
  lines(f2, pch = 10, col="#000000")
  lines(f3, pch = 10, col="#000000")
  mtext('(b)', side=3, at=-0.1,line=0.7,col = "black", cex=1.5)
  
  legend("topright",
         c("Gambin model",expression("Total("*alpha*"=3.3, CI:2.4-4.4)"), 
           expression("Inedible("*alpha*"=3.1, CI:2.2-4.4)"), 
           expression("Edible("*alpha*"=4.0, CI:2.1-7.1)")), 
         lty=c(1, 0, 0, 0), lwd =c(3,3,2.5,3), col=c("#000000","#0072B2","#D55E00","#009E73"), pch = c(16,15,15,15))
}


###
#######----computing the standardized alpha for inedible category---###########

#First work out which of the total communities (edible or inedible) has the 
#smallest number of individuals, and set this to M, e.g.


count(inedible_species)
M <- 306

#Then set abundances to be the abundance vector of the community which has the 
#largest number of individuals. For example, if M came from edible, abundances 
#would be the abundance vector for inedible

abundances <- ab3

rep <- replicate(999, gambin::fit_abundances(abundances, subsample = M)$alpha)

mean(rep)
sd(rep)

summary(rep)

# get CI at 95%
quantile(rep,.025)
quantile(rep,.975)

# 1.37 - 2.41 

#This then gives you mean alpha (and its SD) for the SAD for the community 
#with more individuals standardised for sample size
#labels = c("1", "2-3", "4-7", "8-15", "16-31", "32-63")


#-------------------------------------------------------------------------------
### 3.3	THE TREND IN BETA DIVERSITY WITH INCREASING ELEVATION
#-------------------------------------------------------------------------------

#################################################################################################
##beta-diversity plots - trickier as first we need to make presence-absence matrix
#where rows are sites (plots) and columns are species

#### ---how sads was computed---

## 1. make PA matrix Total forest, Edible, and Inedible species

# un <- unique(tree_data$Spp_code)#species
# pn <- unique(tree_data$Plot_no)
# m <- matrix(0, ncol = length(un), nrow = length(pn))
# colnames(m) <- un
# rownames(m) <- pn

## 2. now fill in PA matrix, with a 1 for a species presence and 0 for absence

# for (i in 1:length(pn)){
# d <- filter(tree_data, Plot_no == pn[i])
# dun <- unique(d$Spp_code)
# w <- which(colnames(m) %in% dun)
# m[i,w] <- 1
# }
# any(rowSums(m) == 0)
# any(colSums(m) == 0)

## 3. compute elevation means per plot

# mean_alt2 <- tree_data %>% 
# group_by(Plot_no) %>% 
# summarise(mean(Alt))

### --- test plot using data provided ----- ####
### import required data 

m <- read.csv("Presence_absence_plot_matrix.csv", header = TRUE, stringsAsFactors = FALSE)
me <- read.csv("Presence_absence_plot_matrix_edible.csv", header = TRUE, stringsAsFactors = FALSE)
mi <- read.csv("Presence_absence_plot_matrix_inedible.csv", header = TRUE, stringsAsFactors = FALSE)


# import elevation data
elevation <- read.table("elevation and plot distance.csv", header = TRUE, sep = ",", row.names = 1)

mean_alt <- unlist(elevation[,2]) %>% as.vector()
mean_alte <- unlist(elevation[,3]) %>% as.vector()
mean_alti <- unlist(elevation[,4]) %>% as.vector()


#calculate total beta-div
#install.packages("betapart")
library(betapart)
bm <- beta.pair(m, index.family = "sorensen") 
dist.DI=dist(mean_alt, method = "euclidean", diag = FALSE, upper = FALSE)


###############################################################
###############################################################
##beta-diversity plots - edible and in edible
## 1. make PA matrix Total forest, Edible, and Inedible species

## 2. now fill in PA matrix, with a 1 for a species presence and 0 for absence

## 3. compute elevation means per plot

#calculate beta-div for edible
#install.packages("betapart")
library(betapart)
bme <- beta.pair(me, index.family = "sorensen") 
dist.DIe=dist(mean_alte, method = "euclidean", diag = FALSE, upper = FALSE)


#calculate beta-div for inedible
#install.packages("betapart")
bmi <- beta.pair(mi, index.family = "sorensen") 
dist.DIi=dist(mean_alti, method = "euclidean", diag = FALSE, upper = FALSE)


library(ecodist)
##M antel correlation test for total forest
## sor component
ecodist::mantel(bm$beta.sor ~ dist.DI, nperm = 1000, nboot = 500, pboot = 0.9, cboot = 0.95)
## sim component
ecodist::mantel(bm$beta.sim ~ dist.DI, nperm = 1000, nboot = 500, pboot = 0.9, cboot = 0.95)


##Mantel correlation test for edible
## sor component
ecodist::mantel(bme$beta.sor ~ dist.DIi, nperm = 1000, nboot = 500, pboot = 0.9, cboot = 0.95)
## sim component
ecodist::mantel(bme$beta.sim ~ dist.DIi, nperm = 1000, nboot = 500, pboot = 0.9, cboot = 0.95)


##Mantel correlation test for inedible
## sor component
ecodist::mantel(bmi$beta.sor ~ dist.DIi, nperm = 1000, nboot = 500, pboot = 0.9, cboot = 0.95)
## sim component
ecodist::mantel(bmi$beta.sim ~ dist.DIi, nperm = 1000, nboot = 500, pboot = 0.9, cboot = 0.95)


##################################################
## plots beta diversity graphs for all attributes
{par(mfrow = c(2,3))
  par(mai = c(0.5, 0.7, 0.3, 0.17)) 
  
  windowsFont("Times New Roman")
  
  ###all forest BSor
  plot(bm$beta.sor~dist.DI, pch=16, main = "All species",
       ylim = c(0.3, 0.8), xlab = NA, ylab = "BSor", cex.lab=2.5,
       cex.main = 2.5, cex.axis = 2, cex=1.5, col="brown")
  sor.lm=lm(bm$beta.sor~dist.DI)
  abline(sor.lm, lwd=2)
  text(130, 0.77, "r = 0.43*** (95%CI=0.36-0.49)", cex=2)
  mtext('(a)', side=3, line=0.6, col = "black", at=last(-55), cex=1.5)
  
  
  ###edible BSor
  plot(bme$beta.sor~dist.DIe, pch=19,  main = "Edible species",
       ylim = c(0.3, 0.9), xlab = NA, ylab = NA, cex.lab=2.5,
       cex.main = 2.5, cex.axis = 2, cex=1.5, col="dark green")
  sor.lme=lm(bme$beta.sor~dist.DIe)
  abline(sor.lme, lwd=2)
  text(140, 0.86, "r = 0.04(95% CI = -0.06-0.15)", cex=2)
  mtext('(b)', side=3, line=0.6, col = "black", at=last(-55), cex=1.5)
  
  
  ###inedible BSOR
  plot(bmi$beta.sor~dist.DIi, pch=19,  main = "Inedible species",
       ylim = c(0.3, 0.9), xlab = NA, ylab = NA, cex.lab=2.5,
       cex.main = 2.5, cex.axis = 2, cex=1.5)
  sor.lmi=lm(bmi$beta.sor~dist.DIi)
  abline(sor.lmi, lwd=2)
  text(140, 0.86, "r = 0.47**(95% CI = 0.36 - 0.55)", cex=2)
  mtext('(c)', side=3, line=0.6, col = "black", at=last(-55), cex=1.5)
  
  par(mai = c(0.72, 0.7, 0.1, 0.17)) 
  
  ### all forest BSim
  plot(bm$beta.sim~dist.DI, pch=16, 
       ylim = c(0.3, 0.8), xlab = "Difference in Elevation (m)", 
       ylab = "BSim", cex.lab=2.5,
       cex.main = 2.5, cex.axis = 2, cex=1.5, col="brown")
  sor.lm=lm(bm$beta.sim~dist.DI)
  abline(sor.lm, lwd=2)
  text(130, 0.77, "r = 0.42** (95%CI=0.35-0.50)", cex=2)
  mtext('(d)', side=3, line=0.6, col = "black", at=last(-55), cex=1.5)
  
  
  ###edible BSim
  plot(bme$beta.sim~dist.DIe, pch=19, ylim = c(0.3, 0.9), 
       xlab = "Difference in Elevation (m)", ylab = NA, cex.lab=2.5,
       cex.main = 2.5, cex.axis = 2, cex=1.5, col="dark green")
  sor.lme=lm(bme$beta.sim~dist.DIe)
  abline(sor.lme, lwd=2)
  text(140, 0.86, "r = -0.002(95% CI = -0.12 - 0.14)", cex=2)
  mtext('(e)', side=3, line=0.6, col = "black", at=last(-55), cex=1.5)
  
  
  ###inedible BSim
  plot(bmi$beta.sim~dist.DIi, pch=19, ylim = c(0.3, 0.9),
       xlab = "Difference in Elevation (m)", ylab = NA, cex.lab=2.5,
       cex.main = 2.5, cex.axis = 2, cex=1.5)
  sor.lmi=lm(bmi$beta.sim~dist.DIi)
  abline(sor.lmi, lwd=2)
  text(140, 0.86, "r = 0.45***(95% CI = 0.36 - 0.53)", cex=2)
  mtext('(f)', side=3, line=0.6, col = "black", at=last(-55), cex=1.5)
}

#-------------------------------------------------------------------------------
### 3.4.3 	NULL MODEL ANALYSIS

################################################################
###NULL MODELLING##################################
##########################################################

#NULL MODEL 1: standard fixed-fixed null model.

#this is custom function which takes a presence-absence matrix, calculates beta-diversty and then does the
#mantel correlation
bd_funsor <- function(x, mean_alt){
  bme <- beta.pair(x, index.family = "sorensen")
  dist.DI=dist(mean_alt, method = "euclidean", diag = FALSE, upper = FALSE)
  vv <- vegan::mantel(bm$beta.sor,dist.DI)
  vv2 <- vv$statistic
  names(vv2) <- "r"
  vv2
}

bd_funsim <- function(x, mean_alt){
  bme <- beta.pair(x, index.family = "sorensen")
  dist.DI=dist(mean_alt, method = "euclidean", diag = FALSE, upper = FALSE)
  vv <- vegan::mantel(bm$beta.sim,dist.DI)
  vv2 <- vv$statistic
  names(vv2) <- "r"
  vv2
}

#this then does the null modelling (999 iterations)
#you will want to report the SES and the P-vale (PR(sim.))
#a positive SEM indicates you have a larger statistic (here correlation coefficient r) than
#expected, and if P is < 0.05, this is significantly more than, and vice versa.
#have a read about how quasiswap works
rvsor <- vegan::oecosimu(m, bd_funsor, mean_alt = mean_alt, method = "quasiswap", nsimul = 999, statistic = "r")
rvsor

rvsim <- vegan::oecosimu(m, bd_funsim, mean_alt = mean_alt, method = "quasiswap", nsimul = 999, statistic = "r")
rvsim


######################################################
##Fixed-fixed null model for EDIBLE Species
bd_funsore <- function(x, mean_alt){
  bme <- beta.pair(x, index.family = "sorensen")
  dist.DIe=dist(mean_alt, method = "euclidean", diag = FALSE, upper = FALSE)
  vve <- vegan::mantel(bme$beta.sor,dist.DIe)
  vve2 <- vve$statistic
  names(vve2) <- "r"
  vve2
}

bd_funsime <- function(x, mean_alt){
  bme <- beta.pair(x, index.family = "sorensen")
  dist.DIe=dist(mean_alt, method = "euclidean", diag = FALSE, upper = FALSE)
  vve <- vegan::mantel(bme$beta.sim,dist.DIe)
  vve2 <- vve$statistic
  names(vve2) <- "r"
  vve2
}

#this then does the null modelling (999 iterations)
rvsore <- vegan::oecosimu(me, bd_funsore, mean_alt = mean_alt, method = "quasiswap", nsimul = 999, statistic = "r")
rvsore

rvsime <- vegan::oecosimu(me, bd_funsime, mean_alt = mean_alt, method = "quasiswap", nsimul = 999, statistic = "r")
rvsime

#############################################
##Fixed-fixed null model for INEDIBLE Species
bd_funsori <- function(x, mean_alt){
  bmi <- beta.pair(x, index.family = "sorensen")
  dist.DIi=dist(mean_alt, method = "euclidean", diag = FALSE, upper = FALSE)
  vvi <- vegan::mantel(bmi$beta.sor,dist.DIi)
  vvi2 <- vvi$statistic
  names(vvi2) <- "r"
  vvi2
}

bd_funsimi <- function(x, mean_alt){
  bmi <- beta.pair(x, index.family = "sorensen")
  dist.DIi=dist(mean_alt, method = "euclidean", diag = FALSE, upper = FALSE)
  vvi <- vegan::mantel(bmi$beta.sim,dist.DIi)
  vvi2 <- vvi$statistic
  names(vvi2) <- "r"
  vvi2
}

#this then does the null modelling (999 iterations)
rvsori <- vegan::oecosimu(mi, bd_funsori, mean_alt = mean_alt, method = "quasiswap", nsimul = 999, statistic = "r")
rvsori

rvsimi <- vegan::oecosimu(mi, bd_funsimi, mean_alt = mean_alt, method = "quasiswap", nsimul = 999, statistic = "r")
rvsimi


######################################
#NULL MODEL 2: randomising the edible - inedible classifcaitons
######################################################

####
###EDIBLE
####################
RANDSore <- function(tree_data){
  
  #here we randomise the h, a and x values across space
  all_spp2 <- tree_data
  all_spp2$Community <- sample(all_spp2$Community, nrow(all_spp2), replace = FALSE)
  
  #remember to change this for inedible species and vice versa
  edible <- filter(all_spp2, Community == "Edible") #edible
  #edible <- filter(all_spp, Community %in% c("a", "x")) #inedible
  
  #make PA matrix
  une <- unique(edible$Spp_code)#species
  pne <- unique(edible$Plot_no)
  me <- matrix(0, ncol = length(une), nrow = length(pne))
  colnames(me) <- une
  rownames(me) <- pne
  
  #now fill in PA matrix, with a 1 for a species presence
  for (i in 1:length(pne)){
    de <- filter(edible, Plot_no == pne[i])
    dune <- unique(de$Spp_code)
    we <- which(colnames(me) %in% dune)
    me[i,we] <- 1
  }
  
  #get mean plot altitude
  mean_alt2 <- tree_data %>% 
    group_by(Plot_no) %>% 
    summarise(mean(Alt))
  
  mean_alt <- unlist(mean_alt2[,2]) %>% as.vector()
  
  bme <- beta.pair(me, index.family = "sorensen")
  dist.DIe=dist(mean_alt, method = "euclidean", diag = FALSE, upper = FALSE)
  vve <- vegan::mantel(bme$beta.sor, dist.DIe) #change to Bsim
  vve2 <- vve$statistic
  return(vve2)
  
}

##then run that function 999 times and calculate the SES and P-value
REPSore <- replicate(999, RANDSore(all_spp2), simplify = "vector")
obsSore <- 0.327 #THIS IS THE OBSERVED MANTEL CORRELATION FOR edible species and overall beta (you will need to change it)
SESSore <- (obsSore - mean(REPSore)) / sd(REPSore)
SESSore
2*pnorm(-abs(SESSore)) #p-value


#again you will need to edit the function and repeat this for bm$beta.sim (just change the one line), 
#and then for inedible species (both sor and sim). For the latter, you only have to change
#that one line in the function - edible <- ....


###############################################################################################################

RANDSime <- function(tree_data){
  
  #here we randomise the h, a and x values across space
  all_spp2 <- tree_data
  all_spp2$Community <- sample(all_spp2$Community, nrow(all_spp2), replace = FALSE)
  
  #remember to change this for inedible species and vice versa
  edible <- filter(tree_data, Community == "Edible") #edible
  #edible <- filter(all_spp, Community %in% c("a", "x")) #inedible
  
  #make PA matrix
  une <- unique(edible$Spp_code)#species
  pne <- unique(edible$Plot_no)
  me <- matrix(0, ncol = length(une), nrow = length(pne))
  colnames(me) <- une
  rownames(me) <- pne
  
  #now fill in PA matrix, with a 1 for a species presence
  for (i in 1:length(pne)){
    de <- filter(edible, Plot_no == pne[i])
    dune <- unique(de$Spp_code)
    we <- which(colnames(me) %in% dune)
    me[i,we] <- 1
  }
  
  #get mean plot altitude
  mean_alt2 <- all_spp2 %>% 
    group_by(Plot_no) %>% 
    summarise(mean(Alt))
  
  mean_alt <- unlist(mean_alt2[,2]) %>% as.vector()
  
  bme <- beta.pair(me, index.family = "sorensen")
  dist.DIe=dist(mean_alt, method = "euclidean", diag = FALSE, upper = FALSE)
  vve <- vegan::mantel(bme$beta.sim, dist.DIe) #change to Bsim
  vve2 <- vve$statistic
  return(vve2)
  
}

##then run that function 999 times and calculate the SES and P-value
REPSime <- replicate(999, RANDSime(all_spp2), simplify = "vector")
obsSime <- 0.3218  #THIS IS THE OBSERVED MANTEL CORRELATION FOR edible species and overall beta (you will need to change it)
SESSime <- (obsSime - mean(REPSime)) / sd(REPSime)
SESSime
2*pnorm(-abs(SESSime)) #p-value


##################################################
#INEDIBLE
#####################################################

RANDSori <- function(tree_data){
  
  #here we randomise the y and n values across space
  all_spp2 <- tree_data
  all_spp2$Community <- sample(all_spp2$Community, nrow(all_spp2), replace = FALSE)
  
  #remember to change this for inedible species and vice versa
  inedible <- filter(all_spp2, Community == "Inedible") #edible
  #edible <- filter(all_spp, Community %in% c("a", "x")) #inedible
  
  #make PA matrix
  uni <- unique(inedible$Spp_code)#species
  pni <- unique(inedible$Plot_no)
  mi <- matrix(0, ncol = length(uni), nrow = length(pni))
  colnames(mi) <- uni
  rownames(mi) <- pni
  
  #now fill in PA matrix, with a 1 for a species presence
  for (i in 1:length(pni)){
    di <- filter(inedible, Plot_no == pni[i])
    duni <- unique(di$Spp_code)
    wi <- which(colnames(mi) %in% duni)
    mi[i,wi] <- 1
  }
  
  #get mean plot altitude
  mean_alt2 <- all_spp2 %>% 
    group_by(Plot_no) %>% 
    summarise(mean(Alt))
  mean_alt <- unlist(mean_alt2[,2]) %>% as.vector()
  
  bmi <- beta.pair(mi, index.family = "sorensen")
  dist.DIi=dist(mean_alt, method = "euclidean", diag = FALSE, upper = FALSE)
  vvi <- vegan::mantel(bmi$beta.sor, dist.DIi) #change to Bsim
  vvi2 <- vvi$statistic
  return(vvi2)
  
}

##then run that function 999 times and calculate the SES and P-value
REPSori <- replicate(999, RANDSori(all_spp2), simplify = "vector")
obsSori <- 0.47083 #THIS IS THE OBSERVED MANTEL CORRELATION FOR edible species and overall beta (you will need to change it)
SESSori <- (obsSori - mean(REPSori)) / sd(REPSori)
SESSori
2*pnorm(-abs(SESSori)) #p-value

RANDSimi <- function(tree_data){
  
  #here we randomise the h, a and x values across space
  all_spp2 <- tree_data
  all_spp2$Community <- sample(all_spp2$Community, nrow(all_spp2), replace = FALSE)
  
  #remember to change this for inedible species and vice versa
  inedible <- filter(all_spp2, Community == "Inedible") #edible
  #edible <- filter(all_spp, Community %in% c("a", "x")) #inedible
  
  #make PA matrix
  uni <- unique(inedible$Spp_code)#species
  pni <- unique(inedible$Plot_no)
  mi <- matrix(0, ncol = length(uni), nrow = length(pni))
  colnames(mi) <- uni
  rownames(mi) <- pni
  
  #now fill in PA matrix, with a 1 for a species presence
  for (i in 1:length(pni)){
    di <- filter(inedible, Plot_no == pni[i])
    duni <- unique(di$Spp_code)
    wi <- which(colnames(mi) %in% duni)
    mi[i,wi] <- 1
  }
  
  #get mean plot altitude
  mean_alt2 <- all_spp2 %>% 
    group_by(Plot_no) %>% 
    summarise(mean(Alt))
  mean_alt <- unlist(mean_alt2[,2]) %>% as.vector()
  
  bmi <- beta.pair(mi, index.family = "sorensen")
  dist.DIi=dist(mean_alt, method = "euclidean", diag = FALSE, upper = FALSE)
  vvi <- vegan::mantel(bmi$beta.sim, dist.DIi) #change to Bsim
  vvi2 <- vvi$statistic
  return(vvi2)
  
}

##then run that function 999 times and calculate the SES and P-value
REPSimi <- replicate(999, RANDSimi(inedible), simplify = "vector")
obsSimi <- 0.45007 #THIS IS THE OBSERVED MANTEL CORRELATION FOR edible species and overall beta (you will need to change it)
SESSimi <- (obsSimi - mean(REPSimi)) / sd(REPSimi)
SESSori
2*pnorm(-abs(SESSimi)) #p-value



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## SENSITIVITY TEST PLOTS
## relationship between beta diversity and altitude for 12 plots
### presence-absence data (for m_12, me_12 and mi_12) has been provided.

## we filtered the first 3 plots out of the data to check for distance effect
# data_name <- filter(all_data, Plot_no >= 4)

## we then filtered edible and inedible species for the 12 plots 
# edible <- filter(data_name, Community == "Edible")
# inedible <- filter(data_name, Community != "Edible") #### or use this 

#################################################################################################
##beta-diversity plots - trickier as first we need to make presence-absence (PA) matrix
#where rows are sites (plots) and columns are species

### PA matrix explanation

### 1. we created the PA matrix using the code below

# un <- unique(data_name$Spp_code)#species
# pn <- unique(data_name$Plot_no)
# m_12 <- matrix(0, ncol = length(un), nrow = length(pn))
# colnames(m_12) <- un
# rownames(m_12) <- pn

# 2. We now filled in PA matrix, with 1 for a species presence and 0 for absence

# for (i in 1:length(pn)){
# d <- filter(data_name, Plot_no == pn[i])
# dun <- unique(d$Spp_code)
# w <- which(colnames(m_12) %in% dun)
# m_12[i,w] <- 1
# }
# any(rowSums(m_12) == 0)
# any(colSums(m_12) == 0)

## 3. compute elevation means per plot for each category using 
## the code below
# mean_alt2 <- data_name %>% 
#  group_by(variable) %>% 
#  summarise(mean(elevation))

# mean_alt <- unlist(mean_alt2[,2]) %>% as.vector()


##### ----- PLEASE, START HERE  -------- ########
#### download and import presence absence data that has been provided in the 
#### supplementary information

m_12 <- read.csv("Presence_absence_plot_matrix_12.csv", header = TRUE, stringsAsFactors = FALSE)
me_12 <- read.csv("Presence_absence_plot_matrix_12_edible.csv", header = TRUE, stringsAsFactors = FALSE)
mi_12 <- read.csv("Presence_absence_plot_matrix_12_inedible.csv", header = TRUE, stringsAsFactors = FALSE)

# get mean plot altitude
elevation_12 <- filter(elevation, Plot_no >= 4)

mean_alt_12 <- unlist(elevation_12[,2]) %>% as.vector()
mean_alte_12 <- unlist(elevation_12[,3]) %>% as.vector()
mean_alti_12 <- unlist(elevation_12[,4]) %>% as.vector()

#calculate beta-div
#install.packages("betapart")
library(betapart)

## total 
bm_12 <- beta.pair(m_12, index.family = "sorensen") 
dist.DI_12=dist(mean_alt_12, method = "euclidean", diag = FALSE, upper = FALSE)

## edible
bme_12 <- beta.pair(me_12, index.family = "sorensen") 
dist.DIe_12=dist(mean_alte_12, method = "euclidean", diag = FALSE, upper = FALSE)

## Inedible 
bmi_12 <- beta.pair(mi_12, index.family = "sorensen") 
dist.DIi_12=dist(mean_alti_12, method = "euclidean", diag = FALSE, upper = FALSE)


##Mantel correlation test total forest
ecodist::mantel(bm_12$beta.sor ~ dist.DI_12, nperm = 1000, nboot = 500, pboot = 0.9, cboot = 0.95)
ecodist::mantel(bm_12$beta.sim ~ dist.DI_12, nperm = 1000, nboot = 500, pboot = 0.9, cboot = 0.95)

##Mantel correlation test edible
ecodist::mantel(bme_12$beta.sor ~ dist.DIe_12, nperm = 1000, nboot = 500, pboot = 0.9, cboot = 0.95)
ecodist::mantel(bme_12$beta.sim ~ dist.DIe_12, nperm = 1000, nboot = 500, pboot = 0.9, cboot = 0.95)

##Mantel correlation test inedible
ecodist::mantel(bmi_12$beta.sor ~ dist.DIi_12, nperm = 1000, nboot = 500, pboot = 0.9, cboot = 0.95)
ecodist::mantel(bmi_12$beta.sim ~ dist.DIi_12, nperm = 1000, nboot = 500, pboot = 0.9, cboot = 0.95)


{par(mfrow = c(2,3))
  par(mai = c(0.5, 0.7, 0.3, 0.17)) 
  
  windowsFont("Times New Roman")
  
  ###all forest BSor
  plot(bm_12$beta.sor~dist.DI_12, pch=16, main = "All species",
       ylim = c(0.3, 0.8), xlab = NA, ylab = "BSor", cex.lab=2.5,
       cex.main = 2.5, cex.axis = 2, cex=1.5, col="brown")
  sor.lm=lm(bm_12$beta.sor~dist.DI_12)
  abline(sor.lm, lwd=2)
  text(130, 0.77, "r = 0.47**(95% CI = 0.35 - 0.56)", cex=2)
  mtext('(a)', side=3, line=0.6, col = "black", at=last(-42), cex=1.5)
  
  ###edible BSor
  plot(bme_12$beta.sor~dist.DIe_12, pch=19,  main = "Edible species",
       ylim = c(0.3, 0.9), xlab = NA, ylab = NA, cex.lab=2.5,
       cex.main = 2.5, cex.axis = 2, cex=1.5, col="dark green")
  sor.lme=lm(bme_12$beta.sor~dist.DIe_12)
  abline(sor.lme, lwd=2)
  text(140, 0.86, "r = 0.02(95% CI = -0.15 - 0.22)", cex=2)
  mtext('(b)', side=3, line=0.6, col = "black", at=last(-42), cex=1.5)
  
  ###inedible BSOR
  plot(bmi_12$beta.sor~dist.DIi_12, pch=19,  main = "Inedible species",
       ylim = c(0.3, 0.9), xlab = NA, ylab = NA, cex.lab=2.5,
       cex.main = 2.5, cex.axis = 2, cex=1.5)
  sor.lmi=lm(bmi_12$beta.sor~dist.DIi_12)
  abline(sor.lmi, lwd=2)
  text(140, 0.86, "r = 0.56***(95% CI = 0.49 - 0.64)", cex=2)
  mtext('(c)', side=3, line=0.6, col = "black", at=last(-42), cex=1.5)
  
  
  
  par(mai = c(0.72, 0.7, 0.1, 0.17)) 
  
  ### all forest BSim
  plot(bm_12$beta.sim~dist.DI_12, pch=16, 
       ylim = c(0.3, 0.8), xlab = "Difference in Elevation (m)", 
       ylab = "BSim", cex.lab=2.5,
       cex.main = 2.5, cex.axis = 2, cex=1.5, col="brown")
  sor.lm=lm(bm_12$beta.sim~dist.DI_12)
  abline(sor.lm, lwd=2)
  text(130, 0.77, "r = 0.48***(95%CI = 0.36 - 0.59)", cex=2)
  mtext('(d)', side=3, line=0.6, col = "black", at=last(-42), cex=1.5)
  
  
  ###edible BSim
  plot(bme_12$beta.sim~dist.DIe_12, pch=19, ylim = c(0.3, 0.9), 
       xlab = "Difference in Elevation (m)", ylab = NA, cex.lab=2.5,
       cex.main = 2.5, cex.axis = 2, cex=1.5, col="dark green")
  sor.lme=lm(bme_12$beta.sim~dist.DIe_12)
  abline(sor.lme, lwd=2)
  text(140, 0.86, "r = -0.03(95% CI = -0.18 - 0.19)", cex=2)
  mtext('(e)', side=3, line=0.6, col = "black", at=last(-42), cex=1.5)
  
  
  ###inedible BSim
  plot(bmi_12$beta.sim~dist.DIi_12, pch=19, ylim = c(0.3, 0.9),
       xlab = "Difference in Elevation (m)", ylab = NA, cex.lab=2.5,
       cex.main = 2.5, cex.axis = 2, cex=1.5)
  sor.lmi=lm(bmi_12$beta.sim~dist.DIi_12)
  abline(sor.lmi, lwd=2)
  text(140, 0.86, "r = 0.55***(95% CI = 0.43 - 0.65)", cex=2)
  mtext('(f)', side=3, line=0.6, col = "black", at=last(-42), cex=1.5)
}


#-------------------------------------------------------------------------------
## 
#-------------------------------------------------------------------------------
# check distance effect on beta diversity
# use m_12, me_12 and mi_12 because the plots area still in the same order
distance_12 <- unlist(elevation_12[,5]) %>% as.vector()

#calculate beta-div
#install.packages("betapart")

# Total
bm_12_d <- beta.pair(m_12, index.family = "sorensen") 
dist.DI_12_d=dist(distance_12, method = "euclidean", diag = FALSE, upper = FALSE)

# edible
bme_12_d <- beta.pair(me_12, index.family = "sorensen") 
dist.DIe_12_d=dist(distance_12, method = "euclidean", diag = FALSE, upper = FALSE)

# inedible
bmi_12_d <- beta.pair(mi_12, index.family = "sorensen") 
dist.DIi_12_d=dist(distance_12, method = "euclidean", diag = FALSE, upper = FALSE)


##Mantel correlation test
library(ecodist)
# total
ecodist::mantel(bm_12_d$beta.sor ~ dist.DI_12_d, nperm = 1000, nboot = 500, pboot = 0.9, cboot = 0.95)
ecodist::mantel(bm_12_d$beta.sim ~ dist.DI_12_d, nperm = 1000, nboot = 500, pboot = 0.9, cboot = 0.95)

# edible
ecodist::mantel(bme_12_d$beta.sor ~ dist.DI_12_d, nperm = 1000, nboot = 500, pboot = 0.9, cboot = 0.95)
ecodist::mantel(bme_12_d$beta.sim ~ dist.DI_12_d, nperm = 1000, nboot = 500, pboot = 0.9, cboot = 0.95)

# inedible
ecodist::mantel(bmi_12_d$beta.sor ~ dist.DI_12_d, nperm = 1000, nboot = 500, pboot = 0.9, cboot = 0.95)
ecodist::mantel(bmi_12_d$beta.sim ~ dist.DI_12_d, nperm = 1000, nboot = 500, pboot = 0.9, cboot = 0.95)



{par(mfrow = c(2,3))
  par(mai = c(0.5, 0.7, 0.3, 0.17)) 
  
  windowsFont("Times New Roman")
  
  ###all forest BSor
  plot(bm_12_d$beta.sor~dist.DI_12_d, pch=16, main = "All species",
       ylim = c(0.3, 0.8), xlab = NA, ylab = "BSor", cex.lab=2,
       cex.main = 2.5, cex.axis = 2, cex=1.5, col="brown")
  sor.lm=lm(bm_12_d$beta.sor~dist.DI_12_d)
  abline(sor.lm, lwd=2)
  text(1100, 0.77, "r = 0.43**(95% CI = 0.32 - 0.52)", cex=1.7)
  mtext('(a)', side=3, line=0.6, col = "black", at=last(-150), cex=1.5)
  
  ###edible BSor
  plot(bme_12_d$beta.sor~dist.DIe_12_d, pch=19,  main = "Edible species",
       ylim = c(0.3, 0.9), xlab = NA, ylab = NA, cex.lab=2,
       cex.main = 2.5, cex.axis = 2, cex=1.5, col="dark green")
  sor.lme=lm(bme_12_d$beta.sor~dist.DIe_12_d)
  abline(sor.lme, lwd=2)
  text(1100, 0.86, "r = 0.11(95% CI = -0.02 - 0.31)", cex=1.7)
  mtext('(b)', side=3, line=0.6, col = "black", at=last(-150), cex=1.5)
  
  ###inedible BSOR
  plot(bmi_12_d$beta.sor~dist.DIi_12_d, pch=19,  main = "Inedible species",
       ylim = c(0.3, 0.9), xlab = NA, ylab = NA, cex.lab=2,
       cex.main = 2.5, cex.axis = 2, cex=1.5)
  sor.lmi=lm(bmi_12_d$beta.sor~dist.DIi_12_d)
  abline(sor.lmi, lwd=2)
  text(1100, 0.86, "r = 0.48***(95% CI = 0.39 - 0.55)", cex=1.7)
  mtext('(c)', side=3, line=0.6, col = "black", at=last(-150), cex=1.5)
  
  
  
  par(mai = c(0.72, 0.7, 0.1, 0.17)) 
  
  ### all forest BSim
  plot(bm_12_d$beta.sim~dist.DI_12_d, pch=16, 
       ylim = c(0.3, 0.8), xlab = "Difference in ground distance (m)", 
       ylab = "BSim", cex.lab=2,
       cex.main = 2.5, cex.axis = 2, cex=1.5, col="brown")
  sor.lm=lm(bm_12_d$beta.sim~dist.DI_12_d)
  abline(sor.lm, lwd=2)
  text(1100, 0.77, "r = 0.46**(95% CI = 0.36 - 0.57)", cex=1.7)
  mtext('(d)', side=3, line=1, col = "black", at=last(-150), cex=1.5)
  
  ###edible BSim
  plot(bme_12_d$beta.sim~dist.DIe_12_d, pch=19, ylim = c(0.3, 0.9), 
       xlab = "Difference in ground distance (m)", ylab = NA, cex.lab=2,
       cex.main = 2.5, cex.axis = 2, cex=1.5, col="dark green")
  sor.lme=lm(bme_12_d$beta.sim~dist.DIe_12_d)
  abline(sor.lme, lwd=2)
  text(1100, 0.86, "r = 0.06(95% CI = -0.09 - 0.27)", cex=1.7)
  mtext('(e)', side=3, line=1, col = "black", at=last(-150), cex=1.5)
  
  ###inedible BSim
  plot(bmi_12_d$beta.sim~dist.DIi_12_d, pch=19, ylim = c(0.3, 0.9),
       xlab = "Difference in ground distance (m)", ylab = NA, cex.lab=2,
       cex.main = 2.5, cex.axis = 2, cex=1.5)
  sor.lmi=lm(bmi_12_d$beta.sim~dist.DIi_12_d)
  abline(sor.lmi, lwd=2)
  text(1100, 0.86, "r = 0.50**(95% CI = 0.42 - 0.57)", cex=1.7)
  mtext('(f)', side=3, line=1, col = "black", at=last(-150), cex=1.5)
}

#-------------------------------------------------------------------------------
### 3.4	DIAMETER SIZE DISTRIBUTION
#-------------------------------------------------------------------------------

##################################################
## Computing dbh distribution curve for the forest

DBH_new <- read.csv("dbh_size_distribution.csv", header = TRUE, stringsAsFactors = FALSE)

## Assign plot diameter size abundance distribution for forest to a variable
dbh_densityplot =ggplot(DBH_new, aes(x=log(DBH), color=Species, group = Species),
                        cex=2, lwd=5, size=2) + 
  geom_density() +
  labs(color = "Category", y = "Relative tree density", x=expression('Log'['e']*'(DBH (cm))'))+
  theme(legend.position = c(.9, .8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(size=18,  family="Times New Roman"),
        axis.text = element_text(size=18, color="black",  family="Times New Roman"))

## create variable for plot names
Plots.labs <- c("Plot 1", "Plot 2", "Plot 3", "Plot 4", "Plot 5", "Plot 6", "Plot 7", "plot 8",
                "Plot 9", "Plot 10", "Plot 11", "Plot 12", "Plot 13", "Plot 14", "Plot 15")
names(Plots.labs) <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)

## assign plot of diameter size abundance distribution for each plot to variable
facet_plot <- ggplot(data=DBH_new, aes(x=log(DBH), color=Species, group=Species)) +
  geom_density(alpha = 0.6) +
  labs(y = "Relative tree density", x=expression('Log'['e']*'(DBH (cm))'), color= "Category")+
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(size=18,  family="Times New Roman"),
        axis.text = element_text(size=13, color="black",  family="Times New Roman")) +
  facet_wrap(~Plots, labeller = labeller(Plots = Plots.labs)) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x=element_blank()
  )

## merge diameter size abundnace distribution for forests an plot in one graph/plot
density_plot =ggarrange(dbh_densityplot, facet_plot, ncol = 2, nrow = 1 , labels = c("(a)","(b)"))

density_plot


#-------------------------------------------------------------------------------
### 3.5	TREE STAND AND BASAL AREA DENSITIES ALONG THE ELEVATIONAL BAND
#-------------------------------------------------------------------------------

##############################################################################
#### generating stand density (tree/ha) and mean basal area per plot for all,
### edible and inedible tree species
#################

#### plot density of trees per hectare and mean basal area per plot 
#### along the elevational gradie

col <- c("dark green", "brown")
palette.colors(palette = "Okabe-Ito")

BA_Den <- read.csv("basal area and plot density.csv", header = TRUE, stringsAsFactors = FALSE)

library(Hmisc)
plotdensity <- ggplot(BA_Den, aes(x=Elevation, y=Density, colour=Category))+ 
  stat_summary(fun.data= mean_cl_normal) + geom_point()+
  xlab("Elevation (m)") + ylab("Density (Trees per ha)")+
  geom_smooth(method = "lm", se=FALSE) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(size=18,  family="Times New Roman"))+
  annotate(geom="text", x=224, y=460, 
           label=expression("Edible: " *R^2*" = 0.02 ("*italic('p value')*" = 0.65)"),
           color="black", cex=5)+
  annotate(geom="text", x=224, y=440, 
           label=expression("Inedible: " *R^2*" = 0.34* ("*italic('p value')*" = 0.02)"),
           color="black", cex=5)


xl <- expression(Mean ~ Basal ~ Area ~ (m^2 ~ per ~ha))


plotBA<- ggplot(BA_Den, aes(x=Elevation, y=Total_BA, colour = Category))+ 
  stat_summary(fun.data= mean_cl_normal) +
  geom_point()+
  xlab("Elevation (m)") + ylab(xl)+
  geom_smooth(method = "lm", se=FALSE)+
  theme(legend.position = c(.9, .9),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(size=18,  family="Times New Roman"))+
  annotate(geom="text", x=224, y=42, 
           label=expression("Edible: " *R^2*" = 0.29* ("*italic('p value')*" = 0.04)"),
           color="black", cex=5)+
  annotate(geom="text", x=224, y=40, 
           label=expression("Inedible: " *R^2*" = 0.05 ("*italic('p value')*" = 0.43)"),
           color="black", cex=5)


library(ggpubr)
ppp <- ggarrange(plotdensity, plotBA, nrow = 1, labels = c("(a)","(b)"))
ppp


#### stand density computation
edible_d <- filter(BA_Den, Category=="Edible")
inedible_d <- filter(BA_Den, Category=="Inedible")

summary(lm(Density~Elevation, edible_d))
summary(lm(Density~Elevation, inedible_d))

#### basa area computation 
edible_BA <- filter(BA_Den, Category=="Edible")
inedible_BA <- filter(BA_Den, Category=="Inedible")

summary(lm(Total_BA~Elevation, edible_BA))
summary(lm(Total_BA~Elevation, inedible_BA))


###############################################
#### ANCOVA
###############################################
######
library(car)
# ancova for density
fit_d1 <- aov(Density~Elevation+Category, BA_Den)
Anova(fit_d1, type="III")

fit_d2 <- aov(Density~Elevation*Category, BA_Den)
Anova(fit_d2, type="III")

anova(fit_d1,fit_d2)

reg.todo_d <- lm(Density~Category/Elevation - 1, data=BA_Den)
anova(reg.todo_d)
summary(reg.todo_d)


# ancova for Basal Area
fit_BA1 <- aov(Total_BA~Elevation+Category, BA_Den)
Anova(fit_BA1, type="III")

fit_BA2 <- aov(Total_BA~Elevation*Category, BA_Den)
Anova(fit_BA2, type="III")

anova(fit_BA1,fit_BA2)

# nested
reg.todo_BA <- lm(Total_BA~Category/Elevation - 1, data=BA_Den)
anova(reg.todo_BA)
summary(reg.todo_BA)



