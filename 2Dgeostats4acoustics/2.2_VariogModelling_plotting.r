###########
# This script plots in a single .pdf file a variable amount of graphs
# obtained during the variogram fitting
# anchovy in the surroundings of Ebro river
# using different methodological options using gstat package
# 
# Created by: M.P. Tugores Ferra
###########

library(fields)
library(RGeoS)
library(gstat)
library(maptools)
library(rgdal)

######
### Data
######

# Anchovy 2003&2004: Delta de l'Ebre

#-------------------------------------------------------------------------------
# Delta de l'Ebre
#-------------------------------------------------------------------------------
rm (list=ls())
# dir.create ("../../../GeoStatistics_toPublish/2_variogram_modelling/differing_lagwidths/")
# wd <- "D:/ACOUSMED_contrato/GeoStatistics/GeoStatistics_toPublish/2_variogram_modelling/differing_lagwidths/"

wd <- "D:/PhD_tesis SEPTIEMBRE workingFiles/2_Trabajo_Tesis/AbundancePrecision/3_GEO_intrinseca_PhD/2_VariogModelling"
setwd (wd)

load("2_VariogModelling_20131009_154033_line415.RData")


#------------------------------------------------------------------------------#
#                         Variogram plots                                      #
#------------------------------------------------------------------------------#

#----------------------
# Nice variogram plots
#----------------------
windowsFonts(
  f1=windowsFont("Arial Black"),
  f2=windowsFont("Bookman Old Style"),
  f3=windowsFont("Comic Sans MS"),
  f4=windowsFont("Symbol"),
  f5=windowsFont("Calibri"),
  f6=windowsFont("Trebuchet")
)

# pdf ("variog.plots.pdf", onefile=TRUE)
# par (family="f5", mfrow=c(3,4), oma=c(0.1,0.1,0.1,0.1), mar=c(2,3,2,2))
par (mfrow=c(3,4), oma=c(0.1,0.1,0.1,0.1), mar=c(2,3,2,2))
for (i in 1:nrow(milista)){
  tmp.variog <- get (paste(milista$variogs [i]))
  tmp.fit <- get(gsub("pred.", "fit.", milista$preds) [i])
  title.plot <- substring(milista$preds[i], 9, nchar(as.vector(milista$preds[i])))
  plot (gamma~dist, tmp.variog, xlab="", ylab="",col="black", pch=19, fg="grey", axes=FALSE, main=title.plot)
  gamma.at <- seq(0, max(tmp.variog$gamma), 8000)
  gamma.labels <- gamma.at/1000
  axis(2,at=gamma.at, labels=gamma.labels, cex.axis=1.5, las=1)
  axis(1,at=seq(0,83340,18520), labels=seq(0,45,10), cex.axis=1.5, las=1)
  box(col="black")
  lines (variogramLine (tmp.fit, 220000, 200), col="grey")
}
# dev.off()






#------------------------------------------------------------------------------#
#---------------------------- End of script!!! --------------------------------#
#------------------------------------------------------------------------------#




