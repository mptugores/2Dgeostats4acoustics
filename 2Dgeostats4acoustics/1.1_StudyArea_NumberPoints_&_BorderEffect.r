
rm(list=ls())

library(gstat)
library(maptools)
library(rgdal)
library(lattice)
trellis.par.set(sp.theme())

# setting working directory
setwd("D:/PhD_tesis SEPTIEMBRE workingFiles/2_Trabajo_Tesis/AbundancePrecision/3_GEO_intrinseca_PhD/data_v1/gen_zones_utm") # you have to choose the unit you'll work in

source("../../functions/setproj.r")

# reading data
for (i in dir(pattern=".csv")){
assign (i, read.csv(i))
}

#-----------------------#
# Subset areas          #
#-----------------------#
# subset RB area
for (i in ls(pattern="gen")){
assign (paste("NS", substr(i,4,5), sep=""), subset(get(i), subzone=="RB"))
}

# subset DE area
for (i in ls(pattern="gen")){
assign (paste("SS", substr(i, 4,5), sep=""), subset(get(i), subzone=="DE"))
}

# setting projection for RB
for (j in ls(pattern="NS")){
assign(j, setproj (get(j), utm="31N"))
}

# setting projection for DE
for (j in ls(pattern="SS")){
assign(j, setproj (get(j), utm="31N"))
}

#-----------------------#
# Construct variograms  #
#-----------------------#
# variogram RB
w.ns <- c(1852, 4*1*1852, 4*2*1852)
for (k in w.ns){
  for (j in ls(pattern="NS03")){
  if (class(get(j))[1]=="SpatialPointsDataFrame"){
  assign (paste("vg.", j, "_", k, sep=""), variogram (m2ee~1, get(j), width=k, cutoff=500000))
  }
  }
}
rm (k, j)

# variogram DE
w.ss <- c(1852,8*1*1852, 8*2*1852)
for (k in w.ss){
  for (j in ls(pattern="SS03")){
  if (class(get(j))[1]=="SpatialPointsDataFrame"){
  assign (paste("vg.", j, "_", k, sep=""), variogram (m2ee~1, get(j), width=k, cutoff=500000))
  }
  }
}
rm (k, j)

#----------------------------#
# Plotting number of pairs of points
#----------------------------#
windowsFonts(
  f1=windowsFont("Arial Black"),
  f2=windowsFont("Bookman Old Style"),
  f3=windowsFont("Comic Sans MS"),
  f4=windowsFont("Times New Roman"),
  f5=windowsFont("Calibri"),
  f6=windowsFont("Trebuchet")
)

par (family="f4", mrow=c(2,3), mar=c(4,4,3,0.5), oma=c(0.5,0.5,0.5,0.5), cex=1.1, font.lab=1)
vg.order <- c("vg.NS03_1852", "vg.NS03_7408", "vg.NS03_14816","vg.SS03_1852","vg.SS03_14816","vg.SS03_29632")
for (j in vg.order){
  name.area <- substr (j,4,5)
  #plot (get(j)$dist, get(j)$np, main=name.area, xlab="distance(m)", ylab="nº pairs")
  plot (round((get(j)$dist/1852), 0), get(j)$np, main=name.area, xlab="Distance (nm)", ylab="nº pairs", pch=20, col="darkgrey")
  lines (round((get(j)$dist/1852), 0), get(j)$np)
}


#--------------
# Northern Subarea
#--------------
# Plots one by one  
par (family="f4", mar=c(4,4,3,0.5), oma=c(0.5,0.5,0.5,0.5), cex=1.1, font.lab=1)

j <- c("vg.NS03_1852")
name.area <- substr (j,4,5)
#plot (get(j)$dist, get(j)$np, main=name.area, xlab="distance(m)", ylab="nº pairs")
plot (round((get(j)$dist/1852), 0), get(j)$np, main=name.area, xlab="Distance (nm)", ylab="nº pairs", pch=20, col="grey20", ylim=c(0,5500))
lines (round((get(j)$dist/1852), 0), get(j)$np, col="grey20")

j <- c("vg.NS03_7408")
points (round((get(j)$dist/1852), 0), get(j)$np, pch=20, col="darkred")
lines (round((get(j)$dist/1852), 0), get(j)$np, col="darkred")

j <- c("vg.NS03_14816")
points (round((get(j)$dist/1852), 0), get(j)$np, pch=20, col="darkgreen")
lines (round((get(j)$dist/1852), 0), get(j)$np, col="darkgreen")

legend(95, 5480, c("EDSU", "ID", "2*ID"), col = c("grey20", "darkred", "darkgreen"),
       text.col = "black", lty = c(1, 1, 1), pch = c(20, 20, 20),
       merge = TRUE, bg = "white")

#--------------
# Southern Subarea
#--------------
# Plots one by one  
par (family="f4", mar=c(4,4,3,0.5), oma=c(0.5,0.5,0.5,0.5), cex=1.1, font.lab=1)

j <- c("vg.SS03_1852")
name.area <- substr (j,4,5)
#plot (get(j)$dist, get(j)$np, main=name.area, xlab="distance(m)", ylab="nº pairs")
plot (round((get(j)$dist/1852), 0), get(j)$np, main=name.area, xlab="Distance (nm)", ylab="nº pairs", pch=20, col="grey20", ylim=c(0,23000))
lines (round((get(j)$dist/1852), 0), get(j)$np, col="grey20")

j <- c("vg.SS03_14816")
points (round((get(j)$dist/1852), 0), get(j)$np, pch=20, col="darkred")
lines (round((get(j)$dist/1852), 0), get(j)$np, col="darkred")

j <- c("vg.SS03_29632")
points (round((get(j)$dist/1852), 0), get(j)$np, pch=20, col="darkgreen")
lines (round((get(j)$dist/1852), 0), get(j)$np, col="darkgreen")

legend(84, 23000, c("EDSU", "ID", "2*ID"), col = c("grey20", "darkred", "darkgreen"),
       text.col = "black", lty = c(1, 1, 1), pch = c(20, 20, 20),
       merge = TRUE, bg = "white")


#----------------------------#
# Border effect  - Kernel 2D #
#----------------------------#

# Density of points - Kernel 2D
# Script "KernelDensity_source_1.R" needs to be run to have the functions kde2d and kde2dplot available
library (MASS)
# script plot for kde2d
kde2dplot <- function(d,                # a 2d density computed by kde2D
                      ncol=num,          # the number of colors to use
                      zlim=c(0,max(z)), # limits in z coordinates
                      nlevels=20,       # see option nlevels in contour
		      theta=30,         # see option theta in persp
		      phi=30)           # see option phi in persp
		      {
z   <- d$z
nrz <- nrow(z)
ncz <- ncol(z)

couleurs  <- tail(topo.colors(trunc(1.4 * ncol)),ncol)
fcol      <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol      <- fcol[-nrz,-ncz]

par(mfrow=c(1,2),mar=c(0.5,0.5,0.5,0.5))
persp(d,col=fcol,zlim=zlim,theta=theta,phi=phi,zlab="density")

par(mar=c(2,2,2,2))
image(d,col=couleurs)
contour(d,add=T,nlevels=nlevels)
box()
}

#--------
# Northern Subarea
#--------

d <- kde2d(NS03$x31, NS03$y31, n=100)
kde2dplot(d)

#--------
# Southern Subarea
#--------
d <- kde2d(SS03$x31, SS03$y31, n=100)
kde2dplot(d)




################################################################################
# End of script!!!
################################################################################
