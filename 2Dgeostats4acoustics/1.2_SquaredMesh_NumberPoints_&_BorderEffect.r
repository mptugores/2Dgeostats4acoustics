library(gstat)
library(rgdal)

xmin <- 220000
ymin <- 4360000
range <- 140000

mygrid <- expand.grid(x=seq(xmin, xmin+range, by=1000), y=seq(ymin, ymin+range, by=1000))
mygrid$constant <- rep(1, nrow(mygrid))
coordinates (mygrid) <- ~x+y
proj4string (mygrid) <- CRS ("+proj=utm +zone=31 +ellps=WGS84")

myvariog <- variogram (constant~1, mygrid, cutoff=300000, width=1852)

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

par (family="f4", mfrow=c(1,1), mar=c(4,4,3,0.5), oma=c(0.5,0.5,0.5,0.5), cex=1.1, font.lab=1)

plot (round((myvariog$dist/1852), 0), myvariog$np, xlab="Distance (nm)", ylab="nº pairs", pch=20, col="darkgrey")
lines (round((myvariog$dist/1852), 0), myvariog$np)


# End of script!!!

