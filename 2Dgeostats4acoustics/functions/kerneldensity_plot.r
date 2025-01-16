# Density of points - Kernel 2D
# Script "KernelDensity_source_1.R" needs to be run to have the functions kde2d and kde2dplot available
library (MASS)
d <- kde2d(mydata$x,mydata$y,n=num)

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


kde2dplot(d)

# End of script
