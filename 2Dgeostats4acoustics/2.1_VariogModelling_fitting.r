###########
# This code is intended to be my 2nd paper in geostatistics and my 2nd chapter
# of my PhD thesis
#
# Computes and fit semivariograms for anchovy in the surroundings of Ebro river
# using different methodological options using gstat package
# 
# It also computes the CVgeo derived from the sampling design using RGeoS package
#
# Created by: M.P. Tugores Ferra
###########


library(fields)
library(RGeoS)
library(gstat)
library(maptools)
library (rgdal)

######
### Data
######

# Anchovy 2003&2004: Delta de l'Ebre

#-------------------------------------------------------------------------------
# Delta de l'Ebre
#-------------------------------------------------------------------------------
rm (list=ls())

# read data file and polygon representing the study area
#setwd ("D:/Handy_work/WKACUGEO/data_ecomeds/")
#setwd ("C:/Handy_copy/WKACUGEO/data_ecomeds/")

# wd <- "D:/ACOUSMED_contrato/GeoStatistics/data/data_utm31/"
# wd <- "D:/PhD_tesis SEPTIEMBRE workingFiles/2_Trabajo_Tesis/AbundancePrecision/3_GEO_intrinseca_PhD/data/data_utm31/"

setwd (wd)

# vector containing the names of the files accomplishing two different pattern criteria
v <- grep(".csv", list.files(pattern="geo_ecomed"), value=TRUE)


# Reading data in a loop
for (i in v){
  objname <- paste("ecomed", substring (i, 14, 15), sep="")
  assign (objname, read.csv(i, header=TRUE))
}


#--------------------------------------------------------------------
# Subsetting to the Study Area area (DE in this case)         
#--------------------------------------------------------------------
ecomed03 <- ecomed03[substr(ecomed03$radial,3,5)>33,]
ecomed04 <- ecomed04[substr(ecomed04$RADIAL,3,5)>33,]
ecomed05 <- ecomed05[substr(ecomed05$RADIAL,3,5)>33,]
ecomed06 <- ecomed06[substr(ecomed06$RADIAL,3,5)>33,]

# Creating a log transformed variable of sA.ancovy
ecomed03$ln.sA.anchovy <- log1p(ecomed03$m2ee)
ecomed04$ln.sA.anchovy <- log1p(ecomed04$m2.ee)
ecomed05$ln.sA.anchovy <- log1p(ecomed05$M2.boquerón)
ecomed06$ln.sA.anchovy <- log1p(ecomed06$M2_EE)

names (ecomed03) <- c("interval", "Total.sA", "date_m", "time_m", "lat", "lon", "transect", "sA.sardine", "sA.anchovy", "x", "y", "ln.sA.anchovy")
names (ecomed04) <- c("interval", "Total.sA", "sA.sardine", "sA.anchovy", "transect","lat", "lon", "x", "y", "ln.sA.anchovy")
names (ecomed05) <- c("interval", "Total.sA", "sA.sardine", "sA.anchovy", "transect", "date_m", "time_m", "lat", "lon", "x", "y", "ln.sA.anchovy")
names (ecomed06) <- c("interval", "Total.sA", "sA.sardine", "sA.anchovy", "lat", "lon", "transect", "x", "y", "ln.sA.anchovy")

#-------------------------
# Thresholded: R Objects thresholded 
#-------------------------
# rm (geo_ecomed_2006_utm)
for (i in ls(pattern="ecomed")){
assign (paste("th.", i, sep=""), get (i) [get(i)$sA.anchovy<526,])
}
# threshold 400 for 2004 & 2005
th.ecomed04 <- th.ecomed04[th.ecomed04$sA.anchovy<400,]
th.ecomed05 <- th.ecomed05[th.ecomed05$sA.anchovy<200,]

# Projecting data to plot in gstat
#-------------------------------------------------------------------------------
# Set coords & project
#-------------------------------------------------------------------------------
setproj <- function (mydata, utm) {
  library (rgdal)
  if (utm=="31N"){
  coordinates(mydata) <- c("x","y")
  proj4string(mydata) <- CRS("+proj=utm +zone=31 +ellps=WGS84")
  }
  if (utm=="30N"){
  coordinates(mydata) <- c("x30","y30")
  proj4string(mydata) <- CRS("+proj=utm +zone=30 +ellps=WGS84")
  }
  mydata
}

ecomed03 <- setproj (ecomed03, "31N")
ecomed04 <- setproj (ecomed04, "31N")
ecomed05 <- setproj (ecomed05, "31N")
ecomed06 <- setproj (ecomed06, "31N")

th.ecomed03 <- setproj (th.ecomed03, "31N")
th.ecomed04 <- setproj (th.ecomed04, "31N")
th.ecomed05 <- setproj (th.ecomed05, "31N")
th.ecomed06 <- setproj (th.ecomed06, "31N")
#-------------------------------------------------------------------------------

ecomed03 <- remove.duplicates (ecomed03)
th.ecomed03 <- remove.duplicates (th.ecomed03)

#----------------------------------
# Also ecomed2003 decimal
#----------------------------------
#setwd ("D:/ACOUSMED_contrato/GeoStatistics/ECOMED2003_general01nm/analisis_espacial/data/")
#ecomed03decimal <- read.table("generalmillasdecimal_utm.csv", sep=",",header=TRUE)

#ecomed03decimal$ln.sA.anchovy <- log1p(ecomed03decimal$ANE)
#names (ecomed03decimal)[4] <- "sA.anchovy"

# Set coords & project and remove duplicates
#ecomed03decimal <- setproj (ecomed03decimal, "31N")
#ecomed03decimal <- remove.duplicates (ecomed03decimal)
#-------------------------------------------------------------------------------

#-----------------------------------------------
# Variograms and plots
#-----------------------------------------------

#---------------------
# RAW and RAW.th dataset
#---------------------
w <- 1852  # lag with equal to EDSU
cutf <- 80000
years <- c("03","04","05","06")

for (y in years){
  for (j in ls (pattern=y)){ # perquè n'hi haurà un the la BD normal i l'altre, thresholded
    if (class(get(j))[1]=="SpatialPointsDataFrame"){
    assign (paste("vg.", j, "_", w, sep=""), variogram (sA.anchovy~1, get(j), width=w))
    }
  }
}


w <- 14816 # lag with equal to 2*inter-transect
cutf <- 80000
years <- c("03","04","05","06")

for (y in years){
  for (j in ls (pattern=y)){
    if (class(get(j))[1]=="SpatialPointsDataFrame"){
    assign (paste("vg.", j, "_", w, sep=""), variogram (sA.anchovy~1, get(j), width=w))
    }
  }
}

#---------------------
# LN
#---------------------
w <- 1852
# datasets <- c("ecomed03", "ecomed04", "ecomed05", "ecomed06", "ecomed03decimal")
datasets <- c("ecomed03", "ecomed04", "ecomed05", "ecomed06")
for (d in datasets){
  assign (paste("vg.", d, ".ln", "_", w, sep=""), variogram (ln.sA.anchovy~1, get(d), width=w))
}

w <- 14816
datasets <- c("ecomed03", "ecomed04", "ecomed05", "ecomed06")
for (d in datasets){
  assign (paste("vg.", d, ".ln", "_", w, sep=""), variogram (ln.sA.anchovy~1, get(d), width=w))
}



ls()
# [1] "cutf"                  "d"                     "datasets"
# [4] "ecomed03"              "ecomed03decimal"       "ecomed04"
# [7] "ecomed05"              "ecomed06"              "i"
#[10] "j"                     "objname"               "setproj"
#[13] "th.ecomed03"           "th.ecomed04"           "th.ecomed05"
#[16] "th.ecomed06"           "v"                     "vg.ecomed03"
#[19] "vg.ecomed03.ln"        "vg.ecomed03decimal"    "vg.ecomed03decimal.ln"
#[22] "vg.ecomed04"           "vg.ecomed04.ln"        "vg.ecomed05"
#[25] "vg.ecomed05.ln"        "vg.ecomed06"           "vg.ecomed06.ln"
#[28] "vg.th.ecomed03"        "vg.th.ecomed04"        "vg.th.ecomed05"
#[31] "vg.th.ecomed06"        "w"                     "wd"
#[34] "x"                     "y"                     "years"
#

# LN Backtransformed has to be done after fitting

#---------------------
# Fitting
#---------------------
# vgm.vector <- c("Sph", "Lin", "Exp") # variogram models
vgm.vector <- c("Sph") # variogram models
fm.vector <- c(7) # fitting methods

#---------------------
for (i in ls(pattern="vg.")){
  for (vgm in vgm.vector){
    for (fm in fm.vector){
      if (class(get(i))[1] == "gstatVariogram"){
      f.nug <- 0
      f.sill <- mean (get(i)[,3])
      f.range <- 2/3*max(get(i)[,2])
      fit.name <- paste("fit.",i, ".", vgm, ".", fm, sep="")
      print(fit.name)
      fitting <- fit.variogram (get(i), vgm (f.sill, vgm, f.range, f.nug), fit.method=fm)
      #warn <- warnings (fitting)
      assign (fit.name, fitting)
  #capture.output (fit.name, warn, "output.fitting.variograms.txt")
      }
    }
  }
}
remove (fitting)

#-------------------------------------------------------------------------------
# Saving warnings in a .csv
#-------------------------------------------------------------------------------
loop.text <- capture.output (
  for (i in ls(pattern="vg.")){
    for (vgm in vgm.vector){
      for (fm in fm.vector){
        if (class(get(i))[1] == "gstatVariogram"){
        f.nug <- 0
        f.sill <- mean (get(i)[,3])
        f.range <- 2/3*max(get(i)[,2])
        fit.name <- paste("fit.",i, ".", vgm, ".", fm, sep="")
        print(fit.name)
        fitting <- fit.variogram (get(i), vgm (f.sill, vgm, f.range, f.nug), fit.method=fm)
        assign (fit.name, fitting)
        }
      }
    }
  }
)

remove (fitting)

# Working directory para guardar datos
# wd <- "D:/ACOUSMED_contrato/GeoStatistics/GeoStatistics_toPublish/2_variogram_modelling/differing_lagwidths/"
#wd <- "C:/Users/user/Documents/geostatistics/"

wd <- "D:/PhD_tesis SEPTIEMBRE workingFiles/2_Trabajo_Tesis/AbundancePrecision/3_GEO_intrinseca_PhD/2_VariogModelling"
setwd (wd)

# write.csv (loop.text, "warnings.fitting.csv")


#---------------------
# ln backtransformed
#---------------------
#--------------------------------------
# Backtransforming the fitted values - ARREGLADA
#--------------------------------------
fits.ln.vector <- agrep("fit.", ls(pattern=".ln"), value=TRUE) #mi set de datos a backtransformar
fits.lnbt.vector <- gsub ("ln", "lnbt", fits.ln.vector)  #copio los .ln a .lnbt; el contenido es de momento el mismo que
# en los ln y sólo cambia el nombre

# Cambia el contenido de los fits lnbt (que eran los mismos que los log) x los lnbt reales
for (i in 1:length(fits.ln.vector)){
  j <- fits.ln.vector [i]
  k <- fits.lnbt.vector [i]
  if (class(get(j))[1]=="variogramModel"){
    assign (k, get(j))
    gk <- get (k)
    dataset <-  strsplit (strsplit(j, ".ln")[[1]], "vg.")[[1]][2]
    variable.id <- as.numeric(which (names (slot(get(dataset), "data"))=="sA.anchovy"))
    variable.ln.id <- as.numeric(which (names (slot(get(dataset), "data"))=="ln.sA.anchovy"))
    variable <- as.data.frame (get(dataset)[,variable.id])[,1]
    variable.ln <- as.data.frame (get(dataset)[,variable.ln.id])[,1]
    A <- (mean (variable)^2)+var(variable)
    sigma.sqr <- log (1+(var(variable)/(mean(variable)^2)))
    B.psill <- 1 - exp (-(sigma.sqr*get(j)[2,2])/var(variable.ln))
    B.pnugget <- 1 - exp (-(sigma.sqr*get(j) [1,2])/var(variable.ln))
    semivar.h.psill <- A*B.psill
    semivar.h.pnugget <-A*B.pnugget
    gk [2,2] <- semivar.h.psill
    gk [1,2] <- semivar.h.pnugget
    assign (k, gk)
  }
}

rm (fits.ln.vector, fits.lnbt.vector, i, j, k, gk, dataset, variable.id, variable.ln.id,
variable, variable.ln, A, sigma.sqr, B.psill, B.pnugget, semivar.h.psill, semivar.h.pnugget)

#----------------------------------
# Retrieving fitted parameters
#----------------------------------
id <- 0
x <- ls(pattern="fit.vg")
for (i in x){
  v <- c(i, get(i)[2,3], get(i)[2,2], get (i)[1,2])
  if (id==0){
    new.v <- v
  }
  if (id>0) {
    new.v <- rbind(new.v, v)
  }
id <- id+1
}

fitted.param <- data.frame (model=new.v[,1], range=as.numeric(new.v[,2]), 
psill=as.numeric(new.v[,3]), nugget=as.numeric(new.v[,4]))
#fitted.param <- fitted.param[-which(fitted.param$model == "fit.vg.ecomed03decimal_1852.Exp.7"),]

# Falta "pred.vg.ecomed03decimal.Exp.7" xq el rango es negativo
rm (id, x, i, v, new.v)

# dir.create ("../../../GeoStatistics_toPublish/2_variogram_modelling/differing_lagwidths/")
# wd <- "D:/ACOUSMED_contrato/GeoStatistics/GeoStatistics_toPublish/2_variogram_modelling/differing_lagwidths/"

wd <- "D:/PhD_tesis SEPTIEMBRE workingFiles/2_Trabajo_Tesis/AbundancePrecision/3_GEO_intrinseca_PhD/2_VariogModelling"
setwd (wd)

# write.csv(fitted.param, "fitted.parameters.csv")

#--------------------------------------
# Loop for variogramLine predictions  
#--------------------------------------

id <- 0
for (f in ls(pattern="fit.")){
  if (class(get(f))[1]=="variogramModel"){
  id <- id+1
  print (id)
    if (get(f)[2,3]>0){
      print (f)
      k <- strsplit (f, "fit.")[[1]][2]
      name.line <- paste ("pred.",k, sep="")
      if (length(grep(".Sph.", k))>0){
      m <- strsplit(k, ".Sph.")[[1]][1]
      }
      if (length(grep(".Exp.", k))>0){
      m <- strsplit(k, ".Exp.")[[1]] [1]
      }
      if (length(grep("lnbt", m))>0){
      m <- gsub("lnbt", "ln", m)
      }
    v.dist <- get(m)[,2]
    prediction <- variogramLine (get(f), dist_vector=v.dist)
    print (name.line)
    print (prediction)
    assign(name.line, prediction)
    preds.variogs <- c (name.line, m)  # this lines produces a table with preds and variogs
    if (id==1){
      new.preds.variogs <- preds.variogs
    }
    if (id>1) {
      new.preds.variogs <- rbind(new.preds.variogs, preds.variogs)
    }  
  }
  if (get(f)[2,3]<0){
  print(f)
  print(get(f)[2,3])
  }
}
}

# Tabla de equivalencias
milista <- data.frame (preds=new.preds.variogs[,1], variogs=new.preds.variogs[,2])
rm (new.preds.variogs)

#----------------------------------------
# Loop for computing fitting parameters 
#----------------------------------------
for (i in 1:length(milista[,1]))  {
  # retrieving necessary information
  pred.i <- get (paste(milista$preds[i])) [,2]
  vg.i <- get (paste(milista$variogs[i])) [,3]
  npairs.i <- get (paste(milista$variogs[1])) [,1]
  dist.i <- get (paste(milista$variogs[1])) [,2]
  # fitting indices
  bias <- sum (pred.i - vg.i ) # bias
  sse <- sum ((pred.i - vg.i )^2) # sum of squared errors
  abs.bias <- sum(abs(pred.i  - vg.i )) # absolute bias
  gof <- sse/sum(vg.i^2) # goodness of fit
  wgof2 <- sum (npairs.i * (pred.i - vg.i )^2)/sum (npairs.i^2)# weighted goodness of fit 2
  wgof3 <- sum ((1/dist.i)* (pred.i - vg.i )^2)/sum ((1/dist.i^2))# wgof 3
  wgof4 <- sum (npairs.i *(1/dist.i)* (pred.i - vg.i )^2)/sum (npairs.i^2*(1/dist.i^2))# wgof 4
  lm.x <- lm (vg.i ~ pred.i)
  r2 <- 1-(sum(residuals(lm.x)^2)/sum((vg.i-mean(vg.i))^2)) # r^2
  # gathering all indices in a vector
  fit.vector <- c (bias, sse, abs.bias, r2, gof, wgof2, wgof3, wgof4)
  print (i)
  print (fit.vector)
  if (i==1){
  new.fit.vector <- fit.vector
  }
  if (i>1){
  new.fit.vector <- as.data.frame(rbind (new.fit.vector, fit.vector))
  }
  #print (warnings())
}

names (new.fit.vector) <- c("bias", "sse", "abs.bias", "r2", "gof", "wgof2", "wgof3", "wgof4")

results.fitting  <- cbind (milista,fitted.param, new.fit.vector)

# Saving fittings in a .csv file
# wd <- "D:/ACOUSMED_contrato/GeoStatistics/GeoStatistics_toPublish/2_variogram_modelling/differing_lagwidths/"
wd <- "D:/PhD_tesis SEPTIEMBRE workingFiles/2_Trabajo_Tesis/AbundancePrecision/3_GEO_intrinseca_PhD/2_VariogModelling"
setwd (wd)

# write.csv (results.fitting , "goodness.of.fit.csv", row.names=FALSE)


# Save image
date.time <- paste(substring(gsub(":", "", gsub("-", "", Sys.time())), 1, 8), 
substring(gsub(":", "", gsub("-", "", Sys.time())), 10, 15), sep="_")

# save.image (paste("2_VariogModelling", date.time, "_line415.RData", sep=""))

rm (date.time)

#------------------------------------------------------------------------------#
#                         STOPPED HERE!!! 16/04/2012                           #
#------------------------------------------------------------------------------#

# Ahora hay que pasar a la versión 32 bits porque RGeoS no está disponible en 
# versión 64 bits

#------------------------------------------------------------------------------#
#                          CV geostatistic                                     #
#------------------------------------------------------------------------------#

# Abrimos R 32 bits
load("C:\\Users\\user\\Documents\\geostatistics\\1_variogmodelling_ndisc250_bt.arreglada20120416_180924_line467.RData")


# Importing a .shp polygon file 
# and getting vertix coordinates
library (maptools)
library (RGeoS)
# for old version of maptools
# a.DE <- read.shape ("pols/sampled_area_DEutm31N.shp")
# polDE <- a.DE[[1]][[1]]$verts
a.DE <- readShapeSpatial ("geostatistics/pols/sampled_area_DEutm31N_0308.shp")
polDE <- a.DE@polygons [[1]]@Polygons [[1]]@coords
db.polDE <- polygon.create(x=polDE[,1],y=polDE[,2],polygon=NA)


# Creating a RGeoS database
db.DE03 <- db.create(ecomed03,ndim=2,nvar=1,nech=dim(ecomed03)[1],autoname=TRUE)
# anchovy
# define locators
db.DE03 <- db.locate(db.DE03,11:12,loctype="x")
db.DE03 <- db.locate(db.DE03,10,loctype="z")
# select points inside polygon DE
db.DE03 <- db.polygon(db.DE03,db.polDE)
cat("nb points: ",db.DE03$nech," ; outside polygon: ",sum(!db.DE03@items$sel),"\n")

# Variogram
projec.toggle(mode=0)
Lag <- 8000; Nlag<-10
vg.eeDE03 <- vario.calc(db.DE03, dirvect=0,lag=Lag, tolang=90, nlag=Nlag, breaks = NA,
           calcul="vg", by.sample=FALSE, opt.code=0, tolcode=0, means=NA)
vg.eeDE03
plot(vg.eeDE03,xlab="Distance",ylab="variogram")

nv <- 1
vg1.init03 <- model.input(model.default=NA,ndim = 2, nvar = nv, order=0,
	  flag.sill=TRUE, flag.plot=FALSE)
2      # count of basic structures
1      # nugget
5691.298   # nugget value
3      # 2:exponentiam, 3:spherical
2176.652  # sill
n      # anisotropy
58709.83  #range

vg1.fit03 <- model.fit(vario=vg.eeDE03, model=vg1.init03, niter = 100, wmode = 2, draw =T,
	flag.ask = FALSE, flag.fit = TRUE)
vg1.fit03

#Covariance Part
#---------------
#- Nugget Effect
#  Sill        =  11726.215
#- Spherical
#  Range       =  58709.830
#  Sill        =      0.000
#Total Sill    =  11726.215

# Datasets
#mydata <- c("ecomed03", "ecomed04", "ecomed05", "ecomed06", "ecomed03decimal", 
#"th.ecomed03", "th.ecomed04", "th.ecomed05", "th.ecomed06", "ecomed03.ln", "ecomed04.ln",
#"ecomed05.ln", "ecomed06.ln", "ecomed03decimal.ln")


# sin decimal
mydata <- c("ecomed03", "ecomed04", "ecomed05", "ecomed06", 
"th.ecomed03", "th.ecomed04", "th.ecomed05", "th.ecomed06", "ecomed03.ln", "ecomed04.ln",
"ecomed05.ln", "ecomed06.ln")

# mydata.ln <- c(ecomed03.ln, ecomed04.ln, ecomed05.ln, ecomed06.ln, ecomed03decimal.ln)
# mydata2 <- c("ecomed03", "ecomed04")

#-------------------------------------------------------------------------------
# Loop for creating RGeoS databases
#-------------------------------------------------------------------------------
for (i in mydata) {
# create RGeoS data structure
  print (i)
  if (length(grep (".ln", i))!=0){
  data.tmp <- as.data.frame(get(paste(gsub(".ln", "", i))))
  db.tmp <- db.create(data.tmp, ndim=2, nvar=1, nech=dim(data.tmp)[1], autoname=TRUE)
  # define locators
  names.coord.x <- which(names(data.tmp)=="x")+1
  names.coord.y <- which(names(data.tmp)=="y")+1
  db.tmp <- db.locate(db.tmp,names.coord.x:names.coord.y,loctype="x")
  db.tmp <- db.locate(db.tmp,which(names(data.tmp)=="ln.sA.anchovy")+1,loctype="z")
  }
  if (length(grep(".ln", i))==0) {
  data.tmp <- as.data.frame(get(paste(i)))
  db.tmp <- db.create(data.tmp, ndim=2, nvar=1, nech=dim(data.tmp)[1], autoname=TRUE)
  # define locators
  names.coord.x <- which(names(data.tmp)=="x")+1
  names.coord.y <- which(names(data.tmp)=="y")+1
  db.tmp <- db.locate(db.tmp,names.coord.x:names.coord.y,loctype="x")
  db.tmp <- db.locate(db.tmp,which(names(data.tmp)=="sA.anchovy")+1,loctype="z")
  }
  #print (db.tmp)
  # sets the limits of study area : limits of plots
  #if (length(which(names(data.tmp)=="ln.sA.anchovy"))!=0){
  #print ("Existe la columna ln.sA.anchovy. ¿Qué quieres que haga con esto?"
  #}
  x1 <- min(data.tmp[,which(names(data.tmp)=="x")])
  x2 <- max(data.tmp[,which(names(data.tmp)=="x")]) 
  y1 <- min(data.tmp[,which(names(data.tmp)=="y")])
  y2 <- max(data.tmp[,which(names(data.tmp)=="y")])
  # select points inside polygon DE
  db.tmp <- db.polygon(db.tmp,db.polDE)
  cat("nb points: ",db.tmp$nech," ; outside polygon: ",sum(!db.tmp@items$sel),"\n")  
  assign (paste("db.", i, sep=""), db.tmp) 
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Creando modelos de fittings para RGeoS con los fits de gstat
#-------------------------------------------------------------------------------
rm(milista)

for (k in results.fitting$preds){
    m.name <- gsub ("pred", "fit", k)
    RGeoS.mname <- paste ("RGeoS.", m.name, sep="")
    assign (RGeoS.mname , vg1.fit03)
    RGeoS.obj <- get(RGeoS.mname)
    # slot(slot(RGeoS.obj, "basics")[[1]], "sill") <- as.matrix(11111111)
    # slot(slot(RGeoS.obj, "basics")[[2]], "range") <- 333333
    # slot(slot(RGeoS.obj, "basics")[[2]], "sill") <- as.matrix(999999999)
    tmp.nugg <- get (m.name) [1,2] # nugget
    tmp.range <- get (m.name) [2,3] # range
    print (get(m.name))
    tmp.sill <- get (m.name) [2,2] # sill
    slot(slot(RGeoS.obj, "basics")[[1]], "sill") <- as.matrix(tmp.nugg)
    slot(slot(RGeoS.obj, "basics")[[2]], "range") <- tmp.range
    slot(slot(RGeoS.obj, "basics")[[2]], "sill") <- as.matrix(tmp.sill)
    if (length(grep ("Exp", m.name))>0){
    slot(slot(RGeoS.obj, "basics")[[2]], "vartype") <- "Exponential"
    }
    if (length(grep ("Sph", m.name))>0){
    slot(slot(RGeoS.obj, "basics")[[2]], "vartype") <- "Spherical"
    }
    if (length(grep ("Lin", m.name))>0){
    slot(slot(RGeoS.obj, "basics")[[2]], "vartype") <- "Linear"
    }
    assign (RGeoS.mname, RGeoS.obj)
}
#-------------------------------------------------------------------------------

#----------------------------------------------------------
# Añadiendo en milista una columna con equivalencias para 
# las bases de datos de RGeoS y para los fittings de RGeoS
#----------------------------------------------------------
results.fitting$db <- gsub ("vg.", "db.", results.fitting$variogs) # variable intermedia; no tener en cuenta
results.fitting$db.def <- grep("db.", unlist(strsplit(results.fitting$db, "_")), value=TRUE)
# corrigiendo bases de datos
results.fitting$db.def [grep("lnbt", results.fitting$model)] <- grep("db.", unlist(strsplit(gsub(".ln_", "_", results.fitting$db[grep("lnbt", results.fitting$model)]), "_")), value=TRUE)
results.fitting$RGeoS.fits <- gsub("pred.", "RGeoS.fit.", results.fitting$preds)

# Como ecomed03decimal ralentiza mucho el ordenador, lo excluyo temporalmente de mis
# cálculos
# milista.notdecimal <- results.fitting [-c(grep("decimal", results.fitting$RGeoS.fits)),]
date.time <- paste(substring(gsub(":", "", gsub("-", "", Sys.time())), 1, 8), 
substring(gsub(":", "", gsub("-", "", Sys.time())), 10, 15), sep="_")

wd <- "C:/Users/user/Documents/geostatistics/"
setwd (wd)

#save.image (paste("1_variogmodelling_ndisc250_bt.arreglada", date.time, "_line648.RData", sep=""))
rm (date.time)

#wd <- "D:/ACOUSMED_contrato/GeoStatistics/GeoStatistics_toPublish/2_variogram_modelling/differing_lagwidths/"


#------------------------
# Loop computing CVgeo
#------------------------
# rm (list=ls (pattern="fit.vg")[-grep("RGeoS", ls(pattern="fit"))])
# rm (list=ls(pattern="pred.vg"))
# rm (wgof2,wgof3,wgof4,x1,x2,y,y1,y2,years)

# results.fitting1 <- results.fitting[1:100,]
# results.fitting2 <- results.fitting[101:200,]
# results.fitting3 <- results.fitting[201:300,]
# results.fitting4 <- results.fitting[301:400,]
# results.fitting5 <- results.fitting[401:500,]
# results.fitting6 <- results.fitting[501:600,]
# results.fitting7 <- results.fitting[601:700,]
# results.fitting8 <- results.fitting[701:759,]

# memory.limit(3000)
# memory.limit (2800)

results.fitting$year <- rep(1, nrow(results.fitting))
results.fitting$model.type <- rep(1, nrow(results.fitting))
results.fitting$transf <- rep(1, nrow(results.fitting))
results.fitting$fit.method <- rep (1, nrow(results.fitting))

results.fitting$year[grep ("ecomed03", results.fitting$model, fixed=TRUE)] <- "2003"
results.fitting$year[grep ("ecomed04", results.fitting$model, fixed=TRUE)] <- "2004"
results.fitting$year[grep ("ecomed05", results.fitting$model, fixed=TRUE)] <- "2005"
results.fitting$year[grep ("ecomed06", results.fitting$model, fixed=TRUE)] <- "2006"
# results.fitting$year[grep("ecomed03decimal", results.fitting$model, fixed=TRUE)] <- "2003d"

results.fitting$model.type [grep("Exp.", results.fitting$model, fixed=TRUE)] <- "Exp"
results.fitting$model.type [grep("Sph.", results.fitting$model, fixed=TRUE)] <- "Sph"

results.fitting$transf [grep(".th.", results.fitting$model, fixed=TRUE)] <- "th"
results.fitting$transf [grep(".ln_", results.fitting$model, fixed=TRUE)] <- "ln"
results.fitting$transf [grep(".lnbt_", results.fitting$model, fixed=TRUE)] <- "lnbt"
results.fitting$transf [results.fitting$transf=="1"] <- "none"

results.fitting$fit.method [grep(".1", results.fitting$model, fixed=TRUE)] <- "1"
results.fitting$fit.method [grep(".2", results.fitting$model, fixed=TRUE)] <- "2"
results.fitting$fit.method [grep(".6", results.fitting$model, fixed=TRUE)] <- "6"
results.fitting$fit.method [grep(".7", results.fitting$model, fixed=TRUE)] <- "7"

date.time <- paste(substring(gsub(":", "", gsub("-", "", Sys.time())), 1, 8), 
substring(gsub(":", "", gsub("-", "", Sys.time())), 10, 15), sep="_")

#save.image (paste("1_variogmodelling_ndisc250_bt.arreglada_", date.time, "_line703_forCVestimation.RData", sep=""))
rm (date.time)

# write.csv (results.fitting , "results.fitting.csv", row.names=FALSE)

#-------------------------------------------------------------------------------
# Todo junto                                                                   
#-------------------------------------------------------------------------------
for (i in 1:length(results.fitting$RGeoS.fits)) {
  # print (i)
  # print (results.fitting$RGeoS.fits[i])
  # CVgeo by kriged mean (RGeoS)
  tmp <- global(get(results.fitting$db.def[i]), model = get(results.fitting$RGeoS.fits[i]), uc=c("1"), 
  polygon = db.polDE , calcul = "krige", ndisc=250, flag.polin=TRUE, flag.wgt=TRUE, ivar = 1, verbose = 1, )
  assign (paste("CVgeo.", results.fitting$RGeoS.fits[i], sep=""), tmp)
  tmp.fit <- results.fitting$RGeoS.fits[i]
  CV <- tmp$cv
  if (i==1){
  CV.new <- CV
  }
  if (i>1){
  CV.new <- c(CV.new, CV)
  }
} 

date.time <- paste(substring(gsub(":", "", gsub("-", "", Sys.time())), 1, 8), 
substring(gsub(":", "", gsub("-", "", Sys.time())), 10, 15), sep="_")

save.image (paste("1_variogmodelling_ndisc250_bt.arreglada_", date.time, "_line731_CVestimation.RData", sep=""))

rm (CV, tmp, tmp.fit)
rm (date.time)
#-------------------------------------------------------------------------------

wd <- "D:/ACOUSMED_contrato/GeoStatistics/GeoStatistics_toPublish/2_variogram_modelling/differing_lagwidths/"
setwd (wd)
# write.csv (CV.new, "CVgeo_ndisc250.csv")



#-------------------------------------------------------------------------------
# Result Outputs
#-------------------------------------------------------------------------------

# my.results1$range <- as.numeric(levels(my.results1$range))[my.results1$range]
# my.results1$psill <- as.numeric(levels(my.results1$psill)) [my.results1$psill]
# my.results1$nugget <- as.numeric(levels(my.results1$nugget)) [my.results1$nugget]

results.variogs <- cbind(results.fitting, CV.new)

# Range in nm and practical range
results.variogs$range.nm <- round (results.variogs$range/1852, 0)
results.variogs$pract.range <- results.variogs$range.nm
results.variogs$pract.range [grep("Exp", results.variogs$model, fixed=TRUE)] <- (results.variogs$pract.range[grep("Exp", results.variogs$model, fixed=TRUE)]*3)

# Modificado 23/09/2011
# results.variogs$filt <- rep ("1", nrow (results.variogs))
# results.variogs$filt [results.variogs$range.nm >60] <- "bad"
# results$filt [results$range.nm < 1] <- "bad"
# results$filt [results$range.nm < 2 & results$nugget==0] <- "bad"
# results.variogs$filt [results.variogs$range.nm < 2] <- "bad"
# results.variogs$filt [results.variogs$range.nm < 3 & results.variogs$nugget==0] <- "bad"

# Filter on the practical range
results.variogs$filt <- rep ("1", nrow (results.variogs))
results.variogs$filt [results.variogs$pract.range >60] <- "bad"
results.variogs$filt [results.variogs$pract.range < 2] <- "bad"
results.variogs$filt [results.variogs$pract.range < 3 & results.variogs$nugget==0] <- "bad"

results.variogs$filt [results.variogs$psill==0] <- "bad"
results.variogs$filt [results.variogs$model=="fit.vg.ecomed05.Exp.1"] <- "bad"
results.variogs$filt [results.variogs$filt=="1"] <- "good"

# write.csv(results.variogs, "results.variogs.csv") 
# write.csv(results, "../../../GeoStatistics_toPublish/2_variogram_modelling/results.variogs.csv") 

#------------#
# Filtering  #
#------------#
results.variogs.filt <- results.variogs [results.variogs$filt=="good",]
results.variogs.filt$total.sill <- (results.variogs.filt$psill)+(results.variogs.filt$nugget)
results.variogs.filt$SpD <- (results.variogs.filt$psill)/(results.variogs.filt$total.sill)

# write.csv(results.variogs.filt, "results.variogs.filt.csv") 

# ver a ver si funcionó bien!!!

#----------------------------------
# Aggregate per transf i year
#----------------------------------

#####################
# Function "count"
count <- function(x) {
length(na.omit(x))
} 
####################

names(results.variogs.filt)
[1] "preds"       "variogs"     "model"       "range"       "psill"      
[6] "nugget"      "bias"        "sse"         "abs.bias"    "r2"         
[11] "gof"         "wgof2"       "wgof3"       "wgof4"       "db"         
[16] "db.def"      "RGeoS.fits"  "year"        "model.type"  "transf"     
[21] "fit.method"  "CV.new"      "range.nm"    "pract.range" "filt"       
[26] "total.sill"  "SpD"       
names(results.variogs.filt) [22] <- "CVgeo"

output.mean <- aggregate(cbind(pract.range, CVgeo, psill, nugget) ~ transf + year, data = results.variogs.filt, mean)
output.sd <- aggregate(cbind(pract.range, CVgeo, psill, nugget) ~ transf + year, data = results.variogs.filt, sd)
output.count <- aggregate(cbind(model)~ transf + year, data = results.variogs.filt, count)

names (output.mean) [3:6] <- c("mean.pract.range.nm", "mean.CVgeo", "mean.psill", "mean.nugget")
names (output.sd) [3:6] <- c("sd.pract.range.nm", "sd.CVgeo", "sd.psill", "sd.nugget")
names (output.count) [3] <- c("count")

output <- cbind (output.mean, output.sd[,3:6], output.count[,3])

output$cv.range <- round(output$sd.pract.range.nm/output$mean.pract.range.nm, 4)
output$cv.CVgeo <- round(output$sd.CVgeo/output$mean.CVgeo, 4)
output$cv.psill <- round(output$sd.psill/output$mean.psill, 4)
output$cv.nugget <- round(output$sd.nugget/output$mean.nugget, 4)

output$mean.range.nm <- round(output$mean.pract.range.nm, 4)
output$mean.CVgeo <- round(output$mean.CVgeo, 4)
output$sd.range.nm <- round(output$sd.pract.range.nm, 4)
output$sd.CVgeo <- round(output$sd.CVgeo, 4)
output$mean.pract.range.nm <- round(output$mean.pract.range.nm, 4)
output$mean.CVgeo <- round(output$mean.CVgeo, 4)
output$sd.pract.range.nm <- round(output$sd.pract.range.nm, 4)
output$sd.CVgeo <- round(output$sd.CVgeo, 4)
output$sd.psill <- round(output$sd.psill, 4) 
output$sd.nugget <- round(output$sd.nugget, 4) 
output$mean.psill <- round(output$mean.psill, 4) 
output$mean.nugget <- round(output$mean.nugget, 4) 

# write.csv (output, "aggregate_CV_range.CVgeo.fitted.param_ndisc250.csv", row.names=TRUE)
# write.csv (output, "D:/ACOUSMED_contrato/GeoStatistics/GeoStatistics_toPublish/2_variogram_modelling/CV_range.CVgeo.fitted.param_ndisc250.csv", col.names=TRUE, row.names=TRUE)

#---------------------------------------
# Aggregate per transf, year i model2
#---------------------------------------

output.mean2 <- aggregate(cbind(pract.range, CVgeo, psill, nugget) ~ transf + year + model.type, data = results.variogs.filt, mean)
output.sd2 <- aggregate(cbind(pract.range, CVgeo, psill, nugget) ~ transf + year + model.type, data = results.variogs.filt, sd)
output.count2 <- aggregate(cbind(model)~ transf + year + model.type, data = results.variogs.filt, count)

names (output.mean2) [4:7] <- c("mean.pract.range.nm", "mean.CVgeo","mean.psill", "mean.nugget")
names (output.sd2) [4:7] <- c("sd.pract.range.nm", "sd.CVgeo", "sd.psill", "sd.nugget")
names (output.count2) [4] <- c("count")

output2 <- cbind (output.mean2, output.sd2[,4:7], output.count2[,4])

output2$cv.range <- round(output2$sd.pract.range.nm/output2$mean.pract.range.nm, 4)
output2$cv.CVgeo <- round(output2$sd.CVgeo/output2$mean.CVgeo, 4)

output2$mean.range.nm <- round(output2$mean.pract.range.nm, 4)
output2$mean.CVgeo <- round(output2$mean.CVgeo, 4)
output2$sd.range.nm <- round(output2$sd.pract.range.nm, 4)
output2$sd.CVgeo <- round(output2$sd.CVgeo, 4)
output2$mean.pract.range.nm <- round(output2$mean.pract.range.nm, 4)
output2$mean.CVgeo <- round(output2$mean.CVgeo, 4)
output2$sd.pract.range.nm <- round(output2$sd.pract.range.nm, 4)
output2$sd.CVgeo <- round(output2$sd.CVgeo, 4)

output2$sd.psill <- round(output2$sd.psill, 4) 
output2$sd.nugget <- round(output2$sd.nugget, 4) 
output2$mean.psill <- round(output2$mean.psill, 4) 
output2$mean.nugget <- round(output2$mean.nugget, 4) 

# write.csv (output2, "aggregate2_CV_range.CVgeo.fitted.param_ndisc250_bytransf.year.model.csv", row.names=TRUE)
# write.csv (output2, "D:/ACOUSMED_contrato/GeoStatistics/GeoStatistics_toPublish/2_variogram_modelling/CV_range.CVgeo.fitted.param_ndisc250_bytransf.year.model.csv", col.names=TRUE, row.names=TRUE)

# Sólo por model type no tiene sentido.

#  model.type mean.pract.range.nm mean.CVgeo sd.pract.range.nm sd.CVgeo cv.range cv.CVgeo
# 1        Exp               35.62      16.95             16.28    10.54    45.70    62.16
# 2        Sph               24.00      12.96              8.82     7.81    36.73    60.22

date.time <- paste(substring(gsub(":", "", gsub("-", "", Sys.time())), 1, 8), 
substring(gsub(":", "", gsub("-", "", Sys.time())), 10, 15), sep="_")

save.image (paste("1_variogmodelling_ndisc250_bt.arreglada_",date.time, "_line906_CVestimation_ENDSCRIPT.RData", sep=""))
rm (date.time)


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

# pdf ("../../../GeoStatistics_toPublish/2_variogram_modelling/variog.plots.pdf", onefile=TRUE)
# par (family="f5", mfrow=c(3,4), oma=c(0.1,0.1,0.1,0.1), mar=c(2,3,2,2))
par (mfrow=c(3,4), oma=c(0.1,0.1,0.1,0.1), mar=c(2,3,2,2))
for (i in 1:nrow(milista)){
  tmp.variog <- get (paste(milista$variogs [i]))
  tmp.fit <- get(gsub("pred.", "fit.", milista$preds) [i])
  plot (gamma~dist, tmp.variog, xlab="", ylab="",col="black", pch=19, fg="grey", axes=FALSE, main=milista$preds[i])
  gamma.at <- seq(0, max(tmp.variog$gamma), 8000)
  gamma.labels <- gamma.at/1000
  axis(2,at=gamma.at, labels=gamma.labels, cex.axis=1.5, las=1)
  axis(1,at=seq(0,83340,18520), labels=seq(0,45,10), cex.axis=1.5, las=1)
  box(col="black")
  lines (variogramLine (tmp.fit, 220000, 200), col="grey")
  #a <- 12*c(1:10)
  #Para sacarlo en .pdf NO hace falta lo siguiente
  #if (i==a[1]|i==a[2]|i==a[3]|i==a[4]|i==a[5]|i==a[6]|i==a[7]|i==a[8]|i==a[9]|i==a[10]){
  #windows()
  #par (family="f5", mfrow=c(3,4), oma=c(0.1,0.1,0.1,0.1), mar=c(2,3,2,2))
  #}
}
# dev.off()

#------------------------------------------------------------------------------#
#---------------------------- End of script!!! --------------------------------#
#------------------------------------------------------------------------------#




