#-------------------------------------------------------------------------------
# Set coords & project
#-------------------------------------------------------------------------------
setproj <- function (mydata, utm) {
  library (rgdal)
  if (utm=="31N"){
    coordinates(mydata) <- ~x31+y31
    proj4string(mydata) <- CRS("+proj=utm +zone=31 +ellps=WGS84")
  }
  if (utm=="30N"){
    coordinates(mydata) <- ~x30+y30
    proj4string(mydata) <- CRS("+proj=utm +zone=30 +ellps=WGS84")
  }
  mydata
}
#-------------------------------------------------------------------------------