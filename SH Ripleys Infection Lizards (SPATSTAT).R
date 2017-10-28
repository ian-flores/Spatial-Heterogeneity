#Spatial Heterogeneity Ripley's K Infenction LIZARDS

#Sets the working directory (Monophyllus)
setwd("C:/Users/Ian/Dropbox/Spatial Heterogeneity Anoles/Data")

#Sets the working directory (Mac)
setwd("~/Dropbox/Spatial Heterogeneity Anoles/Data")

#loads spatstat
library(spatstat)

#loads the data
anolis <- read.csv("SH CSV.csv", header=T)

#plots the data
plot(anolis$xcoord, anolis$ycoord, 
     main="A. gundlachi spatial distribution",
     xlab="X", ylab="Y")

#positive xcoords
positivexcoord <- anolis$xcoord[anolis$Infeccion==1]
positiveycoord <- anolis$ycoord[anolis$Infeccion==1]
positivesvl <- anolis$SVL[anolis$Infeccion==1]

#creates a matrix with the coordinates
positive <- cbind(positivexcoord, positiveycoord)

#creates the window in which the values are calculated
window_positive <- owin(xrange=c(min(positivexcoord, na.rm=T), 
                               max(positivexcoord, na.rm=T)), 
                      yrange=c(min(positiveycoord, na.rm=T), 
                               max(positiveycoord, na.rm=T)))

# Transforms the data to a ppp object
positiveppp <- as.ppp(positive, W=window_anolis)

#Transforms the data to a ppp object, using SVL as a covariate

positive_svl <- as.ppp(positive, W=window_positive, marks=positivesvl)
# Estimates the K value for each value in the radius
K_inf <- Kest(positiveppp, correction="best")


# Plots the K values, without the covariate SVL
plot(envelope(positiveppp, fun=Kest, nsim=1000), 
     main="K Function for Infections", 
     xlab="Distance (m)", ylab="K(r)", sub="Without SVL")

# Plots the K values, with covariate SVL
plot(envelope(positive_svl, fun=Kest, nsim=1000), 
     main="K Function for Infections", 
     xlab="Distance (m)", ylab="K(r)", sub="Including SVL")

# Plots the L transformation, without the covariate SVL
plot(envelope(positiveppp, fun=Kest, nsim=1000, 
              transform = expression(sqrt(./pi))), 
     main="L Function for Infections", ylab="L(r)", sub="Without SVL")

# Plots the L transformation, with covariate SVL
plot(envelope(positive_svl, fun=Kest, nsim=1000, 
              transform = expression(sqrt(./pi))), 
     main="L Function for Infections", ylab="L(r)", sub="Including SVL")
