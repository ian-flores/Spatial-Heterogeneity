#Spatial Heterogeneity Ripley's K Female LIZARDS

#Sets the working directory (Monophyllus)
setwd("C:/Users/Ian/Dropbox/Spatial Heterogeneity Anoles/Data")

#Sets the working directory (Mac)
setwd("~/Dropbox/Spatial Heterogeneity Anoles/Data")

#loads spatstat
library(spatstat)

#loads the data
anolis <- read.csv("SH CSV.csv", header=T)

#plots the general data
plot(anolis$xcoord, anolis$ycoord, 
     main="A. gundlachi spatial distribution",
     xlab="X", ylab="Y")

#subsets the data for females 
hembrasxcoord <- anolis$xcoord[anolis$Sex=="Female"]
hembrasycoord <- anolis$ycoord[anolis$Sex=="Female"]
hembras_svl <- anolis$SVL[anolis$Sex=="Female"]

#plots the data for females
plot(hembrasxcoord, hembrasycoord, 
     main="Female A. gundlachi spatial distribution", 
     xlab="X", ylab="Y")

#creates a matrix with the coordinates
lagartijos_hembras <- cbind(hembrasxcoord, hembrasycoord)

#creates the window in which the values are calculated
window_anolis_hembras <- owin(xrange=c(min(hembrasxcoord, na.rm=T), 
                                       max(hembrasxcoord, na.rm=T)), 
                             yrange=c(min(hembrasycoord, na.rm=T), 
                                      max(hembrasycoord, na.rm=T)))

# Transforms the data to a ppp object
anolis_hembrasppp <- as.ppp(lagartijos_hembras, W=window_anolis_hembras)
anolis_hembras_svlppp <- as.ppp(lagartijos_hembras, 
                                W=window_anolis_hembras, 
                                marks=hembras_svl)

# Estimates the K value for each value in the radius
K <- Kest(anolis_hembrasppp, correction="best")

#Creates a window to show both plots in a single image
par(mfrow=c(2,2))

# Plots the K values, without the covariate SVL
plot(envelope(anolis_hembrasppp, fun=Kest, nsim=1000), 
     main="K Function for female A. gundlachi", 
     xlab="Distance (m)", ylab="K(r)", sub="Without SVL")

# Plots the K values, with covariate SVL
plot(envelope(anolis_hembras_svlppp, fun=Kest, nsim=1000), 
     main="K Function for female A. gundlachi", 
     xlab="Distance (m)", ylab="K(r)", sub="Including SVL")

# Plots the L transformation, without the covariate SVL
plot(envelope(anolis_hembrasppp, fun=Kest, nsim=1000, 
              transform = expression(sqrt(./pi))), 
     main="L Function for female A. gundlachi", ylab="L(r)", 
     xlab="Distance (m)",sub="Without SVL")

# Plots the L transformation, with covariate SVL
plot(envelope(anolis_hembras_svlppp, fun=Kest, nsim=1000, 
              transform = expression(sqrt(./pi))), 
     main="L Function for female A. gundlachi", ylab="L(r)", 
     xlab="Distance(m)", sub="Including SVL")
