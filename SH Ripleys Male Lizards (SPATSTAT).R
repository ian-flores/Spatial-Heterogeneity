#Spatial Heterogeneity Ripley's K Male LIZARDS

#Sets the working directory (Monophyllus)
setwd("C:/Users/Ian/Dropbox/Spatial Heterogeneity Anoles/Data")

#Sets the working directory (Mac)
setwd("~/Dropbox/Spatial Heterogeneity Anoles/Data")

#loads spatstat
library(spatstat)

#loads the data
anolis <- read.csv("SH CSV.csv", header=T)

#plots the general data
plot(anolis$xcoord, anolis$ycoord, main="A. gundlachi spatial distribution",
     xlab="X", ylab="Y")

#subsets the data for males 
machosxcoord <- anolis$xcoord[anolis$Sex=="Male"]
machosycoord <- anolis$ycoord[anolis$Sex=="Male"]
machos_svl <- anolis$SVL[anolis$Sex=="Male"]
#plots the data for males
plot(machosxcoord, machosycoord, main="Male A. gundlachi spatial distribution", 
     xlab="X", ylab="Y")

#creates a matrix with the coordinates
lagartijos_machos <- cbind(machosxcoord, machosycoord)

#creates the window in which the values are calculated
window_anolis_machos <- owin(xrange=c(min(machosxcoord, na.rm=T), max(machosxcoord, na.rm=T)), 
                      yrange=c(min(machosycoord, na.rm=T), max(machosycoord, na.rm=T)))

# Transforms the data to a ppp object
anolis_machosppp <- as.ppp(lagartijos_machos, W=window_anolis_machos)
anolis_machos_svlppp <- as.ppp(lagartijos_machos, W=window_anolis_machos, marks=machos_svl)

# Estimates the K value for each value in the radius
K <- Kest(anolis_machosppp, correction="best")

#Creates a window to show both plots in a single image
par(mfrow=c(2,2))

# Plots the K values, without the covariate SVL
plot(envelope(anolis_machosppp, fun=Kest, nsim=1000), 
     main="K Function for male A. gundlachi", 
     xlab="Distance (m)", ylab="K(r)", sub="Without SVL")

# Plots the K values, with covariate SVL
plot(envelope(anolis_machos_svlppp, fun=Kest, nsim=1000), 
     main="K Function for male A. gundlachi", 
     xlab="Distance (m)", ylab="K(r)", sub="Including SVL")

# Plots the L transformation, without the covariate SVL
plot(envelope(anolis_machosppp, fun=Kest, nsim=1000, 
              transform = expression(sqrt(./pi))), 
     main="L Function for male A. gundlachi", ylab="L(r)", 
     xlab="Distance (m)",sub="Without SVL")

# Plots the L transformation, with covariate SVL
plot(envelope(anolis_machos_svlppp, fun=Kest, nsim=1000, 
              transform = expression(sqrt(./pi))), 
     main="L Function for male A. gundlachi", ylab="L(r)", 
     xlab="Distance(m)", sub="Including SVL")