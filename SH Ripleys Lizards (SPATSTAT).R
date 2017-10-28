#Spatial Heterogeneity Ripley's K LIZARDS

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

#creates a matrix with the coordinates
lagartijos <- cbind(anolis$xcoord, anolis$ycoord)

#creates the window in which the values are calculated
window_anolis <- owin(xrange=c(min(anolis$xcoord, na.rm=T), 
                               max(anolis$xcoord, na.rm=T)), 
             yrange=c(min(anolis$ycoord, na.rm=T), 
                      max(anolis$ycoord, na.rm=T)))

# Transforms the data to a ppp object
anolisppp <- as.ppp(lagartijos, W=window_anolis)

#Transforms the data to a ppp object, using SVL as a covariate

anolis_grandes <- as.ppp(lagartijos, W=window_anolis, marks=anolis$SVL)
# Estimates the K value for each value in the radius
# No se puede expandir el r??
K <- Kest(anolisppp, correction="best")

#Creates a window to show both plots in a single image
par(mfrow=c(1,2))

z <- envelope(anolisppp, fun=Kest, nsim=1000)

# Plots the K values, without the covariate SVL
plot(z, main="Distribucion Espacial de los Lagartijos" , 
     xlab="Distancia (m)", ylab="K(r)", legend=F)

legend(0,600, c("K Observado", "K Teorico", 
                "Intervalo de Confianza"), 
       fill=c("black", "red", "grey"), box.col="white")

dev.copy2pdf(file="Distribucion Espacial de los Lagartijos PDF")
# Plots the K values, with covariate SVL
plot(envelope(anolis_grandes, fun=Kest, nsim=1000), 
     main="K Function for A. gundlachi", 
     xlab="Distance (m)", ylab="K(r)", sub="Including SVL")

# Plots the L transformation, without the covariate SVL
plot(envelope(anolisppp, fun=Kest, nsim=1000, 
              transform = expression(sqrt(./pi))), 
     main="L Function for A. gundlachi", ylab="L(r)", sub="Without SVL")

# Plots the L transformation, with covariate SVL
plot(envelope(anolis_grandes, fun=Kest, nsim=1000, 
              transform = expression(sqrt(./pi))), 
     main="L Function for A. gundlachi", ylab="L(r)", sub="Including SVL")
