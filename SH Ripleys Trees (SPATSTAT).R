#Spatial Heterogeneity Ripley's K Trees

# set working directory (Monophyllus)
setwd("C:/Users/Ian/Dropbox/Spatial Heterogeneity Anoles/Data")

#Sets the working directory (Mac)
setwd("~/Dropbox/Spatial Heterogeneity Anoles/Data")

#loads spatstat
library(spatstat)

# loads the data of the trees
trees <- read.csv("LFDP_Census5.csv", header=T)

# makes a subset of the trees sampled
arboles <- subset(trees, GX < 181.26 & GX > 129.88 & GY < 137.99 
                  & GY > 74.73)

# Plots the pattern
plot(arboles$GX, arboles$GY)

#Creates a matrix with the coordinates
arbolito <- cbind(arboles$GX, arboles$GY)

#creates the window in which the values are calculated
window_trees <- owin(xrange=c(min(arboles$GX, na.rm=T), max(arboles$GX, na.rm=T)), 
                      yrange=c(min(arboles$GY, na.rm=T), max(arboles$GY, na.rm=T)))

# Transforms the tree data to a ppp object
arbolesppp <- as.ppp(arbolito, W=window_trees)

# Transforms the tree data to a ppp object, adding DBH as a co-variate
arbolesgrandes <- as.ppp(arbolito, W=window_trees, marks=arboles$DBH)

# Estimates the K value for each value in the radius
K_trees <- Kest(arbolesppp, correction="best")

#Creates a window to show both plots in a single image
par(mfrow=c(1,2))

#Plots the K values, without the covariate
plot(envelope(arbolesppp, fun=Kest, nsim=100), 
     main="Distribucion Espacial de los Arboles", 
              xlab="Distancia (m)", ylab="K(r)", legend=F)

legend(0,525, c("K Observado", "K Teorico", 
                "Intervalo de Confianza"), 
       fill=c("black", "red", "grey"), box.col="white")

dev.copy2pdf(file="Distribucion Espacial de los arboles PDF")


#Plots the K values, with DBH
plot(envelope(arbolesgrandes, fun=Kest, nsim=100), main="K Function for Trees", 
     xlab="Distance (m)", ylab="K(r)", sub="Including DBH")

# Plots the L transformation, without DBH 
plot(envelope(arbolesppp, fun=Kest, nsim=100, transform = expression(sqrt(./pi))), 
     main="L Function for Trees", ylab="L(r)", sub="Without DBH")

# Plots the L transformation, with DBH 
plot(envelope(arbolesgrandes, fun=Kest, nsim=100, transform = expression(sqrt(./pi))), 
     main="L Function for Trees", ylab="L(r)", sub="Including DBH")
