#Spatial Heterogeneity for Anoles (Spatial EDA)

# set working directory 
setwd("C:/Users/Ian/Dropbox/Spatial Heterogeneity Anoles/Data")
#load splancs
library(splancs)
#load data
malaria <- read.csv("Spatial Heterogeneity CSV 11132015.csv", header=T)

#transforms the x
xtrans <- malaria$xcoord - min(malaria$xcoord)

#transforms the y 
ytrans <- malaria$ycoord - min(malaria$ycoord)

#converts the coordinates into points 
malaria_points <- as.points(xtrans, ytrans)

xmax <- max(xtrans)
xmin <- min(xtrans)
ymax <- max(ytrans)
ymin <- min(ytrans)

mbox_points <- as.points(c(xmin,xmin,xmax,xmax), c(ymin,ymax,ymin,ymax)) 
malaria_polygon <- sbox(mbox_points, xfrac=0.025,yfrac=0.025)
plot(malaria_polygon, asp=1, type='n')

image(kernel2d(malaria_points, malaria_polygon, h0=1, nx=100, ny=100), 
      add=TRUE, col=terrain.colors(15))
pointmap(malaria_points, add=T)

null.csr<-csr(malaria_polygon, 116) 
polymap(malaria_polygon) 
pointmap(null.csr, pch=20, col="red") 
pointmap(malaria_points,add=T)


lambda.m <-pdense(malaria_points, malaria_polygon) 

distances <- seq(from=0, to=22.5, by=0.5) 
iterations <- 1000 

malariaghat<-Ghat(malaria_points, distances) 
plot(distances, malariaghat, type="line", xlab="distance", ylab="Probability of Ghat") 

par(new = TRUE) 
null.csr.hat <-Ghat(null.csr, distances) 
plot(distances,null.csr.hat, type="line", col = "red", xlab=" ", ylab=" ") 
# este comando permitira poner en la misma "pagina" dos graficas 
# par(mfrow = c(1,2))

# K - Ripley's K (kHat1) for continuous forest 1 -c1 - second order distribution 
kHatnull <-khat (null.csr, malaria_polygon,distances) 
kHat1 <- khat(malaria_points, malaria_polygon, distances) 
# generate 95% CI for Ripley's K and modified Ripley's K 
kHat_CI <- Kenv.csr(length(malaria_points[,1]), malaria_polygon, iterations, distances) 
kHat_upper_CI <-kHat_CI$upper 
kHat_lower_CI <-kHat_CI$lower 
# plot Ripley's K observed - K expected 
plot(distances, kHat1-kHatnull, type="l", col="blue", xlab="Radius (m)", 
     ylab="K observed - K expected", main="Lizard's K Observed - Expected") 
abline(h=0, lty="dashed") 
lines(distances, kHat1-kHat_upper_CI, lty=2) 
lines(distances, kHat1-kHat_lower_CI, lty=2) 
# K - Modified Ripley's K (lhat1) for continuous forest 1 -c1 - second order distribution 
lHat1 <- sqrt(kHat1/pi) - distances
# generate 95% CI for Ripley's K and modified Ripley's K 
kHat_CI <- Kenv.csr(length(malaria_points[,1]), malaria_polygon, iterations, distances) 
lHat_upper_CI <- sqrt(kHat_CI$upper/pi) - distances
lHat_lower_CI <- sqrt(kHat_CI$lower/pi) - distances
##### plot the modified Ripley's K with its CI ####
plot(distances, lHat1, type="l", col="red", xlab="Radius (m)", 
     ylab="Modified Ripley's K", main="Spatial Pattern of Anole Lizards", 
     sub="LFDP", ylim=c(-2, 8)) 
abline(h=0, lty="dashed")
lines(distances, lHat_upper_CI, lty=2) 
lines(distances, lHat_lower_CI, lty=2) 

