# set working directory 
setwd("C:/Users/Ian/Dropbox/Spatial Heterogeneity Anoles/Data")
#load splancs
library(splancs)
#load data
infecciones <- read.csv("Spatial Heterogeneity CSV 11132015.csv", header=T)

#Positive 
positive <- subset(infecciones, Infeccion==1)

#transforms x
xtrans <- positive$xcoord- min(positive$xcoord)

#transforms the y
ytrans <- positive$ycoord - min(positive$ycoord)

#converts the coordinates into points 
positive_points <- as.points(xtrans, ytrans)

xmaxi <- max(xtrans)
xmini <- min(xtrans)
ymaxi <- max(ytrans)
ymini <- min(ytrans)

ibox_points <- as.points(c(xmini,xmini,xmaxi,xmaxi), 
                         c(ymini,ymaxi,ymini,ymaxi)) 
positive_polygon <- sbox(ibox_points, xfrac=0.025,yfrac=0.025)
plot(positive_polygon, asp=1, type='n')

image(kernel2d(positive_points, positive_polygon, h0=1, nx=100, ny=100), 
      add=TRUE, col=terrain.colors(15))
pointmap(positive_points, add=T)

null.csr<-csr(positive_polygon, 40) 
polymap(positive_polygon) 
pointmap(null.csr, pch=20, col="red") 
pointmap(positive_points,add=T)


lambda.m <-pdense(positive_points, positive_polygon) 

distances <- seq(from=0, to=22.5, by=0.25) 
iterations <- 1000 

positiveghat<-Ghat(positive_points, distances) 
plot(distances, positiveghat, type="line", xlab="distance", 
     ylab="Probability of Ghat") 

par(new = TRUE) 
null.csr.hat <-Ghat(null.csr, distances) 
plot(distances,null.csr.hat, type="line", col = "red", xlab=" ", ylab=" ") 
# este comando permitira poner en la misma "pagina" dos graficas 
# par(mfrow = c(1,2))

# K - Ripley's K (kHat1) for continuous forest 1 -c1 - second order 
#distribution 
kHatnull <-khat (null.csr, positive_polygon,distances) 
kHat1 <- khat(positive_points, positive_polygon, distances) 
# generate 95% CI for Ripley's K and modified Ripley's K 
kHat_CI <- Kenv.csr(length(positive_points[,1]), positive_polygon, 
                    iterations, distances) 
kHat_upper_CI <-kHat_CI$upper 
kHat_lower_CI <-kHat_CI$lower 
# plot Ripley's K observed - K expected 
plot(distances, kHat1-kHatnull, type="l", col="blue", xlab="Radius (m)", 
     ylab="K observed - K expected", main="Infection's K Observed - Expected", 
     ylim=c(-400,400)) 
abline(h=0, lty="dashed") 
lines(distances, kHat1-kHat_upper_CI, lty=2) 
lines(distances, kHat1-kHat_lower_CI, lty=2) 
# K - Modified Ripley's K (lhat1) for continuous forest 1 -c1 - second order distribution 
lHat1 <- sqrt(kHat1/pi) - distances
# generate 95% CI for Ripley's K and modified Ripley's K 
kHat_CI <- Kenv.csr(length(positive_points[,1]), positive_polygon, iterations, 
                    distances) 
lHat_upper_CI <- sqrt(kHat_CI$upper/pi) - distances
lHat_lower_CI <- sqrt(kHat_CI$lower/pi) - distances
##### plot the modified Ripley's K with its CI ####
plot(distances, lHat1, type="l", col="red", xlab="Radius (m)", 
     ylab="Modified Ripley's K", main="Spatial Pattern of Infections", 
     sub="LFDP", ylim=c(-8,8)) 
abline(h=0, lty="dashed")
lines(distances, lHat_upper_CI, lty=2) 
lines(distances, lHat_lower_CI, lty=2) 