#Spatial Heterogeneity for Trees

# set working directory 
setwd("C:/Users/Ian/Dropbox/Spatial Heterogeneity Anoles/Data")

# load splancs
library(splancs)

# loads the data of the trees
trees <- read.csv("LFDP_Census5.csv", header=T)

# makes a subset of the trees 
arboles <- subset(trees, GX < 204.9 & GX > 129.88 & GY < 137.99 
                  & GY > 74.73)

#transforms the x
xtrans <- arboles$GX - min(arboles$GX)

#transforms the y 
ytrans <- arboles$GY - min(arboles$GY)

#converts the coordinates into points 
tree_points <- as.points(xtrans, ytrans)

xmax <- max(xtrans)
xmin <- min(xtrans)
ymax <- max(ytrans)
ymin <- min(ytrans)

# Generates a bounding polygon 
arbolesbox_points <- as.points(c(xmin,xmin,xmax,xmax), 
                               c(ymin,ymax,ymin,ymax)) 

arboles_polygon <- sbox(arbolesbox_points, xfrac=0.045,
                        yfrac=0.045)
plot(arboles_polygon, asp=1, type='n')

image(kernel2d(tree_points, arboles_polygon, h0=1, nx=100, ny=100), 
      add=TRUE, col=terrain.colors(15))
pointmap(tree_points, add=T)

null.csr<-csr(arboles_polygon, 1294) 
polymap(arboles_polygon) 
pointmap(null.csr, pch=20, col="red") 
pointmap(tree_points,add=T)


lambda.c <- pdense(tree_points, arboles_polygon) 

distances <- seq(from=0, to=22.5, by=0.5) 
iterations <- 1000 

arbolesghat<-Ghat(tree_points, distances) 
plot(distances, arbolesghat, type="line", xlab="distance", 
     ylab="Probability of Ghat") 

par(new = TRUE) 
null.csr.hat <-Ghat(null.csr, distances) 
points(distances,null.csr.hat, type="line", col = "red", xlab=" ", ylab=" ") 
# este comando permitira poner en la misma "pagina" dos graficas 
# par(mfrow = c(1,2))

# K - Ripley's K (kHat1) for continuous forest 1 -c1 - second order distribution 
kHatnull <-khat (null.csr, arboles_polygon,distances) 
kHat1 <- khat(tree_points, arboles_polygon, distances) 
# generate 95% CI for Ripley's K and modified Ripley's K 
kHat_CI <- Kenv.csr(length(arboles$GX), 
                    arboles_polygon, iterations, distances) 
kHat_upper_CI <-kHat_CI$upper 
kHat_lower_CI <-kHat_CI$lower 
# plot Ripley's K observed - K expected 
plot(distances, kHat1-kHatnull, type="l", col="blue", 
     xlab="Radius (m)", 
     ylab="K observed - K expected", 
     main="Trees' K Observed - Expected") 
abline(h=0, lty="dashed") 
lines(distances, kHat1-kHat_upper_CI, lty=2) 
lines(distances, kHat1-kHat_lower_CI, lty=2) 
# K - Modified Ripley's K (lhat1) for continuous forest 1 -c1 - second order distribution 
lHat1 <- sqrt(kHat1/pi)-distances 
# generate 95% CI for Ripley's K and modified Ripley's K 
kHat_CI <- Kenv.csr(length(arboles$GX), arboles_polygon, iterations, 
                    distances) 
lHat_upper_CI <- sqrt(kHat_CI$upper/pi) - distances
lHat_lower_CI <- sqrt(kHat_CI$lower/pi) - distances
##### plot the modified Ripley's K with its CI ####
plot(distances, lHat1, type="l", col="red", xlab="Radius (m)", 
     ylab="Modified Ripley's K", main="Spatial Pattern of Trees", 
     sub="LFDP", ylim=c(-2,2)) 
abline(h=0, lty="dashed")
lines(distances, lHat_upper_CI, lty=2) 
lines(distances, lHat_lower_CI, lty=2) 

