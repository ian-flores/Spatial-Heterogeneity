#Set the directory where the data is
setwd("C:/Users/Ian/Dropbox/Spatial Heterogeneity Anoles/Data")

# loads the data of the trees
trees <- read.csv("LFDP_Census5.csv", header=T)

# makes a subset of the trees 
arboles <- subset(trees, GX < 181.26 & GX > 129.88 & 
                    GY < 137.99 & GY > 74.73)

#Loads the Spatial data
SH <- read.csv("SH CSV.csv", header=T)

#Infected lizards
positive <- subset(SH, Infeccion==1)

#Exports the x coordinates in SH
xcoord <-SH$xcoord

#Exports the y coordinates in SH
ycoord <-SH$ycoord

#Plots all the trees in the LFDP between x=(130,180) and y=(70,120)
plot(arboles$GX, arboles$GY, xlab="",
     ylab="", main="", sub="", pch=21, col="grey", bg="grey", axes=F)
box()

#Adds the trees where the anoles where sampled
#PCH 21 and COL and BG
points(xcoord, ycoord, pch=21, bg="white", col="black")

#Adds the points of the infected lizards
points(positive$xcoord, positive$ycoord, pch=21, bg="black", col="black")

#Adds a legend
legend("topright", c("Arboles", "Lagartijos", "Infecciones"), 
       fill=c("grey", "white", "black"))