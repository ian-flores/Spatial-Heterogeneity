#sets the working directory
setwd("C:/Users/Ian/Dropbox/Spatial Heterogeneity Anoles/Data")

#Loads the data
sh <- read.csv("SH CSV.csv", header=T)

#shows the first 6 rows of data
head(sh)

#adds a column to use as index
sh["NUID"] <- seq(1,115, by=1)

#shows the first 6 rows of data
head(sh)

#creates a matrix using the x and y coordinates
coords <- cbind(sh$xcoord, sh$ycoord)

#loads the Spatial Tools library to calculate euclidian distance 
library("SpatialTools")

#calculates euclidian distance
d <- dist1(coords)

#creates a matrix with the index and infectious status
I <- cbind(sh$NUID, sh$Infeccion)

#Creates a 115x14 empty matrix
con.m = matrix(NA, nrow=115, ncol=13)

#Calculates the connectivity of the lizards 
for (j in 1:13) {
for (i in 1:115) {
con.m[i,j] = length(which(d[i,] < j))
}
}

#Prints the matrix of connectivity
con.m

#Creates a 115x13 empty matrix
coninf.m = matrix(NA, nrow=115, ncol=13)

#Calculates the connectivity of the lizards using only interactions
#with infectious status
for (j in 1:13) {
for (i in 1:115) {
coninf.m[i,j] = length(which(d[i,] < j & I[i,2]==1))
}
}

#Prints the matrix of connectivity with infections
coninf.m

### Logistic Regression ###

#Model for connectivity w SVL and Infections
ilr_1_svl = glm(I[,2]~coninf.m[,1]+sh$SVL)
ilr_1_svl
ilr_2_svl = glm(I[,2]~coninf.m[,2]+sh$SVL)
ilr_2_svl
ilr_3_svl = glm(I[,2]~coninf.m[,3]+sh$SVL)
ilr_3_svl
ilr_4_svl = glm(I[,2]~coninf.m[,4]+sh$SVL)
ilr_4_svl
ilr_5_svl = glm(I[,2]~coninf.m[,5]+sh$SVL)
ilr_5_svl
ilr_6_svl = glm(I[,2]~coninf.m[,6]+sh$SVL)
ilr_6_svl
ilr_7_svl = glm(I[,2]~coninf.m[,7]+sh$SVL)
ilr_7_svl
ilr_8_svl = glm(I[,2]~coninf.m[,8]+sh$SVL)
ilr_8_svl
ilr_9_svl = glm(I[,2]~coninf.m[,9]+sh$SVL)
ilr_9_svl
ilr_10_svl = glm(I[,2]~coninf.m[,10]+sh$SVL)
ilr_10_svl
ilr_11_svl = glm(I[,2]~coninf.m[,11]+sh$SVL)
ilr_11_svl
ilr_12_svl = glm(I[,2]~coninf.m[,12]+sh$SVL)
ilr_12_svl
ilr_13_svl = glm(I[,2]~coninf.m[,13]+sh$SVL)
ilr_13_svl

#Model for SVL only
lr_svl = glm(I[,2]~sh$SVL)
lr_svl

#Models for connectivity w Infections and without svl
ilr_1 = glm(I[,2]~coninf.m[,1])
ilr_1
ilr_2 = glm(I[,2]~coninf.m[,2])
ilr_2
ilr_3 = glm(I[,2]~coninf.m[,3])
ilr_3
ilr_4 = glm(I[,2]~coninf.m[,4])
ilr_4
ilr_5 = glm(I[,2]~coninf.m[,5])
ilr_5
ilr_6 = glm(I[,2]~coninf.m[,6])
ilr_6
ilr_7 = glm(I[,2]~coninf.m[,7])
ilr_7
ilr_8 = glm(I[,2]~coninf.m[,8])
ilr_8
ilr_9 = glm(I[,2]~coninf.m[,9])
ilr_9
ilr_10 = glm(I[,2]~coninf.m[,10])
ilr_10
ilr_11 = glm(I[,2]~coninf.m[,11])
ilr_11
ilr_12 = glm(I[,2]~coninf.m[,12])
ilr_12
ilr_13 = glm(I[,2]~coninf.m[,13])
ilr_13

#Models for connectivity without Infections and with SVL
lr_1_svl = glm(I[,2]~con.m[,1]+sh$SVL)
lr_1_svl
lr_2_svl = glm(I[,2]~con.m[,2]+sh$SVL)
lr_2_svl
lr_3_svl = glm(I[,2]~con.m[,3]+sh$SVL)
lr_3_svl
lr_4_svl = glm(I[,2]~con.m[,4]+sh$SVL)
lr_4_svl
lr_5_svl = glm(I[,2]~con.m[,5]+sh$SVL)
lr_5_svl
lr_6_svl = glm(I[,2]~con.m[,6]+sh$SVL)
lr_6_svl
lr_7_svl = glm(I[,2]~con.m[,7]+sh$SVL)
lr_7_svl
lr_8_svl = glm(I[,2]~con.m[,8]+sh$SVL)
lr_8_svl
lr_9_svl = glm(I[,2]~con.m[,9]+sh$SVL)
lr_9_svl
lr_10_svl = glm(I[,2]~con.m[,10]+sh$SVL)
lr_10_svl
lr_11_svl = glm(I[,2]~con.m[,11]+sh$SVL)
lr_11_svl
lr_12_svl = glm(I[,2]~con.m[,12]+sh$SVL)
lr_12_svl
lr_13_svl = glm(I[,2]~con.m[,13]+sh$SVL)
lr_13_svl

#Models for connectivity without Infections and without SVL
lr_1 = glm(I[,2]~con.m[,1])
lr_1
lr_2 = glm(I[,2]~con.m[,2])
lr_2
lr_3 = glm(I[,2]~con.m[,3])
lr_3
lr_4 = glm(I[,2]~con.m[,4])
lr_4
lr_5 = glm(I[,2]~con.m[,5])
lr_5
lr_6 = glm(I[,2]~con.m[,6])
lr_6
lr_7 = glm(I[,2]~con.m[,7])
lr_7
lr_8 = glm(I[,2]~con.m[,8])
lr_8
lr_9 = glm(I[,2]~con.m[,9])
lr_9
lr_10 = glm(I[,2]~con.m[,10])
lr_10
lr_11 = glm(I[,2]~con.m[,11])
lr_11
lr_12 = glm(I[,2]~con.m[,12])
lr_12
lr_13 = glm(I[,2]~con.m[,13])
lr_13

#loads library "AICcmodavg"
library("AICcmodavg")

#Calculates AIC for ALL the logistic regressions
aictab(list(ilr_1_svl,ilr_2_svl,ilr_3_svl,ilr_4_svl,ilr_5_svl,ilr_6_svl,
ilr_7_svl,ilr_8_svl,ilr_9_svl,ilr_10_svl,ilr_11_svl,ilr_12_svl,ilr_13_svl,
lr_svl,ilr_1,ilr_2,ilr_3,ilr_4,ilr_5,ilr_6,ilr_7,ilr_8,ilr_9,ilr_10,ilr_11,
ilr_12,ilr_13,lr_1_svl,lr_2_svl,lr_3_svl,lr_4_svl,lr_5_svl,lr_6_svl,lr_7_svl,
lr_8_svl,lr_9_svl,lr_10_svl,lr_11_svl,lr_12_svl,lr_13_svl,lr_1,lr_2,
lr_3,lr_4,lr_5,lr_6,lr_7,lr_8,lr_9,lr_10,lr_11,lr_12,lr_13),
modnames=c("Infected Radius 1 + SVL", "Infected Radius 2 + SVL",
"Infected Radius 3 + SVL","Infected Radius 4 + SVL",
"Infected Radius 5 + SVL","Infected Radius 6 + SVL",
"Infected Radius 7 + SVL","Infected Radius 8 + SVL",
"Infected Radius 9 + SVL","Infected Radius 10 + SVL",
"Infected Radius 11 + SVL","Infected Radius 12 + SVL",
"Infected Radius 13 + SVL","SVL", "Infected Radius 1", 
"Infected Radius 2","Infected Radius 3","Infected Radius 4",
"Infected Radius 5","Infected Radius 6","Infected Radius 7",
"Infected Radius 8","Infected Radius 9","Infected Radius 10",
"Infected Radius 11","Infected Radius 12","Infected Radius 13",
"Radius 1 + SVL","Radius 2 + SVL","Radius 3 + SVL","Radius 4 + SVL",
"Radius 5 + SVL","Radius 6 + SVL","Radius 7 + SVL","Radius 8 + SVL",
"Radius 9 + SVL","Radius 10 + SVL","Radius 11 + SVL","Radius 12 + SVL",
"Radius 13 + SVL", "Radius 1", "Radius 2","Radius 3","Radius 4",
"Radius 5","Radius 6","Radius 7","Radius 8","Radius 9","Radius 10",
"Radius 11","Radius 12","Radius 13"))
par(c(mfrow
logi.hist.plot(coninf.m[,13],I[,2])
logi.hist.plot(con.m[,13],I[,2])
logi.hist.plot(confinf.m[,13]+ sh$SVL,I[,2])
logi.hist.plot(con.m[,13]+sh$SVL,I[,2])
