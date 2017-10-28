
#Sets the working directory (Monophyllus)
setwd("C:/Users/Ian/Dropbox/Spatial Heterogeneity Anoles/Data")

#loads spatstat
library(spatstat)

#loads the data
anolis <- read.csv("SH CSV.csv", header=T)

#positive xcoords
positivexcoord <- anolis$xcoord[anolis$Infeccion==1]
positiveycoord <- anolis$ycoord[anolis$Infeccion==1]

positive <- cbind(positivexcoord, positiveycoord)

#creates a matrix with the coordinates
lagartijos <- cbind(anolis$xcoord, anolis$ycoord)

#creates the window in which the values are calculated
window_anolis <- owin(xrange=c(min(anolis$xcoord, na.rm=T), 
                               max(anolis$xcoord, na.rm=T)), 
                      yrange=c(min(anolis$ycoord, na.rm=T), 
                               max(anolis$ycoord, na.rm=T)))

# Transforms the lizard data to a ppp object
anolisppp <- as.ppp(lagartijos, W=window_anolis)

# Transforms the infected data to a ppp object
positiveppp <- as.ppp(positive, W=window_anolis)

# Generates the envelope for the general pop
k_env <- envelope(anolisppp, fun=Kest, nsim=1000)

# Calculates the K values for the infected lizards
K_inf <- Kest(positiveppp, correction="best")

#Plots the envelope
plot.fv(k_env, cbind(lo,hi)~r, ylim=c(0,1000), 
        main="Distribucion Espacial de las Infecciones", 
        legend=F, xlab="Distancia (m)")

# Plots the K values
points(K_inf$r,K_inf$iso, type="l", col="black")

# Plots the theorethical values
points(K_inf$r,K_inf$theo, type="l", col="red", lty=2)
legend(0,950, c("K Observado", "K Teorico", 
                "Intervalo de Confianza"), 
       fill=c("black", "red", "grey"), box.col="white")

dev.copy2pdf(file="Distribucion Espacial de las infecciones PDF")
