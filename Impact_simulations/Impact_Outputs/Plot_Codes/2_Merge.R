#!/usr/bin/Rscript

#Neonatos
target <- 0

combimpacts <- read.table(paste0("Model_Combined_IRR_T_", target, ".txt"), header=TRUE)

#India
dIndia <- read.table(paste0("../Impacts_models_per_country/impact_India_T_", target, ".txt"), header=FALSE)

aux <- subset(combimpacts, country=="India")
dIndia[nrow(dIndia) + 1,] = c("All", 0, 0, 0, aux$IRRdist, aux$IRRdist_low, aux$IRRdist_hi)

colnames(dIndia) <- c("Model", "Imp", "Imp_low", "Imp_hi", "IRR", "IRR_low", "IRR_hi")
dIndia$Country = rep("India", 8)
print(dIndia)




#Indonesia
dIndonesia <- read.table(paste0("../Impacts_models_per_country/impact_Indonesia_T_", target, ".txt"), header=FALSE)
aux <- subset(combimpacts, country=="Indonesia")
dIndonesia[nrow(dIndonesia) + 1,] = c("All", 0, 0, 0, aux$IRRdist, aux$IRRdist_low, aux$IRRdist_hi)
colnames(dIndonesia) <- c("Model", "Imp", "Imp_low", "Imp_hi", "IRR", "IRR_low", "IRR_hi")
dIndonesia$Country = rep("Indonesia", 8)
print(dIndonesia)




#Ethiopia
dEthiopia <- read.table(paste0("../Impacts_models_per_country/impact_Ethiopia_T_", target, ".txt"), header=FALSE)
aux <- subset(combimpacts, country=="Ethiopia")
dEthiopia[nrow(dEthiopia) + 1,] = c("All", 0, 0, 0,aux$IRRdist, aux$IRRdist_low, aux$IRRdist_hi)
colnames(dEthiopia) <- c("Model", "Imp", "Imp_low", "Imp_hi", "IRR", "IRR_low", "IRR_hi")
dEthiopia$Country = rep("Ethiopia", 8)
print(dEthiopia)


#Ahora ya tenemos la información necesaria para cada plot.
library(dplyr)
df_final <- bind_rows(dIndia, dIndonesia, dEthiopia)
write.table(df_final, paste0("Impact_all_countries_T_", target, ".txt"), quote=FALSE, row.names=FALSE)

#--------------------------------------------------------------#

target <- 3

combimpacts <- read.table(paste0("Model_Combined_IRR_T_", target, ".txt"), header=TRUE)

#India
dIndia <- read.table(paste0("../Impacts_models_per_country/impact_India_T_", target, ".txt"), header=FALSE)

aux <- subset(combimpacts, country=="India")
dIndia[nrow(dIndia) + 1,] = c("All", 0, 0, 0, aux$IRRdist, aux$IRRdist_low, aux$IRRdist_hi)

colnames(dIndia) <- c("Model", "Imp", "Imp_low", "Imp_hi", "IRR", "IRR_low", "IRR_hi")
dIndia$Country = rep("India", 8)
print(dIndia)




#Indonesia
dIndonesia <- read.table(paste0("../Impacts_models_per_country/impact_Indonesia_T_", target, ".txt"), header=FALSE)
aux <- subset(combimpacts, country=="Indonesia")
dIndonesia[nrow(dIndonesia) + 1,] = c("All", 0, 0, 0, aux$IRRdist, aux$IRRdist_low, aux$IRRdist_hi)
colnames(dIndonesia) <- c("Model", "Imp", "Imp_low", "Imp_hi", "IRR", "IRR_low", "IRR_hi")
dIndonesia$Country = rep("Indonesia", 8)
print(dIndonesia)




#Ethiopia
dEthiopia <- read.table(paste0("../Impacts_models_per_country/impact_Ethiopia_T_", target, ".txt"), header=FALSE)
aux <- subset(combimpacts, country=="Ethiopia")
dEthiopia[nrow(dEthiopia) + 1,] = c("All", 0, 0, 0,aux$IRRdist, aux$IRRdist_low, aux$IRRdist_hi)
colnames(dEthiopia) <- c("Model", "Imp", "Imp_low", "Imp_hi", "IRR", "IRR_low", "IRR_hi")
dEthiopia$Country = rep("Ethiopia", 8)
print(dEthiopia)


#Ahora ya tenemos la información necesaria para cada plot.
library(dplyr)
df_final <- bind_rows(dIndia, dIndonesia, dEthiopia)
write.table(df_final, paste0("Impact_all_countries_T_", target, ".txt"), quote=FALSE, row.names=FALSE)
