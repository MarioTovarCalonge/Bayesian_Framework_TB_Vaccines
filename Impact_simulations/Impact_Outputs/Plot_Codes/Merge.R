#!/usr/bin/Rscript

combimpacts <- read.table("Model_Combined_IRR.txt", header=TRUE)

#India
dIndia <- read.table("../Impacts_models_per_country/impact_India.txt", header=FALSE)

aux <- subset(combimpacts, country=="India")
dIndia[nrow(dIndia) + 1,] = c("All", 0, 0, 0, aux$IRR, aux$IRR_low, aux$IRR_hi)

colnames(dIndia) <- c("Model", "Imp", "Imp_low", "Imp_hi", "IRR", "IRR_low", "IRR_hi")
dIndia$Country = rep("India", 8)
print(dIndia)




#Indonesia
dIndonesia <- read.table("../Impacts_models_per_country/impact_Indonesia.txt", header=FALSE)
aux <- subset(combimpacts, country=="Indonesia")
dIndonesia[nrow(dIndonesia) + 1,] = c("All", 0, 0, 0, aux$IRR, aux$IRR_low, aux$IRR_hi)
colnames(dIndonesia) <- c("Model", "Imp", "Imp_low", "Imp_hi", "IRR", "IRR_low", "IRR_hi")
dIndonesia$Country = rep("Indonesia", 8)
print(dIndonesia)




#Ethiopia
dEthiopia <- read.table("../Impacts_models_per_country/impact_Ethiopia.txt", header=FALSE)
aux <- subset(combimpacts, country=="Ethiopia")
dEthiopia[nrow(dEthiopia) + 1,] = c("All", 0, 0, 0,aux$IRR, aux$IRR_low, aux$IRR_hi)
colnames(dEthiopia) <- c("Model", "Imp", "Imp_low", "Imp_hi", "IRR", "IRR_low", "IRR_hi")
dEthiopia$Country = rep("Ethiopia", 8)
print(dEthiopia)


#Ahora ya tenemos la información necesaria para cada plot.
library(dplyr)
df_final <- bind_rows(dIndia, dIndonesia, dEthiopia)
write.table(df_final, "Impact_all_countries.txt", quote=FALSE, row.names=FALSE)
