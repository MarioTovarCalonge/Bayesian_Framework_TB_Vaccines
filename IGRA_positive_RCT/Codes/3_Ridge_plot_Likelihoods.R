#!/usr/bin/Rscript


library(ggplot2)
library(ggridges)
library(viridis)

library(dplyr)
bootmax <- 5000
aux <- data.frame(ID=1:bootmax)

vec_fact <- numeric()
vec_data <- numeric()
df_aux <- data.frame(it=1:bootmax, pm1=rep(0, bootmax), pm2=rep(0, bootmax), pm3=rep(0, bootmax), pm4=rep(0, bootmax), pm5=rep(0, bootmax), pm6=rep(0, bootmax), pm7=rep(0, bootmax))

#Load data
for(i in 1:7){
    dat <- read.table(paste0("Outputs/Model_", i, "/like_CI_", i, "_A_Mixed_pop.txt"), header=T)
    df_aux[,i+1] <- dat$Data
}
df_aux$norm=rowSums(df_aux[,2:8])
df_aux[, 2:8] <- df_aux[,2:8]/df_aux$norm

#mix data
for(i in 1:7){
    aux[ , ncol(aux) + 1] <- df_aux[, i+1] #Append new column
    colnames(aux)[ncol(aux)] <- paste0("Model_", i, "_data") #Rename column name
    vec_fact <- c(vec_fact, rep(paste0("M",i), bootmax))
    vec_data <- c(vec_data, df_aux[, i+1])
}

#Ahora aux contiene por columnas los resultados de cada modelo.
#Lo transformo en otro data frame
df_plot <- data.frame(Model=as.factor(vec_fact), Data=vec_data)
head(df_plot)

p <- ggplot(df_plot, aes(x = Data, y = Model, fill = Model)) +
    geom_density_ridges() +
    theme_ridges() +
    labs(x=expression(paste("L(m|V", E["dis"], ")")), y="") +
    scale_fill_viridis(option='cividis', discrete=TRUE) +
    theme(legend.position="none")
    
pdf("model_posterior.pdf", height=4, width=6)
print(p)
dev.off()



    
    

