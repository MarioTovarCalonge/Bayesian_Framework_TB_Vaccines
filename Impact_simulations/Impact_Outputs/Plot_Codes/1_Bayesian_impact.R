#!/usr/bin/Rscript

countries <- c('India', 'Indonesia', 'Ethiopia')
Ls <- read.table("Data/Likelihood_values.txt", header=T)
LM <- Ls$LM
norm <- sum(LM)
enditer <- c(500, 500, 500)
names(enditer) <- countries

#Neonatos
target <- 0

df_out <- data.frame(country=numeric(), Model=numeric(), IRR=numeric(), IRR_low=numeric(), IRR_hi=numeric())
print(df_out)

for (country in countries){
    IRR_M <- list()
    for(model in 1:7){
        df_aux <- read.table(paste0('../iteration_impacts/', country, '_M', model, '_T_', target, '_dataimp.txt'), header=T)
        IRR <- (df_aux$i_con - df_aux$i_vac)/df_aux$i_con * 100
        IRR_M[[model]] <- IRR
    }
    
    combIRR <- numeric(enditer[country])
    for (i in 1:enditer[country]){
        wsum <- IRR_M[[1]][i]*LM[1] + IRR_M[[2]][i]*LM[2] + IRR_M[[3]][i]*LM[3] + IRR_M[[4]][i]*LM[4] + IRR_M[[5]][i]*LM[5] + IRR_M[[6]][i]*LM[6] + IRR_M[[7]][i]*LM[7]
        combIRR[i] <- wsum/norm
    }
    
    combIRR_low <- quantile(combIRR, .025)    
    combIRR_up <- quantile(combIRR, .975)
    combIRR_median <- quantile(combIRR, .5)
    print(c(country, "All", combIRR_median, combIRR_low, combIRR_up))
    
    df_out[nrow(df_out) + 1,] = c(country, "All", combIRR_median, combIRR_low, combIRR_up)
    print(df_out)
}
    
write.table(df_out, paste0("Model_Combined_IRR_T_", target, ".txt"), quote=FALSE, row.names=FALSE)


#Adol
target <- 3

df_out <- data.frame(country=numeric(), Model=numeric(), IRR=numeric(), IRR_low=numeric(), IRR_hi=numeric())
print(df_out)

for (country in countries){
    IRR_M <- list()
    for(model in 1:7){
        df_aux <- read.table(paste0('../iteration_impacts/', country, '_M', model, '_T_', target, '_dataimp.txt'), header=T)
        IRR <- (df_aux$i_con - df_aux$i_vac)/df_aux$i_con * 100
        IRR_M[[model]] <- IRR
    }
    
    combIRR <- numeric(enditer[country])
    for (i in 1:enditer[country]){
        wsum <- IRR_M[[1]][i]*LM[1] + IRR_M[[2]][i]*LM[2] + IRR_M[[3]][i]*LM[3] + IRR_M[[4]][i]*LM[4] + IRR_M[[5]][i]*LM[5] + IRR_M[[6]][i]*LM[6] + IRR_M[[7]][i]*LM[7]
        combIRR[i] <- wsum/norm
    }
    
    combIRR_low <- quantile(combIRR, .025)
    combIRR_up <- quantile(combIRR, .975)
    combIRR_median <- quantile(combIRR, .5)
    print(c(country, "All", combIRR_median, combIRR_low, combIRR_up))
    
    df_out[nrow(df_out) + 1,] = c(country, "All", combIRR_median, combIRR_low, combIRR_up)
    print(df_out)
}
    
write.table(df_out, paste0("Model_Combined_IRR_T_", target, ".txt"), quote=FALSE, row.names=FALSE)
