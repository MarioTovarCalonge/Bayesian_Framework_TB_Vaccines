#!/usr/bin/Rscript

args<-commandArgs(TRUE)
country <- args[1]
mode <- args[2]
target <- as.integer(args[3])


M1 = F; M2 = F; M3 = F; M4 = F; M5 = F; M6 = F; M7 = F
mode_id = 0
if(mode=='M1'){
    M1 <- T
    mode_id = 1
} else if(mode=='M2') {
    M2 <- T
    mode_id = 2
} else if(mode=='M3'){
    M3 <- T
    mode_id = 3
} else if(mode=='M4'){
    M4 <- T
    mode_id = 4
} else if(mode=='M5'){
    M5 <- T
    mode_id = 5
} else if(mode=='M6'){
    M6 <- T
    mode_id = 6
} else if(mode=='M7'){
    M7 <- T
    mode_id = 7
}

#target <- 0
target_group <- 'LR'
waning <- 0.00

aux <- read.table("Eps_estim_A_Mixed_pop.txt", header=FALSE, skip=1)
meds <- as.vector(unlist(aux$V5))
lows <- as.vector(unlist(aux$V6))
his <- as.vector(unlist(aux$V7))

med <- meds[mode_id]
low <- lows[mode_id]
hi <-  his[mode_id]

print(med)
print(low)
print(hi)

sd = (hi - low)/3.92

cov_save = numeric(500)

library(truncnorm)

for(iter in 0:500){
    cov = rtruncnorm(1, a=0, b=1, mean=med, sd = sd)
    cov_save[iter] = cov
    
    if(iter==0){
        output_path <- "../Inputs/General_Inputs/perfil.txt"
        #output_path <- paste0("Profiles/M", mode_id, "profile.txt")
        cov <- med
        if(med>1){
            cov=1
        }
    } else{
        output_path <- paste0("../Inputs_Iter/", country, "/Perfil/perfil_", sprintf("%03d", iter), ".txt")
        output_dir <- paste0("../Inputs_Iter/", country, "/Perfil/")
        if (!dir.exists(output_dir)){
            dir.create(output_dir, recursive=TRUE)
        }
    }
    
    line <- sprintf("target = %d", target)
    write(line, output_path)
    line <- sprintf("waning_h = %f", 10)
    write(line, output_path, append=TRUE)
    line <- "transition age S_L_R efficacy e_low e_hi coverage c_low c_hi waning w_low w_hi blocking b_low b_hi"
    write(line, output_path, append=TRUE)

    if(M1){
        line <- paste0("E_p_AON ", sprintf("%d ", target), target_group, " 1.0 1.0 1.0 ", sprintf("%f %f %f", cov, cov, cov), sprintf(" %.2f %.2f %.2f", waning, waning, waning), " 0.0 0.0 0.0")
        write(line, output_path, append=TRUE)
        
    } 

    if(M2){
        line <- paste0("E_q ", sprintf("%d ", target), target_group, " 1.0 1.0 1.0 ", sprintf("%f %f %f", cov, cov, cov), sprintf(" %.2f %.2f %.2f", waning, waning, waning), " 0.0 0.0 0.0")
        write(line, output_path, append=TRUE)
        
    }

    if(M3){
        line <- paste0("E_rl ", sprintf("%d ", target), target_group, " 1.0 1.0 1.0 ", sprintf("%f %f %f", cov, cov, cov), sprintf(" %.2f %.2f %.2f", waning, waning, waning), " 0.0 0.0 0.0")
        write(line, output_path, append=TRUE)
        
    }

    if(M4){
        line <- paste0("E_p_AON ", sprintf("%d ", target), target_group, " 1.0 1.0 1.0 ", sprintf("%f %f %f", cov, cov, cov), sprintf(" %.2f %.2f %.2f", waning, waning, waning), " 0.0 0.0 0.0")
        write(line, output_path, append=TRUE)
        line <- paste0("E_q ", sprintf("%d ", target), target_group, " 1.0 1.0 1.0 ", sprintf("%f %f %f", cov, cov, cov), sprintf(" %.2f %.2f %.2f", waning, waning, waning), " 0.0 0.0 0.0")
        write(line, output_path, append=TRUE)
        
    }

    if(M5){
        line <- paste0("E_p_AON ", sprintf("%d ", target), target_group, " 1.0 1.0 1.0 ", sprintf("%f %f %f", cov, cov, cov), sprintf(" %.2f %.2f %.2f", waning, waning, waning), " 0.0 0.0 0.0")
        write(line, output_path, append=TRUE)
        line <- paste0("E_rl ", sprintf("%d ", target), target_group, " 1.0 1.0 1.0 ", sprintf("%f %f %f", cov, cov, cov), sprintf(" %.2f %.2f %.2f", waning, waning, waning), " 0.0 0.0 0.0")
        write(line, output_path, append=TRUE)
        
    }

    if(M6){
        line <- paste0("E_q ", sprintf("%d ", target), target_group, " 1.0 1.0 1.0 ", sprintf("%f %f %f", cov, cov, cov), sprintf(" %.2f %.2f %.2f", waning, waning, waning), " 0.0 0.0 0.0")
        write(line, output_path, append=TRUE)
        line <- paste0("E_rl ", sprintf("%d ", target), target_group, " 1.0 1.0 1.0 ", sprintf("%f %f %f", cov, cov, cov), sprintf(" %.2f %.2f %.2f", waning, waning, waning), " 0.0 0.0 0.0")
        write(line, output_path, append=TRUE)
        
    }

    if(M7){
        line <- paste0("E_p_AON ", sprintf("%d ", target), target_group, " 1.0 1.0 1.0 ", sprintf("%f %f %f", cov, cov, cov), sprintf(" %.2f %.2f %.2f", waning, waning, waning), " 0.0 0.0 0.0")
        write(line, output_path, append=TRUE)
        line <- paste0("E_q ", sprintf("%d ", target), target_group, " 1.0 1.0 1.0 ", sprintf("%f %f %f", cov, cov, cov), sprintf(" %.2f %.2f %.2f", waning, waning, waning), " 0.0 0.0 0.0")
        write(line, output_path, append=TRUE)
        line <- paste0("E_rl ", sprintf("%d ", target), target_group, " 1.0 1.0 1.0 ", sprintf("%f %f %f", cov, cov, cov), sprintf(" %.2f %.2f %.2f", waning, waning, waning), " 0.0 0.0 0.0")
        write(line, output_path, append=TRUE)
        
    }
}

print(median(cov_save))




    
    
    



