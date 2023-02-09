#!/usr/bin/Rscript

args<-commandArgs(TRUE)
country <- args[1]
model <- args[2]


available <- c(1, 0, 1, 2, 2, 3, 2, 1, 1, 1, 2, 1, 2, 1, 1, 2, 2, 2, 0, 2, 1, 2, 2, 2)
names(available) <- c('Bangladesh', 'Canada', 'China', 'Dem_Congo', 'Democratic_Republic_of_the_Congo', 'Ecuador', 'Ethiopia', 'India', 'Indonesia', 'Madagascar', 'Morocco', 'Myanmar', 'Nigeria', 'Pakistan', 'Philippines', 'Somalia', 'South_Africa', 'Tanzania', 'Ukraine', 'United_Republic_of_Tanzania', 'Viet_Nam', 'Uganda', 'Kenya', 'Zambia')
if(!(country %in% names(available))){
    print(paste0(country, ' is not available in the database. Raising error...'))
    stop()
}

#--------------------------------------------------------Functions-------------------------------------------------------#
readOutputPop <- function(obj, a_ini, a_fin=NULL){ 
    reserv <- 19
    AG <- 15
    sim_years <- length(obj$V1)/3
    if (!is.null(a_fin)) {
        if(a_ini!=a_fin){
            if(a_ini<a_fin){
                ids <- seq(1+a_ini, AG*reserv+1, 15)
                for(j in (a_ini+1):a_fin){
                    aux <- seq(1 + j, AG*reserv+1, 15)
                    ids <- c(ids, aux)
                }
            } else{
                print("Error: a_ini should be smaller than a_fin\n")
                return(0)
            }
        } else{
            ids <- seq(1+a_ini, AG*reserv+1, 15)
        }
    } else{
        ids <- seq(1+a_ini, AG*reserv+1, 15)
    }
    p <- rowSums(obj[,ids])
    temp <- p[1:sim_years] + p[(sim_years+1):(2*sim_years)] + p[(2*sim_years+1):(3*sim_years)]
    return(temp)
}

readOutputFlux <- function(obj, a_ini, a_fin=NULL, ind){ 
    AG <- 15
    sim_years <- length(obj$V1)/3
    if (!is.null(a_fin)) { 
        if(a_ini!=a_fin){
            if(a_ini<a_fin){
                ind <- 1 + AG*ind
                conc <- seq(ind[1] + a_ini, ind[1] + a_fin)
                for(i in 2:length(ind)){
                    aux <- seq(ind[i] + a_ini, ind[i] + a_fin)
                    conc <- c(conc, aux)
                }
                ids <- conc
            } else{
                print("Error: a_ini should be smaller than a_fin\n")
                return(0)
            }
        } else{
            ids <- 1 + a_ini + AG*ind
        }
    } else{
        ids <- 1 + a_ini + AG*ind
    }
    p <- rowSums(obj[,ids])
    temp <- p[1:sim_years] + p[(sim_years+1):(2*sim_years)] + p[(2*sim_years+1):(3*sim_years)]
    return(temp)
}

calc <- function(v){
    temp <- numeric(length(v))
    for(i in 1:length(v)){
        temp[i] <- sum(v[1:i])
    }
    return(temp)
}


#------------------------------------------------------------------------------------------------------------------------------------#

#----------------------------------------Block for calculate impacts simulation------------------------------------------#

scaling <- 1000000
ind_inc_D <- 19 + c(3, 4, 5, 6, 7, 8, 23, 24, 25, 26, 27, 28, 61, 62, 63, 56, 70, 71)
ind_mort <- 19 + c(36, 37, 38, 57, 58, 59, 60)
AG_ini <- 1
AG_fin <- 15
auxPath <- path.expand("~")

path_base <- paste0(auxPath, "/TB_SIM_OUT/", country, '/Object_Folder/Baseline/Output_object_1.txt')
base <- read.table(path_base, header=F)
i <- readOutputFlux(base, AG_ini, AG_fin, ind = ind_inc_D)
imp <- calc(i)

path_vac <- paste0(auxPath, "/TB_SIM_OUT/", country, '/Object_Folder/Baseline/Output_object_2.txt')
vac <- read.table(path_vac, header=F)
i_vac <- readOutputFlux(vac, AG_ini, AG_fin, ind = ind_inc_D)
imp_vac <- calc(i_vac)

Impact <- (imp[length(imp)] - imp_vac[length(imp)])/imp[length(imp)]*100
IRR <- (i[length(i)] - i_vac[length(i)])/i[length(i)]*100

miniter <- 1
maxiter <- 500
saved_impact <- numeric( maxiter - miniter )
saved_IRR <- numeric( maxiter - miniter )
#########################NEW############################
saved_imp_con <- numeric( maxiter - miniter )
saved_imp_vac <- numeric( maxiter - miniter )
saved_i_con <- numeric( maxiter - miniter )
saved_i_vac <- numeric( maxiter - miniter )

if(TRUE){
print("Reading outputs\n")
pb <- txtProgressBar(min = 1, max = (maxiter - miniter), style = 3)
for(iter in miniter:maxiter){
    index = iter - miniter + 1
    
    path_base <- paste0(auxPath, "/TB_SIM_OUT/", country, '/Object_Folder/Iter_', sprintf("%03d", iter), "/Output_object_1.txt")
    base <- read.table(path_base, header=F)
    i <- readOutputFlux(base, AG_ini, AG_fin, ind = ind_inc_D)
    imp <- calc(i)
    
    path_vac <- paste0(auxPath, "/TB_SIM_OUT/", country, '/Object_Folder/Iter_', sprintf("%03d", iter), "/Output_object_2.txt")
    vac <- read.table(path_vac, header=F)
    i_vac <- readOutputFlux(vac, AG_ini, AG_fin, ind = ind_inc_D)
    imp_vac <- calc(i_vac)
 
    saved_impact[index] = (imp[length(imp)] - imp_vac[length(imp)])/imp[length(imp)]*100
    saved_IRR[index] = (i[length(i)] - i_vac[length(i)])/i[length(i)]*100
    
    ######################NEW##########################
    saved_imp_con[index] = imp[length(imp)]
    saved_imp_vac[index] = imp_vac[length(imp)]
    saved_i_con[index] = i[length(i)]
    saved_i_vac[index] = i_vac[length(i)]
    
    setTxtProgressBar(pb, index)
}
close(pb)


#############################################################
df_ii <- data.frame(ID=1:length(saved_i_con), imp_con = saved_imp_con, imp_vac = saved_imp_vac, i_con = saved_i_con, i_vac = saved_i_vac)
write.table(df_ii, paste0("../Impact_Outputs/iteration_impacts/", country, "_", "M", model, "_dataimp.txt"), quote=FALSE, row.names=FALSE)
#############################################################

impact_low <- quantile(saved_impact, .025)    
impact_up <- quantile(saved_impact, .975)
impact_median <- quantile(saved_impact, .5)
    
IRR_low <- quantile(saved_IRR, .025)    
IRR_up <- quantile(saved_IRR, .975)
IRR_median <- quantile(saved_IRR, .5)
    


if(!file.exists(paste0("../Impact_Outputs/Impacts_models_per_country/impact_", country, ".txt"))){
    line <- sprintf("1 %f %f %f %f %f %f", impact_median, impact_low, impact_up, IRR_median, IRR_low, IRR_up)
    write(line, paste0("../Impact_Outputs/Impacts_models_per_country/impact_", country, ".txt"), append=T)
} else{
    aux <- read.table(paste0("../Impact_Outputs/Impacts_models_per_country/impact_", country, ".txt"), header=FALSE)
    l <- length(aux$V1)
    line <- sprintf("%d %f %f %f %f %f %f", l+1, impact_median, impact_low, impact_up, IRR_median, IRR_low, IRR_up)
    write(line, paste0("../Impact_Outputs/Impacts_models_per_country/impact_", country, ".txt"), append=T)
}
}



    
    
    

