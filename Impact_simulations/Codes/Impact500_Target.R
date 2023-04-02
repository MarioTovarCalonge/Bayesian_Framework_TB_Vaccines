#!/usr/bin/Rscript

args<-commandArgs(TRUE)
country <- args[1]
model <- args[2]
target <- as.integer(args[3])


available <- c(1, 0, 1, 2, 2, 3, 2, 1, 1, 1, 2, 1, 2, 1, 1, 2, 2, 2, 0, 2, 1, 2, 2, 2)
names(available) <- c('Bangladesh', 'Canada', 'China', 'Dem_Congo', 'Democratic_Republic_of_the_Congo', 'Ecuador', 'Ethiopia', 'India', 'Indonesia', 'Madagascar', 'Morocco', 'Myanmar', 'Nigeria', 'Pakistan', 'Philippines', 'Somalia', 'South_Africa', 'Tanzania', 'Ukraine', 'United_Republic_of_Tanzania', 'Viet_Nam', 'Uganda', 'Kenya', 'Zambia')
if(!(country %in% names(available))){
    print(paste0(country, ' is not available in the database. Raising error...'))
    stop()
}

#--------------------------------------------------------Functions-------------------------------------------------------#
get_vaccine_info <- function(path_to_profile){
    data <- read.table(path_to_profile, header=T, skip=2)
    transition_names <- c("E_rl ", "E_p ", "E_rs ", "E_rdef ", "E_rn ", "E_q ")
    get_transitions <- as.vector(unlist(data[,1]))
    age = data[1, ]$age
    target_type <- as.vector(data[1, ]$S_L_R)
    eff <- data[1, ]$efficacy
    cov <- data[1, ]$coverage
    waning <- data[1, ]$waning
    blocking <- data[1, ]$blocking
    transitions = get_transitions
    values <- c("age"= age, "efficacy" = eff, "coverage" = cov, "waning" = waning, "blocking" = blocking)
    ret <- list(transitions, target_type, values)
    return(ret)
}

get_vaccine_info2 <- function(mod){
    data <- read.table(paste0("profiles/M", mod, "profile.txt"), header=T, skip=2)
    transition_names <- c("E_rl ", "E_p ", "E_rs ", "E_rdef ", "E_rn ", "E_q ")
    get_transitions <- as.vector(unlist(data[,1]))
    age = data[1, ]$age
    target_type <- as.vector(data[1, ]$S_L_R)
    eff <- data[1, ]$efficacy
    cov <- data[1, ]$coverage
    waning <- data[1, ]$waning
    blocking <- data[1, ]$blocking
    transitions = get_transitions
    values <- c("age"= age, "efficacy" = eff, "coverage" = cov, "waning" = waning, "blocking" = blocking)
    ret <- list(transitions, target_type, values)
    return(ret)
}





#------------------------------------------------------------------------------------------------------------------------------------#

#----------------------------------------Block for calculate impacts simulation------------------------------------------#

scaling <- 1000000
ind_inc_D <- 19 + c(3, 4, 5, 6, 7, 8, 23, 24, 25, 26, 27, 28, 61, 62, 63, 56, 70, 71)
ind_mort <- 19 + c(36, 37, 38, 57, 58, 59, 60)
AG_ini <- 1
AG_fin <- 15
auxPath <- path.expand("~")

#Baseline
vac_path <- paste0("../Inputs/General_Inputs/perfil.txt")
info <- get_vaccine_info(vac_path)
#info <- get_vaccine_info2(model)
aux_tran <- unlist(info[1])
target_type <- unlist(info[2])
aux_info <- unlist(info[3])

line_tran <- target_type
for(h in 1:length(aux_tran)){
    aux_paste <- paste0(line_tran, "_", aux_tran[h])
    line_tran <- aux_paste
}


pathout_vaccine <- paste0(auxPath, "/TB_SIM_OUT_Bien_KCaR_TrNorm/", country, "/Vaccines/Vac_", sprintf("T_%d-E_%.3f-C_%.3f-W_%.2f", as.integer(aux_info["age"]), aux_info["efficacy"], aux_info["coverage"], aux_info["waning"] ),
                    "/", line_tran, "/Baseline/")

it0dat <- read.table(paste0(pathout_vaccine, "Output_vac.txt"), header=TRUE)

i <- it0dat$inc_annual
imp <- it0dat$inc_acum
i_vac <- it0dat$inc_cov_annual
imp_vac <- it0dat$inc_cov_acum

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

    pathout_vaccine <- paste0(auxPath, "/TB_SIM_OUT_Bien_KCaR_TrNorm/", country, "/Vaccines/Vac_", sprintf("T_%d-E_%.3f-C_%.3f-W_%.2f", as.integer(aux_info["age"]), aux_info["efficacy"], aux_info["coverage"], aux_info["waning"] ),
                    "/", line_tran, "/Iters/")
                    
    itIndexDat <- read.table(paste0(pathout_vaccine, "Output_vac_It_", sprintf("%03d", iter), ".txt"), header=TRUE)
    
    i <- itIndexDat$inc_annual
    imp <- itIndexDat$inc_acum
    i_vac <- itIndexDat$inc_cov_annual
    imp_vac <- itIndexDat$inc_cov_acum
 
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
system("mkdir -p ../Impact_Outputs/iteration_impacts/")
df_ii <- data.frame(ID=1:length(saved_i_con), imp_con = saved_imp_con, imp_vac = saved_imp_vac, i_con = saved_i_con, i_vac = saved_i_vac)
write.table(df_ii, paste0("../Impact_Outputs/iteration_impacts/", country, "_", "M", model, "_T_", target, "_dataimp.txt"), quote=FALSE, row.names=FALSE)
#############################################################

impact_low <- quantile(saved_impact, .025)    
impact_up <- quantile(saved_impact, .975)
impact_median <- quantile(saved_impact, .5)
    
IRR_low <- quantile(saved_IRR, .025)    
IRR_up <- quantile(saved_IRR, .975)
IRR_median <- quantile(saved_IRR, .5)
    

system("mkdir -p ../Impact_Outputs/Impacts_models_per_country/")
if(!file.exists(paste0("../Impact_Outputs/Impacts_models_per_country/impact_", country, "_T_", target, ".txt"))){
    line <- sprintf("1 %f %f %f %f %f %f", impact_median, impact_low, impact_up, IRR_median, IRR_low, IRR_up)
    write(line, paste0("../Impact_Outputs/Impacts_models_per_country/impact_", country, "_T_", target, ".txt"), append=T)
} else{
    aux <- read.table(paste0("../Impact_Outputs/Impacts_models_per_country/impact_", country, "_T_", target, ".txt"), header=FALSE)
    l <- length(aux$V1)
    line <- sprintf("%d %f %f %f %f %f %f", l+1, impact_median, impact_low, impact_up, IRR_median, IRR_low, IRR_up)
    write(line, paste0("../Impact_Outputs/Impacts_models_per_country/impact_", country, "_T_", target, ".txt"), append=T)
}
}



    
    
    

