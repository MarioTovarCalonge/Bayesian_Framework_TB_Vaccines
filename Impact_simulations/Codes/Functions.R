#!/usr/bin/Rscript

#' Function used for creating all the stuff neccesary for the code to run.
#'
#' @param ctr A string containing the name of the country in which the simulation is launched.
#' @param itr An integer indicating the number of the error calculation iteration.
#' @return This function does return the path to the master file created, which allows c functions to locate the inputs contained on the package database.
#' @keywords masterfile workspace
#' @examples
#' f("China", 0)
set_TBsim_workspace <- function(ctr, itr){
    available <- c(1, 0, 1, 2, 2, 3, 2, 1, 1, 1, 2, 1, 2, 1, 1, 2, 2, 2, 0, 2, 1, 2, 2, 2)
    names(available) <- c('Bangladesh', 'Canada', 'China', 'Dem_Congo', 'Democratic_Republic_of_the_Congo', 'Ecuador', 'Ethiopia', 'India', 'Indonesia', 'Madagascar', 'Morocco', 'Myanmar', 'Nigeria', 'Pakistan', 'Philippines', 'Somalia', 'South_Africa', 'Tanzania', 'Ukraine', 'United_Republic_of_Tanzania', 'Viet_Nam', 'Uganda', 'Zambia', 'Kenya')
    if(!(ctr %in% names(available))){
        print(paste0(ctr, ' is not available in the database. Raising error...'))
        stop()
    }
    
    pathToMF <- create_maste_file(ctr, itr, available[ctr])
    
    if(itr==0){
        output_dir <- paste0("../Outputs/", ctr, "/")
        if (!dir.exists(output_dir)){
            dir.create(paste0(output_dir, "Fitted_inputs/Initial_Conditions/"), recursive = TRUE)
            dir.create(paste0(output_dir, "Iter/"), recursive = TRUE)
            dir.create(paste0(output_dir, "Vaccines/"), recursive = TRUE)
            create_initial_fit_param(ctr, 0)
            
            ini_cond_file <- paste0("../Outputs/", ctr, "/Fitted_inputs/Initial_Conditions/initial_conditions.txt")
            if(!file.exists(ini_cond_file)){
                file.copy("../Inputs/General_Inputs/Initial_Conditions/initial_conditions_Base.txt", ini_cond_file)
            }
            
        } else {
            print("Master Dir already exists!")
        }
    }
        
    if(itr>0){
        output_dir <- paste0("../Outputs/", ctr, "/Iter/", ctr, "_", sprintf("%03d", itr), "/Fitted_inputs/Initial_Conditions/")
        if (!dir.exists(output_dir)){
            dir.create(output_dir, recursive = TRUE)
            destfile <- paste0("../Outputs/", ctr, "/Fitted_inputs/Parameters_to_fit.txt")
            if(!file.exists(destfile)){
                create_initial_fit_param(ctr, itr)
            } else{
                file.copy(destfile, paste0("../Outputs/", ctr, "/Iter/", ctr, "_", sprintf("%03d", itr), "/Fitted_inputs/Parameters_to_fit.txt"))
            }
            system(paste0("cp ../Outputs/", ctr, "/Fitted_inputs/Initial_Conditions/initial_conditions.txt ../Outputs/", ctr, "/Iter/", ctr, "_", sprintf("%03d", itr), "/Fitted_inputs/Initial_Conditions/initial_conditions.txt"))
        } else {
            print("SubDir already exists!")
        }
        
    }
    return(pathToMF)
}

#' A Function that build the master file for c codes based on the relative paths to the inputs contained in package database.
#'
#' @param ctr A string containing the name of the country in which the simulation is launched.
#' @param itr An integer indicating the number of the error calculation iteration.
#' @param reg An integer indicating the region of the country (0 Europe, 1 Asia, 2 Africa).
#' @return This function does return the path to the master file created.
#' @keywords masterfile
create_maste_file <- function(ctr, itr, reg){
    path_MF <- tempdir()
    print(path_MF)
    #path_MF <- paste0("./")
    temp_path <- paste0(path_MF, "mf_", sprintf("%03d", itr), ".txt")

    if(itr==0){
	    write(paste0("SIMULATION_NAME ", ctr), temp_path)
        write(paste0("COUNTRY ", ctr, "\n"), temp_path, append=T)
	    
	    write(paste0("EXECUTION_MODES"), temp_path, append=T)
	    write(sprintf("dynamic_type %d", 2), temp_path, append=T)
        write(sprintf("contact_type %d", 1), temp_path, append=T)

	    write(paste0("INPUTS: \n"), temp_path, append=T)

	    path_input_dict <- paste0("../Inputs/General_Inputs/dictionary.txt")
        write(paste0("states_dictionary ", path_input_dict), temp_path, append=T)

        path_global_par <- paste0("../Inputs/General_Inputs/global_parameters_mc_update.txt")
        write(paste0("global_parameters ", path_global_par), temp_path, append=T)

        path_reg_par <- paste0("../Inputs/Input_countries/", ctr, "/regional_parameters_literature.txt")
	    write(paste0("regional_parameters ", path_reg_par), temp_path, append=T)

        if(reg==0){ aux <- sprintf("../Inputs/Contacts/Polymod_contacts_pi.txt") }
        if(reg==1){ aux <- sprintf("../Inputs/Contacts/Asia_contacts_pi.txt") }
        if(reg==2){ aux <- sprintf("../Inputs/Contacts/Africa_contacts_pi.txt") }
	    write(sprintf("contacts %s \n", aux), temp_path, append=T);

        path_demo<- paste0("../Inputs/Input_countries/", ctr,"/Demographic_pyramid.txt")
        write(paste0("demography ", path_demo), temp_path, append=T)

	    path_window <- paste0("../Inputs/General_Inputs/window_2100.txt")
        write(paste0("window ", path_window), temp_path, append=T)

    } else{
	    
        write(paste0("SIMULATION_NAME ", ctr, "_", sprintf("%03d", itr)), temp_path)
        write(paste0("COUNTRY ", ctr, "\n"), temp_path, append=T)
	    
        write(paste0("EXECUTION_MODES"), temp_path, append=T)
	    write(sprintf("dynamic_type %d", 2), temp_path, append=T)
        write(sprintf("contact_type %d", 1), temp_path, append=T)

	    write(paste0("INPUTS: \n"), temp_path, append=T)

	    path_input_dict <- paste0("../Inputs/General_Inputs/dictionary.txt")
        write(paste0("states_dictionary ", path_input_dict), temp_path, append=T)

        path_global_par <- paste0("../Inputs_Iter/", ctr, "/Global/global_parameters_mc_update_", sprintf("%03d", itr),".txt")
	    write(paste0("global_parameters ", path_global_par), temp_path, append=T)

        path_reg_par <- paste0("../Inputs_Iter/", ctr, "/Regional/regional_parameters_literature_", sprintf("%03d", itr),".txt")
	    write(paste0("regional_parameters ", path_reg_par), temp_path, append=T)
	    
        if(reg==0){ aux <- paste0("../Inputs_Iter/", ctr,"/Contacts/Polymod_contacts_pi_", sprintf("%03d", itr), ".txt") }
        if(reg==1){ aux <- paste0("../Inputs_Iter/", ctr,"/Contacts/Asia_contacts_pi_", sprintf("%03d", itr), ".txt") }
        if(reg==2){ aux <- paste0("../Inputs_Iter/", ctr,"/Contacts/Africa_contacts_pi_", sprintf("%03d", itr), ".txt") }
	    write(sprintf("contacts %s \n", aux), temp_path, append=T);

        path_demo <- paste0("../Inputs_Iter/", ctr, "/Demo/Demographic_pyramid_", sprintf("%03d", itr),".txt")
	    write(paste0("demography ", path_demo), temp_path, append=T)

	    path_window <- paste0("../Inputs/General_Inputs/window_2100.txt")
        write(paste0("window ", path_window), temp_path, append=T)	    

    }

    return(temp_path)
}

#' This function is used for creating a default file of parameters that are fitted by the code for minimise fitness.
#'
#' @param ctr A string containing the name of the country in which the simulation is launched.
#' @param itr An integer indicating the number of the error calculation iteration.
#' @return This function does not return anything.
#' @keywords Parameters
create_initial_fit_param <- function(ctr, itr){
    param <- c("d0", "beta0", "beta1", "d1", "dist_est")
    param_up <- c(12, 50, 20, 0.2, 0.5)
    param_low <- c(0, 0, -0.2, -0.2, -0.5)
    if(itr==0){
        param_file <- paste0("../Outputs/", ctr, "/Fitted_inputs/Parameters_to_fit.txt")
    } else{
        param_file <- paste0("../Outputs/", ctr, "/Iter/", ctr, "_", sprintf("%03d", itr), "/Fitted_inputs/Parameters_to_fit.txt")
    }
    write(paste0("entries = 5"), param_file)
    write(paste0("Parameter Index Type Amin Amax Value LB UB"), param_file, append=T)
    for(i in 0:4){
        num = runif(1, 0, 1)
        line <- paste0(param[i+1], sprintf(" %d Diag_force 0 120 %.2f %f %f", i, num, param_low[i+1], param_up[i+1]))
        write(line, param_file, append=T)
    }
}

#' This function is used for creating an user-defined vaccine profile.
#'
#' @param target An integer in the interval [0-14] pointing out the age group in which vaccination is desired.
#' @param target_group A string indicating if vaccine is applied to susceptibles, latent or recovered individuals.
#' @param eff A number that controls the efficacy of the vaccine [0-1]. Defaults to 1 which means perfect protection.
#' @param cov A number that controls the fraction of the targeted population that will be vaccinated [0-1]. Defaults to 1 which means total coverage.
#' @param waning Numeric input that indicate the desired level of annual waning (0.05) stands for a 5\% waning annualy.
#' @return This function does not return anything.
#' @keywords vaccine
create_vaccine_file <- function(target=12, target_group='SLR', eff=1.0, cov=1.0, waning=0.00){
    if (target < 0) stop("'target' must be >= 0")
    if (target > 15) stop("'target' must be <= 15")
    if (eff < 0) stop("'eff' must be >= 0")
    if (eff > 1) stop("'eff' must be <= 1.0")
    if (cov < 0) stop("'cov' must be >= 0")
    if (cov > 1) stop("'cov' must be <= 1.0")
    if (waning < 0) stop("'waning' must be >= 0")

    path_vaccine_file <- "../Inputs/General_Inputs/perfil.txt"
    transition_names <- c("E_rl ", "E_p ", "E_rs ", "E_rdef ", "E_rn ", "E_q ")

    line <- sprintf("target = %d", target)
    write(line, path_vaccine_file)
    line <- sprintf("waning_harris = %f", 10)
    write(line, path_vaccine_file, append=TRUE)
    line <- "transition age S_L_R efficacy e_low e_hi coverage c_low c_hi waning w_low w_hi blocking b_low b_hi"
    write(line, path_vaccine_file, append=TRUE)

    for(i in 1:length(transition_names)){
        line <- paste0(transition_names[i], sprintf("%d ", target), target_group, 
                    " 1.0 1.0 1.0 ", sprintf("%.2f %.2f %.2f", eff*cov, eff*cov, eff*cov), 
                    sprintf(" %.2f %.2f %.2f", waning, waning, waning), " 0.0 0.0 0.0")
        write(line, path_vaccine_file, append=TRUE)
    }
}

#' This function is used for creating new files for parameters and demography by taking the central values and the bounds as the 95\% CI of a normal distribution, and generating a new central value according to this distribution.
#'
#' @param country A string containing the name of the country in which the simulation is launched.
#' @param ini An integer indicating the number of iteration for which we will produce a new contact file for the error calculation procedure.
#' @param fin An integer indicating if more than one file should be generated for the error calculation procedure (i.e) if we want to create all the files corresponding to the itearations all together. It defaults to ini as this is usually done one by one.
#' @return This function does return nothing.
#' @keywords Error
create_recal_parameters <- function(country, ini, fin, owrite=TRUE){
    ap <- function(ag, Zm, mn, d_inf, d_sup){
        amax <- 14
        if(Zm<0){
            aux <- (Zm*(2*ag-amax)*d_sup)/(amax*1.96) + mn
        } else{
            aux <- (Zm*(2*ag-amax)*d_inf)/(amax*1.96) + mn
        }
        return(aux)
    }
    
    overwrite <- owrite

    if(.Platform$OS.type == "unix") {
        system(paste0("mkdir -p ../Inputs_Iter/", country, "/Demo/"))
        system(paste0("mkdir -p ../Inputs_Iter/", country, "/Burden/"))
        system(paste0("mkdir -p ../Inputs_Iter/", country, "/Global/"))
        system(paste0("mkdir -p ../Inputs_Iter/", country, "/Regional/"))
    } else {
        system(paste0("mkdir ../Inputs_Iter/", country, "/Demo/"))
        system(paste0("mkdir ../Inputs_Iter/", country, "/Burden/"))
        system(paste0("mkdir ../Inputs_Iter/", country, "/Global/"))
        system(paste0("mkdir ../Inputs_Iter/", country, "/Regional/"))
    }

    #------------------------------------Recal Demographic pyramids----------------------------------#
    file_in <- paste0('../Inputs/Input_countries/', country, '/Demographic_pyramid.txt')
    Z <- numeric(as.numeric(fin)-as.numeric(ini))

    for(a in ini:fin){
        file_out = paste0('../Inputs_Iter/', country, '/Demo/Demographic_pyramid_', sprintf("%03d", a), '.txt')
        demo <- read.table(file_in, header=F, skip=2)

        mean <- demo$V7
        sd <- (demo$V9 - demo$V8)/3.92
        Zs <- rnorm(1, 0, 1)
        Z[a - as.numeric(ini) + 1] <- Zs

        arg <- rep(seq(0, 14), length(mean)/15)

        demo$V7 <- ap(arg, Zs, mean, mean - demo$V8, demo$V9 - mean)

        if(file.exists(file_out)){
            if(overwrite){
                line <- sprintf("entries=%d", length(mean))
                write(line, file_out)
                line <- sprintf("amin    amax    emin    emax    tmin    tmax    value    interval")
                write(line, file_out, append=T)
                write.table(demo, file_out, row.names=FALSE, col.names=FALSE, quote=FALSE, append=T)
            } else{
                print("Not overwriting")
            }
        } else{
            line <- sprintf("entries=%d", length(mean))
            write(line, file_out)
            line <- sprintf("amin    amax    emin    emax    tmin    tmax    value    interval")
            write(line, file_out, append=T)
            write.table(demo, file_out, row.names=FALSE, col.names=FALSE, quote=FALSE, append=T)
        }
    }

    #------------------------------------------Recal TB burden---------------------------------------#
    file_in <- paste0('../Inputs/Input_countries/', country, '/TB_burden.txt')
    for(a in ini:fin){
        Zs <- rnorm(1, 0, 1)
        
        #WHO INCIDENCE AND MORTALITY
        file_out = paste0('../Inputs_Iter/', country, '/Burden/TB_burden_', sprintf("%03d", a), '.txt')
        burden <- read.table(file_in, header=T)
        
        sd_inc <- (burden$Inc_hi - burden$Inc_lo)/3.92
        mean_inc <- burden$Inc
        burden$Inc <- mean_inc + Zs*sd_inc

        sd_mort <- (burden$Mort_hi - burden$Mort_lo)/3.92
        mean_mort <- burden$Mort
        burden$Mort <- mean_mort + Zs*sd_mort
        
        #Check if file exists and should ve overwriten or not.
        if(file.exists(file_out)){
            if(overwrite){
                write.table(burden, file_out, row.names=FALSE, quote=FALSE)
            } else{
                print("Not overwriting")
            }
        } else{
            write.table(burden, file_out, row.names=FALSE, quote=FALSE)
        }
    }

    #-------------------------------------------Recal Global----------------------------------------#
    filein <- paste0('../Inputs/General_Inputs/global_parameters_mc_update.txt')
    for(a in ini:fin){
        file_out = paste0('../Inputs_Iter/', country, '/Global/global_parameters_mc_update_', sprintf("%03d", a), '.txt')

        par <- read.table(filein, header=F, skip=2)
        mean <- par$V7
        sd <- (par$V9 - par$V8)/3.92
        
        new <- NNegD(mean, par$V8, par$V9)
        par$V7 <- new

        #Check if file exists and should ve overwriten or not.
        if(file.exists(file_out)){
            if(overwrite){
                line <- sprintf("entries = %d", length(mean))
                write(line, file_out)
                line <- sprintf("description index      amin      amax      emin      emax      value      interval")
                write(line, file_out, append=T)
                write.table(par, file_out, row.names=FALSE, col.names=FALSE, quote=FALSE, append=T)
            } else{
                print("Not overwriting")
            }
        } else{
            line <- sprintf("entries = %d", length(mean))
            write(line, file_out)
            line <- sprintf("description index      amin      amax      emin      emax      value      interval")
            write(line, file_out, append=T)
            write.table(par, file_out, row.names=FALSE, col.names=FALSE, quote=FALSE, append=T)
        }
    }

    #------------------------------------------Recal Regional---------------------------------------#
    filein <- paste0('../Inputs/Input_countries/', country, '/regional_parameters_literature.txt')
    for(a in ini:fin){
        file_out = paste0('../Inputs_Iter/', country, '/Regional/regional_parameters_literature_', sprintf("%03d", a), '.txt')

        reg <- read.table(filein, header=F, skip=2)
        mean <- reg$V7
        sd <- (reg$V9 - reg$V8)/3.92
        
        new <- NNegD(mean, reg$V8, reg$V9) #Log-normal
        reg$V7 <- new
        
        #Check if file exists and should ve overwriten or not.
        if(file.exists(file_out)){
            if(overwrite){
                line <- sprintf("entries = %d", length(mean))
                write(line, file_out)
                line <- sprintf("description index      amin      amax      emin      emax      value      interval")
                write(line, file_out, append=T)
                write.table(reg, file_out, row.names=FALSE, col.names=FALSE, quote=FALSE, append=T)
            } else{
                print("Not overwriting")
            }
        } else{
            line <- sprintf("entries = %d", length(mean))
            write(line, file_out)
            line <- sprintf("description index      amin      amax      emin      emax      value      interval")
            write(line, file_out, append=T)
            write.table(reg, file_out, row.names=FALSE, col.names=FALSE, quote=FALSE, append=T)
        }
    }
}


#' This function is used for creating a new contact file by taking the central values and the bounds as the 95\% CI of a normal distribution, and generating a new central value according to this distribution.
#'
#' @param country A string containing the name of the country in which the simulation is launched.
#' @param ini An integer indicating the number of iteration for which we will produce a new contact file for the error calculation procedure.
#' @param fin An integer indicating if more than one file should be generated for the error calculation procedure (i.e) if we want to create all the files corresponding to the itearations all together. It defaults to ini as this is usually done one by one.
#' @param reg An integer indicating the region of the country (0 Europe, 1 Asia, 2 Africa).
#' @return This function does return nothing.
#' @keywords Contacts
create_recal_contacts <- function(country, ini, fin=ini, reg){
    if(reg==0){
        doc <- '../Inputs/Contacts/Polymod_contacts_pi.txt'
    } else if(reg==1){
        doc <- '../Inputs/Contacts/Asia_contacts_pi.txt'
    } else if(reg==2){
        doc <- '../Inputs/Contacts/Africa_contacts_pi.txt'
    } else if(reg==3){
        doc <- '../Inputs/Contacts/LATAM_contacts_pi.txt'
    }

    A <- data.matrix(read.table(doc, nrows=15))
    B <- data.matrix(read.table(doc, skip=15, nrows=15))
    C <- data.matrix(read.table(doc, skip=31, nrows=15))

    M <- matrix(0, length(A[,1]), length(A[,1]))
    
    options(digits=16)

    if(.Platform$OS.type == "unix") {
        system(paste0("mkdir -p ../Inputs_Iter/", country, "/Contacts/"))
    } else {
        system(paste0("mkdir ../Inputs_Iter/", country, "/Contacts/"))
    }

    for(a in ini:fin){
        for(i in 1:length(A[,1])){
            for(j in 1:length(A[,1])){
                mn <- A[i,j]
                sd_pos <- (C[i,j]-mn)/1.96
                sd_neg  <- (mn-B[i,j])/1.96
                
                Z <- rnorm(1, 0, 1)
                while(abs(Z)>1.96){
                    Z <- rnorm(1, 0, 1)
                }
                if(Z>0){
                    r <- mn + Z*sd_pos
                } else{
                    r <- mn + Z*sd_neg
                }
                M[i, j] <- r
            }
        }

        if(reg==1){
            output = paste0('../Inputs_Iter/', country, '/Contacts/Asia_contacts_pi_', sprintf("%03d", a), '.txt')
        } else if(reg==2){
            output = paste0('../Inputs_Iter/', country, '/Contacts/Africa_contacts_pi_', sprintf("%03d", a), '.txt')
        } else if(reg==0){
            output = paste0('../Inputs_Iter/', country, '/Contacts/Polymod_contacts_pi_', sprintf("%03d", a), '.txt')
        } else if(reg==3){
            output = paste0('../Inputs_Iter/', country, '/Contacts/LATAM_contacts_pi_', sprintf("%03d", a), '.txt')
        }
        
        write.table(M, file=output, row.names=FALSE, col.names=FALSE, quote=FALSE)
        write("", output, append=T)
        write.table(B, file=output, row.names=FALSE, col.names=FALSE, quote=FALSE, append=T)
        write("", output, append=T)
        write.table(C, file=output, row.names=FALSE, col.names=FALSE, quote=FALSE, append=T)
        write("", output, append=T)
    }
}

NNegD <- function(mu, alpha, beta){
    c = beta > alpha
    index <- which(F==c)
    for(i in index){
        aux = beta[i]
        beta[i] = alpha[i]
        alpha[i] = aux
    }
        
    Z <- rnorm(length(mu), 0, 1)
    while(sum(abs(Z)>1.96)>0){
        Z <- rnorm(length(mu), 0, 1)
    }
    sigma_pos <- (beta-mu)/1.96
    sigma_neg <- (mu-alpha)/1.96
        
    i_pos <- which(Z>0)
    i_neg <- which(Z<0)
    r <- numeric(length(mu))
    for(i in i_pos){
        r[i] = mu[i] + Z[i]*sigma_pos[i]
    }
    for(i in i_neg){
        r[i] = mu[i] + Z[i]*sigma_neg[i]
    }
    return(r)
}

#' A Function for fitting demographic pyramids included in the database to 9th grade polynomials whose coefficients are used in both fitter and runner for obtaining demographic evolution.
#'
#' @param demo_file A string containing the path to the pyramid file. No default.
#' @param ctr A string containing the name of the country whose demographic evolution will be fitted. No default.
#' @param itr An integer pointing the iteration of the simulation. Defaults to 0, which is the iteration of the baseline.
#' @param y_min_fit An integer for the first year of the fit. Defaults to 2000.
#' @param y_max_fit An integer for the last year of the fit. Defaults to 2018, which is the last year available of the TB Burden in package database.
#' @param y_min_run An integer for the first year of the run. Defaults to 2000.
#' @param y_max_run An integer for the last year of the run. Defaults to 2100, which is the last year of UN database demographic predictions.
#' @return This function does not return anything.
#' @keywords demography, pyramids, demographic evolution
demography_fit <- function(demo_file, ctr, itr, y_min_fit=2000, y_max_fit=2018, y_min_run=2000, y_max_run=2100){
    options(scipen=999)
    if(itr==0){
        fitting_polynomial_file = paste0("../Outputs/", ctr, "/Fitted_inputs/Pyr_coef_fit.txt")
        running_polynomial_file = paste0("../Outputs/", ctr, "/Fitted_inputs/Pyr_coef_run.txt")
    } else{
        fitting_polynomial_file = paste0("../Outputs/", ctr, "/Iter/", ctr, "_", sprintf("%03d", itr), "/Fitted_inputs/Pyr_coef_fit.txt")
        running_polynomial_file = paste0("../Outputs/", ctr, "/Iter/", ctr, "_", sprintf("%03d", itr), "/Fitted_inputs/Pyr_coef_run.txt")
    }

    aux <- read.table(demo_file, header=F, skip=2)
    AG <- 15
    b <- list()
    for(i in 1:AG){
        b[[i]] <- as.vector(unlist(aux$V7[aux$V1==(5*(i-1))]))
    }

    est <- matrix(0, 15, 10)
    est_r <- matrix(0, 15, 10)
    x <- seq(0, length(b[[1]]) )
    for (k in 1:15){       
        y <- as.vector(unlist(b[[k]]))
        df <- data.frame(t = x[1:(y_max_fit - y_min_fit + 1)], P = y[(y_min_fit - 1950 + 1):(y_max_fit - 1950 + 1)])
        df2 <- data.frame(t = x[1:(y_max_run - y_min_run + 1)], P = y[(y_min_run - 1950 + 1):(y_max_run - 1950+1)])

        res <- lm(formula = P ~ poly(t, 9, raw = TRUE), data = df)
        res_2 <- lm(formula = P ~ poly(t, 9, raw = TRUE), data = df2)

        for (i in 1:10){
            est[k, i] <- as.numeric(res$coef[i])
            est_r[k, i] <- as.numeric(res_2$coef[i])
        }
        
    }

    write(sprintf("%d %d %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f", 0, 0, est[1, 1], est[1, 2], est[1, 3], est[1, 4], est[1, 5], est[1, 6], est[1, 7], est[1, 8], est[1, 9], est[1, 10]), fitting_polynomial_file)
    write(sprintf("%d %d %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f", 0, 0, est_r[1, 1], est_r[1, 2], est_r[1, 3], est_r[1, 4], est_r[1, 5], est_r[1, 6], est_r[1, 7], est_r[1, 8], est_r[1, 9], est_r[1, 10]), running_polynomial_file)
    for (k in 1:14){
        write(sprintf("%d %d %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f", k, 0, est[k+1, 1], est[k+1, 2], est[k+1, 3], est[k+1, 4], est[k+1, 5], est[k+1, 6], est[k+1, 7], est[k+1, 8], est[k+1, 9], est[k+1, 10]), fitting_polynomial_file, append=TRUE)
        write(sprintf("%d %d %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f", k, 0, est_r[k+1, 1], est_r[k+1, 2], est_r[k+1, 3], est_r[k+1, 4], est_r[k+1, 5], est_r[k+1, 6], est_r[k+1, 7], est_r[k+1, 8], est_r[k+1, 9], est_r[k+1, 10]), running_polynomial_file, append=TRUE)
    }
}

#' This function is used for writting into a temporal file the fitted parameters given by the LM algorithm. This temporal file will be used when running the simulation.
#'
#' @param x A vector containing the values of the parameters that are being calibrated by the LM algorithm.
#' @param ctr A string containing the name of the country in which the simulation is launched.
#' @param itr An integer indicating the number of the error calculation iteration.
#' @return This function does not return anything.
rewritte_parameters_to_fit <- function(x, ctr, itr=0){
    if(itr==0){
        par_dest <- paste0("../Outputs/", ctr, "/Fitted_inputs/Parameters_to_fit.txt")
    } else{
        par_dest <- paste0("../Outputs/", ctr, '/Iter/', ctr, '_', sprintf("%03d", itr), "/Fitted_inputs/Parameters_to_fit.txt")
    }

    df <- read.table(par_dest, skip=1, header=T)
    df$Value <- x 

    line <- sprintf("entries = %d", length(x))
    write(line, par_dest)
    line <- "Parameter Index Type Amin Amax Value LB UB"
    write(line, par_dest, append=T)
    write.table(format(df, digits=22), par_dest, col.names=F, row.names=FALSE, quote=FALSE, append=TRUE)
}

#' This function is an interface for call the fitting part of the code within R. It allows user defined fitness functions to be passed as arguments to the code. It is used inside the fitter and fitter_nlfb funtions and should not be modified by the user. 
#'
#' @param FUN A function for calculate fitness. It can be user-defined and should return a vector of residuals after evaluation of fitness criterion. Default is set to WHO data fitness.  
#' @param obj An Output_type object produced by fitter. FUN is supposed to calculate residuals from here, as out contains all the info of fluxes and reservoirs during all the years that takes the span of the fitting procedure.
#' @param ... Arguments needed by FuncFitness. Default is none.
#' @return This function does return a vector of the residuals produced by FUN for the LM algorithm.
#' @keywords Fitness IC
FITNESS <- function(FUN = fitness_INC_MORT, obj, ...) { 
    residuals <- FUN(obj, ...) 
    return(residuals) 
}

#' This function is used for call the fitting part of the code within R. It serves as an interface between R and the c codes. It is only invoked once Lm has finished and it produces the initial conditions that will be used when running the simulation.
#'
#' @param ctr A string containing the name of the country in which the simulation is launched.
#' @param itr An integer indicating the number of the error calculation iteration.
#' @param fit_span An integer indicating the number of years that we use for the fitting procedure. Defaults to 19 that are the available years of WHO data in the database, which is used by the default fitness funtion.
#' @param mod_cts An integer indicating how contact matrices should evolve (0-> method 0, 3 -> method 3). Default is 3, which means that contact matrix evolves with demography and preserving connectivity
#' @param MFILE A string with the path of the masterfile of database paths. It is build automatically and should not be modified by user, as code will build and pass it when neccesary.
#' @param FuncFitness A function for calculate fitness. It can be user-defined and should return a vector of residuals after evaluation of fitness criterion. Default is set to WHO data fitness.  
#' @param ... Arguments needed by FuncFitness. Default is none.
#' @return This function does return a vector of the residuals produced by func_fitness for the LM algorithm.
#' @keywords Fitness IC
fitter <- function(ctr, itr, fit_span=19, mod_cts=0, MFILE, func_fitness = fitness_INC_MORT, ...) {
    temp = numeric(3*1531*fit_span)
    temp_C <- .C("model_fit", mf = MFILE, writ = as.integer(1), con_mod = as.integer(mod_cts), objeto = temp)
    lc = 15*19 + 1 + 83*15
    mm <- matrix(temp_C$objeto, nrow = lc)
    mm <- as.data.frame(t(mm))

    auxPath <- path.expand("~")
    if(itr==0){
        pathout <- paste0(auxPath, "/TB_SIM_OUT_Bien_KCaR_TrNorm/", ctr, '/Object_Folder/Baseline')
    } else{
        pathout <- paste0(auxPath, "/TB_SIM_OUT_Bien_KCaR_TrNorm/", ctr, '/Object_Folder/Iter_', sprintf("%03d", itr))
    }
    if (!dir.exists(pathout)){
        dir.create(pathout, recursive = TRUE)
    }
    savepath <- paste0(pathout, "/Fitter_object.txt")
    write.table(mm, savepath, col.names=F, row.names=F, quote=F)
    residuals <- FITNESS(FUN = func_fitness, obj = mm, ...)
    
    line <- sprintf("I P %d %f", itr, sum(residuals))
    write(line, paste0(pathout, "/poids_", ctr, ".txt"))
    return(residuals)
}

#' This function is used for calling the fitting part of the code within R and the LM algorithm. It serves as an interface between R and the c codes. It is only invoked inside LM.
#'
#' @param x A vector containing the values of the parameters that are being calibrated by the LM algorithm.
#' @param ctr A string containing the name of the country in which the simulation is launched.
#' @param itr An integer indicating the number of the error calculation iteration.
#' @param fit_span An integer indicating the number of years that we use for the fitting procedure. Defaults to 19 that are the available years of WHO data in the database, which is used by the default fitness funtion.
#' @param mod_cts An integer indicating how contact matrices should evolve (0-> method 0, 3 -> method 3). Default is 3, which means that contact matrix evolves with demography and preserving connectivity
#' @param MFILE A string with the path of the masterfile of database paths. It is build automatically and should not be modified by user, as code will build and pass it when neccesary.
#' @param FuncFitness A function for calculate fitness. It can be user-defined and should return a vector of residuals after evaluation of fitness criterion. Default is set to WHO data fitness.  
#' @param ... Arguments needed by FuncFitness. Default is none.
#' @return This function does return a vector of the residuals produced by func_fitness for the LM algorithm.
#' @keywords Fitness
fitter_nlfb <- function(x, ctr, itr, fit_span=19, mod_cts=3, MFILE, func_fitness = fitness_INC_MORT, ...) {
    rewritte_parameters_to_fit(x, ctr, itr)
    temp = numeric(3*1531*fit_span)
    temp_C <- .C("model_fit", mf = MFILE, writ = as.integer(0), con_mod = as.integer(mod_cts), objeto = temp)

    lc = 15*19 + 1 + 83*15
    mm <- matrix(temp_C$objeto, nrow = lc)
    mm <- as.data.frame(t(mm))
    
    residuals <- FITNESS(FUN = func_fitness, obj = mm, ...)

    if(sum(residuals) > 2000){
        stop("Bad initial values for parameters, trying again with new values\n")
    }

    print(sum(residuals))
    return(residuals)
}

#' This function is used for calculate the residuals to the all age incidence and mortality data reported by the WHO.
#'
#' @param out An Output_type object produced by both runner and fitter. readOutputPop and readOutputFlux will take this object and calculate incidence and mortality normalized by population.
#' @param ... extra parameters. Here it is passed by other high-level functions the path to the TB burden file where WHO data is.
#' @return This function does return a vector of the residuals for the LM algorithm.
#' @keywords Fitness
fitness_INC_MORT <- function(out, ...) {
    arg <- list(...)
    TB_file <- arg[["TB_burden_file"]][1]
    Burden <- read.table(file=TB_file, header=T)

    reserv <- 19
    scaling <- 1000000
    ind_inc_D <- reserv + c(3, 4, 5, 6, 7, 8, 23, 24, 25, 26, 27, 28, 61, 62, 63, 56, 70, 71)
    ind_mort <- reserv + c(36, 37, 38, 57, 58, 59, 60)
    AG_ini <- 1
    AG_fin <- 15

    inc_new <- readOutputFlux(out, AG_ini, AG_fin, ind = ind_inc_D)
    mort_new <- readOutputFlux(out, AG_ini, AG_fin, ind = ind_mort)
    pop <- readOutputPop(out, AG_ini, AG_fin)   
    
    inc_new <- inc_new/pop*scaling
    mort_new <- mort_new/pop*scaling
    
    ages_fit <- length(Burden$Inc)
    fitness <- numeric(2*ages_fit)
    for(i in 1:ages_fit){
        if((inc_new[i] - Burden$Inc[i]) < 0.0){
            fitness[i] = (inc_new[i] - Burden$Inc[i])/(-mean(Burden$Inc))
        } else{
            fitness[i] = (inc_new[i] - Burden$Inc[i])/(mean(Burden$Inc))
        }
    }
    for(i in (ages_fit+1):(2*ages_fit) ){
        if((mort_new[i-ages_fit] - Burden$Mort[i-ages_fit]) < 0.0){
            fitness[i] = (mort_new[i-ages_fit] - Burden$Mort[i-ages_fit])/(-mean(Burden$Mort))
        } else{
            fitness[i] = (mort_new[i-ages_fit] - Burden$Mort[i-ages_fit])/(mean(Burden$Mort))
        }
    }
    return(fitness)
}

#' A Function for recovering the population asociated to a certain age group or to a sequential combination of age-groups
#'
#' @param obj An Output_type object produced by both runner and fitter.
#' @param a_ini A number indicating the initial age group to start reading, or the group targeted if a_fin defaults to NULL. It has no default.
#' @param a_fin A number indicating the initial age group to stop reading. It defaults to NULL, which make the routine read only the Population asociated to a_ini.
#' @return This function does return the temporal vector of the the population produced by the model. 
#' @keywords Output, Population
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

#' A Function for recovering the fluxes asociated to a certain age group or combination of sequential age-groups
#'
#' @param obj An Output_type object produced by both runner and fitter.
#' @param a_ini A number indicating the initial age group to start reading, or the group targeted if a_fin defaults to NULL. It has no default.
#' @param a_fin A number indicating the initial age group to stop reading. It defaults to NULL, which make the routine read only the fluxes asociated to a_ini.
#' @param ind A vector containing the position of the fluxes in the Output_object, in format 19 (reservories) + indices of flux from 1 to 83. Options are ind_inc_D, ind_inc_T, ind_mort, ind_diag.
#' @return This function does return the temporal vector of the the fluxes asociated to a certain observable that are produced by the model. 
#' @keywords Output, fluxes
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

get_inc_mort_fromObject <- function(obj1, obj2){
    scaling <- 1000000
    ind_inc_D <- 19 + c(3, 4, 5, 6, 7, 8, 23, 24, 25, 26, 27, 28, 61, 62, 63, 56, 70, 71)
    ind_mort <- 19 + c(36, 37, 38, 57, 58, 59, 60)
    AG_ini <- 1
    AG_fin <- 15
    auxPath <- path.expand("~")
    calc <- function(v){
        temp <- numeric(length(v))
        for(i in 1:length(v)){
            temp[i] <- sum(v[1:i])
        }
        return(temp)
    }

    pop <- readOutputPop(obj1, AG_ini, AG_fin)

    i <- readOutputFlux(obj1, AG_ini, AG_fin, ind = ind_inc_D)
    imp <- calc(i)
    m <- readOutputFlux(obj1, AG_ini, AG_fin, ind = ind_mort)
    mort <- calc(m)
    
    i_vac <- readOutputFlux(obj2, AG_ini, AG_fin, ind = ind_inc_D)
    imp_vac <- calc(i_vac)
    m_vac <- readOutputFlux(obj2, AG_ini, AG_fin, ind = ind_mort)
    mort_vac <- calc(m_vac)
    
    df1 <- data.frame(t = 1999 + seq(1:(length(imp))), popu = pop, inc_annual = i, inc_acum = imp, mort_annual = m, mort_acum = mort,
                      inc_cov_annual = i_vac, inc_cov_acum = imp_vac, mort_cov_annual = m_vac, mort_cov_acum = mort_vac)    
    return(df1)
}

#' A Function used for running both fitter and runner simulations.
#'
#' @param ctry A string containing the name of the country in which the simulation is launched.
#' @param iteration An integer indicating the number of the error calculation iteration.
#' @param mode_contact An integer indicating how contact matrices should evolve (0-> method 0, 3 -> method 3). Default is 3, which means that contact matrix evolves with demography and preserving connectivity
#' @param fitting A logical parameter indicating whether fitting step is desired Default is T
#' @param running A logical parameter indicating whether running step is desired. Default is T
#' @param FuncFitness A function for calculate fitness. It can be user-defined and should return a vector of residuals after evaluation of fitness criterion. Default is set to WHO data fitness.  
#' @param ... Arguments needed by FuncFitness. Default is none.
#' @return This function does return two output objects, one for the baseline of the simulation and the other for the vaccine run (Default vaccine of no vaccination profile has been set by user)
simTB_iter <- function(ctry, iteration=0, mode_contact=3, fitting=T, running=T, running_vaccine=F, FuncFitness = fitness_INC_MORT, ...){
    if(iteration==0){
        demography_file = paste0("../Inputs/Input_countries/", ctry, "/Demographic_pyramid.txt")
        TB_file = paste0("../Inputs/Input_countries/", ctry, "/TB_burden.txt")
        fitting_param_file = paste0("../Outputs/", ctry, "/Fitted_inputs/Parameters_to_fit.txt")
    } else{
        demography_file = paste0("../Inputs_Iter/", ctry, "/Demo/Demographic_pyramid_", sprintf("%03d", iteration), ".txt")
        TB_file = paste0("../Inputs_Iter/", ctry, "/Burden/TB_burden_", sprintf("%03d", iteration), ".txt")
        fitting_param_file = paste0("../Outputs/", ctry, "/Iter/", ctry, sprintf("_%03d", iteration), "/Fitted_inputs/Parameters_to_fit.txt")
    }
    
    master_file <- set_TBsim_workspace(ctry, iteration)

    #----------------------------------------Block for time span in both fiter and runner----------------------------------------------#
    Y <- read.table("../Inputs/General_Inputs/window_2100.txt", header=F)
    year_ini_demo <- as.integer(Y$V1) #Simulation's start
    year_fin_demo <- as.integer(Y$V3) #Simulation's end
    years <- year_fin_demo - year_ini_demo + 1 #Simulation's span
    year_min <- as.integer(Y$V1) #Fit's start
    year_max <- as.integer(Y$V2) #Fit's end
    span <- year_max - year_min + 1 #Fit's span
    #----------------------------------------------------------------------------------------------------------------------------------#

    demography_fit(demography_file, ctry, iteration, year_min, year_max, year_ini_demo, year_fin_demo)

    if(fitting){
        data_param <- read.table(fitting_param_file, skip=1, header=T)
        LB <- unlist(data_param$LB)
        UB <- unlist(data_param$UB)
        lb <- numeric(length(data_param$LB))
        ub <- rep(1.0, length(data_param$UB))

        guess <- unlist(data_param$Value)
        p = ifelse( guess < (1 + 1e-10) , 1, 0)
        if(sum(p) < length(p) ){
            print(guess)
            guess <- ( guess - LB) / (UB - LB)
            print(guess)
        }
        names(guess) <- unlist(data_param$Parameter)

        if(require(minpack.lm)){
            library(minpack.lm) 
            auxPath <- path.expand("~")
            if(iteration==0){
                weightPathOut <- paste0(auxPath, "/TB_SIM_OUT_Bien_KCaR_TrNorm/", ctry, '/Object_Folder/Baseline/poids_', ctry, ".txt")
            } else{
                weightPathOut <- paste0(auxPath, "/TB_SIM_OUT_Bien_KCaR_TrNorm/", ctry, '/Object_Folder/Iter_', sprintf("%03d", iteration), '/poids_', ctry, ".txt")
            }

            out <- tryCatch(
            {
                
                modfit <- nls.lm(par = guess, fn = fitter_nlfb, lower = lb, upper = ub, control=c(watch=F, smallsstest=F, rofftest=F, femax=10000),
                                ctr = ctry, itr = iteration, fit_span=span, mod_cts=mode_contact, MFILE=master_file, func_fitness = FuncFitness, ...) #Residual function based
                x <- coef(modfit)
                rewritte_parameters_to_fit(x, ctry, iteration)
                trash <- fitter(ctr = ctry, itr = iteration, fit_span=span, mod_cts=mode_contact, MFILE=master_file, func_fitness = FuncFitness, ...)
                
                wout <- read.table(weightPathOut, header=F)
                x <- as.vector(unlist(wout))
                fitness <- x[4]
                if(fitness > 5.0){
                    line <- sprintf("First fitting stage failure in iteration %03d", iteration)
                    write(line, "Failure_output.txt", append=T)
                    
                    library(nlsr)
                    data_param <- read.table(fitting_param_file, skip=1, header=TRUE)
                    guess <- unlist(data_param$Value)
                    names(guess) <- unlist(data_param$Parameter)
                    
                    modfit <- nlfb(start = guess, resfn = fitter_nlfb, lower = lb, upper = ub, trace = FALSE,
                                ctr = ctry, itr = iteration, fit_span=span, mod_cts=mode_contact, MFILE=master_file, func_fitness = FuncFitness, ...)
                                
                    x <- coef(modfit)
                    rewritte_parameters_to_fit(x, ctry, iteration)
                        
                    trash <- fitter(ctr = ctry, itr = iteration, fit_span=span, mod_cts=mode_contact, MFILE=master_file, func_fitness = FuncFitness, ...)
                    wout2 <- read.table(weightPathOut, header=FALSE)
                    x2 <- as.vector(unlist(wout2))
                    fitness2 <- x2[4]
                    if(fitness2 > 10.0){
                        line <- sprintf("Two stages fitting failure in iteration %03d", iteration)
                        write(line, "Failure_output.txt", append=T)
                        #stop("Two stages fitting failure, recalculating...")
                    }        
                }
                print("All-OK")
            },
            error=function(cond3) {
                stop(cond3)
            }) 

            
        } else{
            stop("No fitting package available. Try insalling minpack.lm and nlsr")
        }
        
        print(sprintf("Iteration %d completed", iteration))
    }


    if(running){
        mcpy = 15*19 + 1 + 83*15
        temp1 = numeric(3*mcpy*years)
        temp2 = numeric(3*mcpy*years)
        error_handler = 0
        if(running_vaccine){
            if(iteration==0){
                adjust <- 0
            } else{
                adjust <- 1
            }
            temp_C <- .C("R2C_runner", mf = master_file, as.integer(adjust), as.integer(0), as.integer(iteration), as.integer(mode_contact), vaccinate=as.integer(1), objeto1 = temp1, objeto2 = temp2, out_handler = as.integer(error_handler))
            if(as.integer(temp_C$out_handler)==1) stop("Running error")
            #save output
            mm1 <- matrix(temp_C$objeto1, nrow = mcpy )
            mm1 <- as.data.frame(t(mm1))
            mm2 <- matrix(temp_C$objeto2, nrow = mcpy )
            mm2 <- as.data.frame(t(mm2))
        } else{
            temp_C <- .C("R2C_runner", mf = master_file, as.integer(0), as.integer(0), as.integer(iteration), as.integer(mode_contact), vaccinate=as.integer(0), objeto1 = temp1, objeto2 = temp2, out_handler = as.integer(error_handler))
            if(as.integer(temp_C$out_handler)==1) stop("Running error")
            #save output
            mm1 <- matrix(temp_C$objeto1, nrow = mcpy )
            mm1 <- as.data.frame(t(mm1))
            
        }
        
        auxPath <- path.expand("~")
        if(iteration==0){
            pathout <- paste0(auxPath, "/TB_SIM_OUT_Bien_KCaR_TrNorm/", ctry, '/Object_Folder/Baseline')
        } else{
            pathout <- paste0(auxPath, "/TB_SIM_OUT_Bien_KCaR_TrNorm/", ctry, '/Object_Folder/Iter_', sprintf("%03d", iteration))
        }
        if (!dir.exists(pathout)){
            dir.create(pathout, recursive = TRUE)
        }
        savepath1 <- paste0(pathout, "/Output_object_1.txt")
        write.table(mm1, savepath1, col.names=F, row.names=F, quote=F)
        if(running_vaccine){
            savepath2 <- paste0(pathout, "/Output_object_2.txt")
            write.table(mm2, savepath2, col.names=F, row.names=F, quote=F)
        }
        
        save_vaccine = TRUE
        if(save_vaccine && running_vaccine){
            vac_path <- paste0("../Inputs/General_Inputs/perfil.txt")
            info <- get_vaccine_info(vac_path)
            aux_tran <- unlist(info[1])
            target_type <- unlist(info[2])
            aux_info <- unlist(info[3])
            
            line_tran <- target_type
            for(h in 1:length(aux_tran)){
                aux_paste <- paste0(line_tran, "_", aux_tran[h])
                line_tran <- aux_paste
            }
            
            if(iteration==0){
                pathout_vaccine <- paste0(auxPath, "/TB_SIM_OUT_Bien_KCaR_TrNorm/", ctry, "/Vaccines/Vac_",
                                sprintf("T_%d-E_%.3f-C_%.3f-W_%.2f", as.integer(aux_info["age"]), aux_info["efficacy"], aux_info["coverage"], aux_info["waning"] ), 
                                "/", line_tran, "/Baseline/")
            } else{
                pathout_vaccine <- paste0(auxPath, "/TB_SIM_OUT_Bien_KCaR_TrNorm/", ctry, "/Vaccines/Vac_",
                                sprintf("T_%d-E_%.3f-C_%.3f-W_%.2f", as.integer(aux_info["age"]), aux_info["efficacy"], aux_info["coverage"], aux_info["waning"] ), 
                                "/", line_tran, "/Iters/")
            }
            
            if (!dir.exists(pathout_vaccine)){
                dir.create(pathout_vaccine, recursive = TRUE)
            }
            out_df <- get_inc_mort_fromObject(mm1, mm2)
            if(iteration==0){
                savepath1 <- paste0(pathout_vaccine, "Output_vac.txt")
                write.table(out_df, savepath1, col.names=T, row.names=F, quote=F)
            } else{
                savepath1 <- paste0(pathout_vaccine, "Output_vac_It_", sprintf("%03d", iteration), ".txt")
                write.table(out_df, savepath1, col.names=T, row.names=F, quote=F)
            }
            
            #if(iteration==0){
            #    savepath1 <- paste0(pathout_vaccine, "Output_object_1.txt")
            #    write.table(mm1, savepath1, col.names=F, row.names=F, quote=F)
            #    savepath2 <- paste0(pathout_vaccine, "Output_object_2.txt")
            #    write.table(mm2, savepath2, col.names=F, row.names=F, quote=F)
            #} else{
            #    savepath1 <- paste0(pathout_vaccine, "Output_object_1_It_", sprintf("%03d", iteration),".txt")
            #    write.table(mm1, savepath1, col.names=F, row.names=F, quote=F)
            #    savepath2 <- paste0(pathout_vaccine, "Output_object_2_It_", sprintf("%03d", iteration),".txt")
            #    write.table(mm2, savepath2, col.names=F, row.names=F, quote=F)
            #}
            
        }
        
    }

}

#' A Function for getting the efficacy, coverage, age group targeted of a given vaccine.
#' @param path_to_profile Path to vaccine profile txt file.
#' @return This function returns the info needed to save the outputs of a given vaccine.
#' @keywords vaccine
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

#' A Function used for getting the maximum number of allowed iterations in order to not overcome R's memory limit.
#' @param objSize A float meaning the size in MB of the output object. Size is obtained with object.size(). Defaults to 3MB
#' @param MemPercentUsage A float meaning the percentage of maximum memory that we want to occupy. Defaults to 20\%
#' @return This function returns an integer.
#' @keywords Nsim memory
#' @examples
#' get_N_sim(3, 0.5)
get_N_sim <- function(objSize=3, MemPercentUsage=0.2){
    if(objSize<=0){
        warning("value of objSize <=0\n Raising error")
        stop("Size of Output Object<=0. Unable to compute Nsim. Check object generation")
    }
    if(MemPercentUsage<0){
        warning("value of MemPercentUsage <0\n MemPercentUsage set =0")
        MemPercentUsage = max(0, MemPercentUsage)
    }
    if(MemPercentUsage>1){
        warning("value of MemPercentUsage > 1\n MemPercentUsage set =1")
        MemPercentUsage = min(1, MemPercentUsage)
    }
    if(.Platform$OS.type == "unix") {
        if(.Platform$pkgType=="source"){
            x <- system('grep MemTotal /proc/meminfo', intern = TRUE)
            x <- unlist(strsplit(x, " "))
            y <- grepl("\\d", x)
            index <- which(T==y)[[1]]
            Mem_R <- as.integer(x[index])
        } else{
            x <- system('sysctl hw.memsize', intern=T)
            x <- unlist(strsplit(x, " "))
            y <- grepl("\\d", x)
            index <- which(T==y)[[1]]
            Mem_R <- as.numeric(x[index])
        }
        
    } else {
        Mem_R <- memory.limit()
    }
    
    N_sim <- Mem_R*MemPercentUsage*0.85/(1024*2*objSize)
    N_sim <- round(N_sim)
    return(N_sim)
}
