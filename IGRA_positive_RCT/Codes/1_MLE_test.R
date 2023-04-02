#!/usr/bin/Rscript

library(bbmle)
library(fanplot)
library(ggplot2)
require(gridExtra)
library(kdensity)
library(KernSmooth)
library(parallel)
library(iterators)
library(foreach)
library(doParallel)

if(!file.exists("trial_frac_cases.so")){
    if(!file.exists("trial_frac_cases.c")){
        print(paste0('C library not available in current directory. Raising error...'))
        stop()
    }
    print("Building shared library...")
    system("R CMD SHLIB trial_frac_cases.c")
}
dyn.load("trial_frac_cases.so")

#-------------------------------------------------------------------------------NUBE-----------------------------------------------------------------

vaccine_pop <- c(706, 339, 581) #4-5 6 7-11
control_pop <- c(724, 321, 594)
pop_vac = integer(15)
pop_cont = integer(15)
probs <- numeric(15)

#South_Africa
a <- read.table("Data_countries/demo_2015_SA.txt", header=TRUE)
x <- a$value

probs[4:5] = x[4:5]/sum(x[4:5])
probs[6] = 1
probs[7:10] = x[7:10]/sum(x[7:10])

pop_vac[4] = round(vaccine_pop[1]*probs[4])
pop_vac[5] = vaccine_pop[1] - pop_vac[4]
pop_vac[6] = vaccine_pop[2]
for(i in 7:9){
    pop_vac[i] = round(vaccine_pop[3]*probs[i])
}
pop_vac[10] = vaccine_pop[3] - (pop_vac[7] + pop_vac[8] + pop_vac[9])

pop_cont[4] = round(control_pop[1]*probs[4])
pop_cont[5] = control_pop[1] - pop_cont[4]
pop_cont[6] = control_pop[2]
for(i in 7:9){
    pop_cont[i] = round(control_pop[3]*probs[i])
}
pop_cont[10] = control_pop[3] - (pop_cont[7] + pop_cont[8] + pop_cont[9])

PR_C_SA = pop_cont/sum(pop_cont)
PR_V_SA = pop_vac/sum(pop_vac)

#Kenya
a <- read.table("Data_countries/demo_2015_Kenya.txt", header=TRUE)
x <- a$value

probs[4:5] = x[4:5]/sum(x[4:5])
probs[6] = 1
probs[7:10] = x[7:10]/sum(x[7:10])

pop_vac[4] = round(vaccine_pop[1]*probs[4])
pop_vac[5] = vaccine_pop[1] - pop_vac[4]
pop_vac[6] = vaccine_pop[2]
for(i in 7:9){
    pop_vac[i] = round(vaccine_pop[3]*probs[i])
}
pop_vac[10] = vaccine_pop[3] - (pop_vac[7] + pop_vac[8] + pop_vac[9])

pop_cont[4] = round(control_pop[1]*probs[4])
pop_cont[5] = control_pop[1] - pop_cont[4]
pop_cont[6] = control_pop[2]
for(i in 7:9){
    pop_cont[i] = round(control_pop[3]*probs[i])
}
pop_cont[10] = control_pop[3] - (pop_cont[7] + pop_cont[8] + pop_cont[9])

PR_C_K = pop_cont/sum(pop_cont)
PR_V_K = pop_vac/sum(pop_vac)

#Zambia
a <- read.table("Data_countries/demo_2015_Zambia.txt", header=TRUE)
x <- a$value

probs[4:5] = x[4:5]/sum(x[4:5])
probs[6] = 1
probs[7:10] = x[7:10]/sum(x[7:10])

pop_vac[4] = round(vaccine_pop[1]*probs[4])
pop_vac[5] = vaccine_pop[1] - pop_vac[4]
pop_vac[6] = vaccine_pop[2]
for(i in 7:9){
    pop_vac[i] = round(vaccine_pop[3]*probs[i])
}
pop_vac[10] = vaccine_pop[3] - (pop_vac[7] + pop_vac[8] + pop_vac[9])

pop_cont[4] = round(control_pop[1]*probs[4])
pop_cont[5] = control_pop[1] - pop_cont[4]
pop_cont[6] = control_pop[2]
for(i in 7:9){
    pop_cont[i] = round(control_pop[3]*probs[i])
}
pop_cont[10] = control_pop[3] - (pop_cont[7] + pop_cont[8] + pop_cont[9])

PR_C_Z = pop_cont/sum(pop_cont)
PR_V_Z = pop_vac/sum(pop_vac)

LL_sto_2 <- function(epsilon) {
    out_handler <- numeric(Ntrials)
    aux <- .C("start_sto_sim", AG_C_SA=as.numeric(PR_C_SA), AG_V_SA=as.numeric(PR_V_SA), AG_C_K=as.numeric(PR_C_K), AG_V_K=as.numeric(PR_V_K), AG_C_Z=as.numeric(PR_C_Z), AG_V_Z=as.numeric(PR_V_Z), mode=as.integer(mode_gen), epsilon=as.numeric(epsilon), VEdis=as.numeric(out_handler))
    VE_dis_mod = aux$VEdis
    return(VE_dis_mod)
}

fitG = function(x,y,mu,sig,scale){
    f = function(p){
        d = p[3]*dnorm(x,mean=p[1],sd=p[2])
        sum((d-y)^2)
    }
    optim(c(mu,sig,scale),f)
}

user_made_bootstrap <- function(x, Replicates=5000, sample_max=1000000){
    
    numCores <- detectCores()
    registerDoParallel(numCores)
    
    total_ind <- 1:length(x)
    
    like <- foreach (h=1:Replicates, .combine = 'c') %dopar% {
        ind <- sample(total_ind, sample_max, replace=TRUE)
        new_x <- x[ind]
        p_vd <- bkde(new_x)
        i <- which.min(abs(p_vd$x - 49.7))
        p_vd$y[i]
    }
    return(like)
}

Ntrials <- 10000 #10000 Hay que cambiarlo también en el .c
eps_axix <- seq(0, 0.995, 0.005)
bootmax <- 5000
modes = 1:7

df_aux <- data.frame(it=1:bootmax, pm1=rep(0, bootmax), pm2=rep(0, bootmax), pm3=rep(0, bootmax), pm4=rep(0, bootmax), pm5=rep(0, bootmax), pm6=rep(0, bootmax), pm7=rep(0, bootmax))
likelihood = numeric(length(modes))
likelihood_low = numeric(length(modes))
likelihood_up = numeric(length(modes))
likelihood_med = numeric(length(modes))
    
age_gen <- "Mixed_pop"
scale_par <- c(.01, .001, 1, .1, .01, .1, .1)

line <- sprintf("Model eps_estim eps_low eps_hi")
write(line, paste0("Outputs/Eps_estim_A_", age_gen, ".txt"))

for(mode_gen in modes){
    output_dir <- paste0("Outputs/Model_", sprintf("%d", mode_gen), "/")
    if (!dir.exists(output_dir)){
        dir.create(output_dir)
    } 

    nube <- matrix(0, Ntrials, length(eps_axix))
    nube_bin <- matrix(0, length(eps_axix), Ntrials)
    nube_data <- matrix(0, length(eps_axix), 4)
    peps <- numeric(length(eps_axix))

    print("Computing\n")
    pb <- txtProgressBar(min = 1, max = length(eps_axix), style = 3)
    sent=1
    for(eps in eps_axix){
        x <- LL_sto_2(eps)
          
        aux <- tryCatch({
            #---------------------------------------# Cálculo de P(eps|VEdis)
            p_vd <- bkde(x)
            ind <- which.min(abs(p_vd$x - 49.7))
            peps[sent] = p_vd$y[ind]
                
        },
        error=function(condition2){
            print('ERROR!')
            print(condition2)
            print(sent)
                
            FUN = kdensity(x, kernel = 'gaussian', na.rm = FALSE, normalized = TRUE)
            peps[sent] = FUN(49.7)
                
                
        },
        warning=function(condition2) {
            print('WARNING!')
            print(condition2)
            print(sent)
        }) 
            
        nube[1:(Ntrials), sent] = x
        nube_bin[sent, 1:(Ntrials)] = x
        nube_data[sent, 1:4] = c(eps, median(x), min(x), max(x))
            
        sent = sent + 1
        setTxtProgressBar(pb, sent)
    }
    close(pb)
    line <- sprintf("eps_estim VEdis")
    write(line, paste0("Outputs/Model_", sprintf("%d", mode_gen), "/test_sto_A_", age_gen, ".txt"))
    write.table(nube, paste0("Outputs/Model_", sprintf("%d", mode_gen), "/test_sto_A_", age_gen, ".txt"), append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
        
    #Normalize peps
    #peps[peps>1]=NaN
    df_out_aux <- data.frame(eps=eps_axix, prob=peps)
    write.table(df_out_aux, paste0("Outputs/Model_", sprintf("%d", mode_gen), "/P(VE_dis|eps)_", age_gen, ".txt"), row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    normalize <- sum(peps)
    peps <- peps/normalize
    plot(eps_axix, peps)
    df_out <- data.frame(eps=eps_axix, prob=peps)
    write.table(df_out, paste0("Outputs/Model_", sprintf("%d", mode_gen), "/P(eps|VE_dis)_", age_gen, ".txt"), row.names=FALSE, col.names=TRUE, quote=FALSE)
        
    #Opción para modelos con valores extremos
    param <- fitG(eps_axix, peps, 0.5, 0.2, scale_par[mode_gen])$par
    print(param)
    if(param[1] < 0 || param[2] < 0){
        param[1] = 0
        param[2] = 1
    }
        
    gauss_pdf <- function(x){
        return(dnorm(x, param[1], param[2]))
    }
        
    #Plot the fit of p(eps) to gaussian
    aux_plot <- ggplot() +
            geom_point(aes(x=eps_axix, y=peps)) +
            geom_line(aes(x=eps_axix, y=param[3]*dnorm(eps_axix, param[1], param[2])), color='red') +
            labs(x = expression(epsilon), y = expression('Pr('~epsilon~'|VEdis)') ) +
            theme(axis.line = element_line(colour = "black", size=0.3), legend.title=element_blank() )
    ggsave(paste0('Outputs/Model_', sprintf("%d", mode_gen), '/eps_dist_fit_A_', age_gen, '.png'), aux_plot)

    p_low_norm <- 0.025*(1.0 - integrate(gauss_pdf, 1, Inf)$value)
    if(qnorm(.975, param[1], param[2])>1){
        eps_hi_G = 1
    } else{
        p_hi_norm <- 0.975*(1.0 - integrate(gauss_pdf, 1, Inf)$value)
        eps_hi_G = qnorm(p_hi_norm, param[1], param[2])
    }

    eps_low_G = qnorm(p_low_norm, param[1], param[2])
    eps_med_G = qnorm(0.5, param[1], param[2])
    
    eps_med = qnorm(0.5, param[1], param[2])
    eps_low = qnorm(0.025, param[1], param[2])
    eps_hi = qnorm(.975, param[1], param[2])

    line <- sprintf("%d %f %f %f %f %f %f", mode_gen, eps_med_G, eps_low_G, eps_hi_G, eps_med, eps_low, eps_hi)
    write(line, paste0("Outputs/Eps_estim_A_", age_gen, ".txt"), append=T)
       
    df <- data.frame(eps = eps_axix, yval = nube_data[, 2], y_min = nube_data[, 3], y_max = nube_data[, 4])
        
    vector_data <- as.vector(nube_bin)
    p_vd <- bkde(vector_data)
    ind <- which.min(abs(p_vd$x - 49.7))
    likelihood[mode_gen] = p_vd$y[ind]
    like_CI <- user_made_bootstrap(vector_data, bootmax)
    
    df_aux[,mode_gen+1] <- like_CI
    
    out_ridge_plot_likelihood <- data.frame(ID=1:length(like_CI), Data=like_CI)
    write.table(out_ridge_plot_likelihood, paste0('Outputs/Model_', sprintf("%d", mode_gen), '/like_CI_', sprintf("%d", mode_gen), '_A_', age_gen, '.txt'), row.names = F, quote = F)
        
}

#Normalize likelihoods
normalize <- sum(likelihood)
likelihood <- likelihood/normalize

df_aux$norm=rowSums(df_aux[,2:8])
df_aux[, 2:8] <- df_aux[,2:8]/df_aux$norm
for(mm in modes){
    like_CI <- df_aux[, mm+1]
    
    likelihood_low[mm] <- quantile(like_CI, .025)
    likelihood_up[mm] <- quantile(like_CI, .975)
    likelihood_med[mm] <- quantile(like_CI, .5)
}

dfL <- data.frame(Model = modes, LM = likelihood, LM_low = likelihood_low, LM_up = likelihood_up, LM_med = likelihood_med)
write.table(dfL, paste0("Outputs/Likelihood_values.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE)

