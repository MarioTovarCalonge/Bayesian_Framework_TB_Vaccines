#!/usr/bin/Rscript

#library
library(ggplot2)

if(!file.exists("trial_frac_cases.so")){
    if(!file.exists("trial_frac_cases.c")){
        print(paste0('C library not available in current directory. Raising error...'))
        stop()
    }
    print("Building shared library...")
    system("R CMD SHLIB trial_frac_cases.c")
}
dyn.load("trial_frac_cases.so")

#-------------------age dist----------------------#

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


launch_control_branch <- function() {
    aux <- .C("start_sto_sim", AG_C_SA=as.numeric(PR_C_SA), AG_V_SA=as.numeric(PR_V_SA), AG_C_K=as.numeric(PR_C_K), AG_V_K=as.numeric(PR_V_K), AG_C_Z=as.numeric(PR_C_Z), AG_V_Z=as.numeric(PR_V_Z))
}

#--------------------Sims----------------------#

output_dir <- paste0("Case_Dist/")
if (!dir.exists(output_dir)){
    dir.create(output_dir)
}

print("Computing\n")

x <- launch_control_branch()
          
#Load cases data and produce plot:
#countries are loaded as 0 -> South Africa, 1 -> Kenya, 2 -> Zambia

#Plot each country
axisx <- c("0-5", "6-10", "11-15", "16-20", "21-25", "26-30", "31-35", "36-40", "41-45", "46-50", "51-55", "56-60", "61-65", "66-70", "70+")
countries <- c("South_Africa", "Kenya", "Zambia")

plist <- list()

for(ctrn in 0:2){
    print(countries[ctrn+1])
    df <- read.table(paste0("Case_Dist/Fracs_", ctrn, ".txt"), header=TRUE)
    fp <- numeric(15)
    rfp <- numeric(15)
    sp <- numeric(15)
        
    for(a in 0:14){
        dfaux <- subset(df, Age==a)
        
        aux <- colSums(dfaux)
        fp_aux <- aux["FastProg"]
        rfp_aux <- aux["ReinfFP"]
        sp_aux <- aux["SlowProg"]
        d_aux <- fp_aux + rfp_aux + sp_aux
            
        fp[a+1] <- fp_aux/d_aux*100
        rfp[a+1] <- rfp_aux/d_aux*100
        sp[a+1] <- sp_aux/d_aux*100
    }

    df_out <- data.frame(type=c(rep("ReinfFP", 15), rep("FastProg", 15), rep("SlowProg", 15)),
                        age=c(rep(axisx, 3)),
                        value=c(rfp, fp, sp))
    df_out$age <- as.factor(df_out$age)

    write.table(df_out, paste0("Outputs/Events_per_age_", countries[ctrn+1], ".txt"), row.names=FALSE, col.names=TRUE, quote=FALSE)

    p_dist <- ggplot(df_out, aes(fill=type, y=value, x=age)) +
       geom_bar(position="stack", stat="identity") + 
       scale_x_discrete(limits = axisx) +
       scale_fill_manual(values = c("red3", "deepskyblue2", "tan1")) +
       theme_minimal() +
       theme(axis.line = element_blank(),
           text = element_text(size=15, family="serif"),
           axis.text.x = element_text(angle = 30),
           axis.title = element_blank(),
           legend.title = element_blank(),
           panel.grid.minor = element_blank())
           
    plist[[ctrn+1]] <- p_dist

    pdf(paste0("../Case_dist_", countries[ctrn+1], ".pdf"))
    print(p_dist)
    dev.off()
}

