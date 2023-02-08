#!/usr/bin/Rscript

#library
library(ggplot2)

#-------------------age dist----------------------#

control_pop <- c(724, 321, 594)
pop_cont = integer(15)
probs <- numeric(15)

#South_Africa
a <- read.table("Data_countries/demo_2015_SA.txt", header=TRUE)
x <- a$value

probs[4:5] = x[4:5]/sum(x[4:5])
probs[6] = 1
probs[7:10] = x[7:10]/sum(x[7:10])

pop_cont[4] = round(control_pop[1]*probs[4])
pop_cont[5] = control_pop[1] - pop_cont[4]
pop_cont[6] = control_pop[2]
for(i in 7:9){
    pop_cont[i] = round(control_pop[3]*probs[i])
}
pop_cont[10] = control_pop[3] - (pop_cont[7] + pop_cont[8] + pop_cont[9])

PR_C_SA = pop_cont/sum(pop_cont)

#Kenya
a <- read.table("Data_countries/demo_2015_Kenya.txt", header=TRUE)
x <- a$value

probs[4:5] = x[4:5]/sum(x[4:5])
probs[6] = 1
probs[7:10] = x[7:10]/sum(x[7:10])

pop_cont[4] = round(control_pop[1]*probs[4])
pop_cont[5] = control_pop[1] - pop_cont[4]
pop_cont[6] = control_pop[2]
for(i in 7:9){
    pop_cont[i] = round(control_pop[3]*probs[i])
}
pop_cont[10] = control_pop[3] - (pop_cont[7] + pop_cont[8] + pop_cont[9])

PR_C_K = pop_cont/sum(pop_cont)

#Zambia
a <- read.table("Data_countries/demo_2015_Zambia.txt", header=TRUE)
x <- a$value

probs[4:5] = x[4:5]/sum(x[4:5])
probs[6] = 1
probs[7:10] = x[7:10]/sum(x[7:10])

pop_cont[4] = round(control_pop[1]*probs[4])
pop_cont[5] = control_pop[1] - pop_cont[4]
pop_cont[6] = control_pop[2]
for(i in 7:9){
    pop_cont[i] = round(control_pop[3]*probs[i])
}
pop_cont[10] = control_pop[3] - (pop_cont[7] + pop_cont[8] + pop_cont[9])

PR_C_Z = pop_cont/sum(pop_cont)
agew <- matrix(c(PR_C_SA, PR_C_K, PR_C_Z), ncol=3)

ctrw <- c(0.808178, 0.1491281, 0.04269393)

#--------------------Sims----------------------#
          
#Load cases data and produce plot:
#countries are loaded as 0 -> South Africa, 1 -> Kenya, 2 -> Zambia

#Plot each country
axisx <- c("0-5", "6-10", "11-15", "16-20", "21-25", "26-30", "31-35", "36-40", "41-45", "46-50", "51-55", "56-60", "61-65", "66-70", "70+")
countries <- c("South_Africa", "Kenya", "Zambia")

plist <- list()

fp <- numeric(3)
rfp <- numeric(3)
sp <- numeric(3)

for(ctrn in 0:2){
    print(countries[ctrn+1])
    df <- read.table(paste0("Case_Dist/Fracs_", ctrn, ".txt"), header=TRUE)
        
    for(a in 3:9){
        dfaux <- subset(df, Age==a)
        
        aux <- colSums(dfaux)
        fp_aux <- aux["FastProg"]
        rfp_aux <- aux["ReinfFP"]
        sp_aux <- aux["SlowProg"]

        fp[ctrn+1] <- fp[ctrn+1] + fp_aux*agew[a+1, ctrn+1]
        rfp[ctrn+1] <- rfp[ctrn+1] + rfp_aux*agew[a+1, ctrn+1]
        sp[ctrn+1] <- sp[ctrn+1] + sp_aux*agew[a+1, ctrn+1]
    }
}

cpa_fp <- fp[1]*ctrw[1] + fp[2]*ctrw[2] + fp[3]*ctrw[3]
cpa_rfp <- rfp[1]*ctrw[1] + rfp[2]*ctrw[2] + rfp[3]*ctrw[3]
cpa_sp <- sp[1]*ctrw[1] + sp[2]*ctrw[2] + sp[3]*ctrw[3]
aD <- cpa_fp + cpa_rfp + cpa_sp

per_cpa_fp <- cpa_fp/aD*100
per_cpa_rfp <- cpa_rfp/aD*100
per_cpa_sp <- cpa_sp/aD*100

line <- sprintf("FP: %f \n SP: %f \n rFP: %f \n", per_cpa_fp, per_cpa_sp, per_cpa_rfp)
print(line)

df_out <- data.frame(tg=c("FastProg", "ReinfFP", "SlowProg"), val=c(per_cpa_fp, per_cpa_rfp, per_cpa_sp))
write.table(df_out, "Outputs/cpa_events.txt", quote=FALSE, row.names=FALSE)

p_cpa <- ggplot(df_out, aes(fill=tg, y=val, x="")) +
   geom_bar(position="stack", stat="identity") +
    labs(x="", y="Fraction of total events (%)") +
   scale_fill_manual(values = c("red3", "deepskyblue2", "tan1")) +
   theme_minimal() +
   theme(axis.line = element_blank(),
       text = element_text(size=15, family="serif"),
       legend.title = element_blank(),
       panel.grid.minor = element_blank())

pdf("../cpa_events.pdf", width=3)
print(p_cpa)
dev.off()
