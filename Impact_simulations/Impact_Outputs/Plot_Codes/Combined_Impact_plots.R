#!/usr/bin/Rscript

library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
pdf(NULL)


df <- read.table("Impact_all_countries.txt", header=TRUE)
df1 <- subset(df, Country=="India")
df1_1 <- subset(df1, Model!='All')
df1_2 <- subset(df1, Model=='All')

p1 <- ggplot(data = df1_1, aes(x=Model, y=IRR)) +
        geom_bar(stat="identity", fill='lightskyblue', width=0.8) +
        
        geom_errorbar(aes(ymax = IRR_hi, ymin = IRR_low), position='dodge', width = .3) +
        
        ylim(0, 6) +
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
        theme_minimal() + 
        theme(legend.title=element_blank(), axis.title.x=element_blank()) + 
        labs(x = "", y = "") +
        ggtitle("A: India")

df1_2$Model <- c("Bayesian")
p1_aux <- ggplot(data = df1_2, aes(x=Model, y=IRR)) +
            geom_bar(stat="identity", fill='wheat', width=1.2) +

            geom_errorbar(aes(ymax = IRR_hi, ymin = IRR_low), position='dodge', width = .4) +
    
            ylim(0, 6) +
            ggtitle("") +
            theme_minimal() +
            theme(legend.title=element_blank(), axis.title.x=element_blank()) +
            labs(x = "", y = "") +
            theme(axis.text.y=element_blank(),
                axis.ticks.y=element_blank())

p1_comb <- plot_grid(p1, p1_aux, rel_widths=c(0.875, .15))

pdf("All_Eff_India.pdf", height=5, width=7)
print(p1_comb)
dev.off()

df2 <- subset(df, Country=="Indonesia")
df2_1 <- subset(df2, Model!='All')
df2_2 <- subset(df2, Model=='All')

p2 <- ggplot(data = df2_1, aes(x=Model, y=IRR)) +
        geom_bar(stat="identity", fill='lightskyblue', width=0.8) +
        
        geom_errorbar(aes(ymax = IRR_hi, ymin = IRR_low), position='dodge', width = .3) +
        
        ylim(0, 8) +
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
        theme_minimal() +
        theme(legend.title=element_blank(), axis.title.x=element_blank()) +
        labs(x = "", y = "") +
        ggtitle("B: Indonesia")

df2_2$Model <- c("Bayesian")
p2_aux <- ggplot(data = df2_2, aes(x=Model, y=IRR)) +
            geom_bar(stat="identity", fill='wheat', width=1.2) +

            geom_errorbar(aes(ymax = IRR_hi, ymin = IRR_low), position='dodge', width = .4) +
    
            ylim(0, 8) +
            ggtitle("") +
            theme_minimal() +
            theme(legend.title=element_blank(), axis.title.x=element_blank()) +
            labs(x = "", y = "") +
            theme(axis.text.y=element_blank(),
                axis.ticks.y=element_blank())

p2_comb <- plot_grid(p2, p2_aux, rel_widths=c(0.875, .15))
pdf("All_Eff_Indonesia.pdf", height=5, width=7)
print(p2_comb)
dev.off()

df3 <- subset(df, Country=="Ethiopia")
df3_1 <- subset(df3, Model!='All')
df3_2 <- subset(df3, Model=='All')

p3 <- ggplot(data = df3_1, aes(x=Model, y=IRR)) +
        geom_bar(stat="identity", fill='lightskyblue', width=0.8) +
        
        geom_errorbar(aes(ymax = IRR_hi, ymin = IRR_low), position='dodge', width = .3) +
        
        ylim(0, 25) +
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
        theme_minimal() +
        theme(legend.title=element_blank(), axis.title.x=element_blank()) +
        labs(x = "", y = "") +
        ggtitle("C: Ethiopia")

df3_2$Model <- c("Bayesian")
p3_aux <- ggplot(data = df3_2, aes(x=Model, y=IRR)) +
            geom_bar(stat="identity", fill='wheat', width=1.2) +

            geom_errorbar(aes(ymax = IRR_hi, ymin = IRR_low), position='dodge', width = .4) +
    
            ylim(0, 25) +
            ggtitle("") +
            theme_minimal() +
            theme(legend.title=element_blank(), axis.title.x=element_blank()) +
            labs(x = "", y = "") +
            theme(axis.text.y=element_blank(),
                axis.ticks.y=element_blank())

p3_comb <- plot_grid(p3, p3_aux, rel_widths=c(0.875, .15))
pdf("All_Eff_Ethiopia.pdf", height=5, width=7)
print(p3_comb)
dev.off()


#Ahora todos juntos con cowplot
p_all <- plot_grid(p1_comb, p2_comb, p3_comb, ncol=1)

#create common x and y labels
y.grob <- textGrob("IRR (%)", gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)
x.grob <- textGrob("Model", gp=gpar(fontface="bold", col="black", fontsize=15))
    
#add to plot and save
pdf("../All_Eff_Combined.pdf", height=10, width=7)
grid.arrange(arrangeGrob(p_all, left = y.grob, bottom = x.grob))
dev.off()

