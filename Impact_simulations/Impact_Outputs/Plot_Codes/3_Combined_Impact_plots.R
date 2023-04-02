#!/usr/bin/Rscript

library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
pdf(NULL)

ylimall <- 15
fsizeText <- 12
fsizeTit <- 20
fsizeColnames <- 25

#neonatos
target <- 0

df <- read.table(paste0("Impact_all_countries_T_", target, ".txt"), header=TRUE)
df1 <- subset(df, Country=="India")
df1_1 <- subset(df1, Model!='All')
df1_2 <- subset(df1, Model=='All')

p1 <- ggplot(data = df1_1, aes(x=Model, y=IRR)) +
        geom_bar(stat="identity", fill='lightskyblue', width=0.8) +
        
        geom_errorbar(aes(ymax = IRR_hi, ymin = IRR_low), position='dodge', width = .3) +
        
        ylim(0, ylimall) +
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
        theme_minimal() + 
        theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.text = element_text(size=fsizeText), plot.title=element_text(size=fsizeTit)) +
        labs(x = "", y = "") +
        ggtitle("A: India")

df1_2$Model <- c("Bayesian")
p1_aux <- ggplot(data = df1_2, aes(x=Model, y=IRR)) +
            geom_bar(stat="identity", fill='wheat', width=1.2) +

            geom_errorbar(aes(ymax = IRR_hi, ymin = IRR_low), position='dodge', width = .4) +
    
            ylim(0, ylimall) +
            ggtitle("") +
            theme_minimal() +
            theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.text = element_text(size=fsizeText)) +
            labs(x = "", y = "") +
            theme(axis.text.y=element_blank(),
                axis.ticks.y=element_blank())

p1_comb <- plot_grid(p1, p1_aux, rel_widths=c(0.875, .15))

pdf(paste0("All_Eff_India_T_", target, ".pdf"), height=5, width=7)
print(p1_comb)
dev.off()

df2 <- subset(df, Country=="Indonesia")
df2_1 <- subset(df2, Model!='All')
df2_2 <- subset(df2, Model=='All')

p2 <- ggplot(data = df2_1, aes(x=Model, y=IRR)) +
        geom_bar(stat="identity", fill='lightskyblue', width=0.8) +
        
        geom_errorbar(aes(ymax = IRR_hi, ymin = IRR_low), position='dodge', width = .3) +
        
        ylim(0, ylimall) +
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
        theme_minimal() +
        theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.text = element_text(size=fsizeText), plot.title=element_text(size=fsizeTit)) +
        labs(x = "", y = "") +
        ggtitle("B: Indonesia")

df2_2$Model <- c("Bayesian")
p2_aux <- ggplot(data = df2_2, aes(x=Model, y=IRR)) +
            geom_bar(stat="identity", fill='wheat', width=1.2) +

            geom_errorbar(aes(ymax = IRR_hi, ymin = IRR_low), position='dodge', width = .4) +
    
            ylim(0, ylimall) +
            ggtitle("") +
            theme_minimal() +
            theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.text = element_text(size=fsizeText)) +
            labs(x = "", y = "") +
            theme(axis.text.y=element_blank(),
                axis.ticks.y=element_blank())

p2_comb <- plot_grid(p2, p2_aux, rel_widths=c(0.875, .15))
pdf(paste0("All_Eff_Indonesia_T_", target, ".pdf"), height=5, width=7)
print(p2_comb)
dev.off()

df3 <- subset(df, Country=="Ethiopia")
df3_1 <- subset(df3, Model!='All')
df3_2 <- subset(df3, Model=='All')

p3 <- ggplot(data = df3_1, aes(x=Model, y=IRR)) +
        geom_bar(stat="identity", fill='lightskyblue', width=0.8) +
        
        geom_errorbar(aes(ymax = IRR_hi, ymin = IRR_low), position='dodge', width = .3) +
        
        ylim(0, ylimall) +
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
        theme_minimal() +
        theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.text = element_text(size=fsizeText), plot.title=element_text(size=fsizeTit)) +
        labs(x = "", y = "") +
        ggtitle("C: Ethiopia")

df3_2$Model <- c("Bayesian")
p3_aux <- ggplot(data = df3_2, aes(x=Model, y=IRR)) +
            geom_bar(stat="identity", fill='wheat', width=1.2) +

            geom_errorbar(aes(ymax = IRR_hi, ymin = IRR_low), position='dodge', width = .4) +
    
            ylim(0, ylimall) +
            ggtitle("") +
            theme_minimal() +
            theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.text = element_text(size=fsizeText)) +
            labs(x = "", y = "") +
            theme(axis.text.y=element_blank(),
                axis.ticks.y=element_blank())

p3_comb <- plot_grid(p3, p3_aux, rel_widths=c(0.875, .15))
pdf(paste0("All_Eff_Ethiopia_T_", target, ".pdf"), height=5, width=7)
print(p3_comb)
dev.off()

#Adol
target <- 3

df <- read.table(paste0("Impact_all_countries_T_", target, ".txt"), header=TRUE)
df1 <- subset(df, Country=="India")
df1_1 <- subset(df1, Model!='All')
df1_2 <- subset(df1, Model=='All')

p4 <- ggplot(data = df1_1, aes(x=Model, y=IRR)) +
        geom_bar(stat="identity", fill='lightskyblue', width=0.8) +
        
        geom_errorbar(aes(ymax = IRR_hi, ymin = IRR_low), position='dodge', width = .3) +
        
        ylim(0, ylimall) +
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
        theme_minimal() +
        theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.text = element_text(size=fsizeText), plot.title=element_text(size=fsizeTit)) +
        labs(x = "", y = "") +
        ggtitle("D: India")

df1_2$Model <- c("Bayesian")
p4_aux <- ggplot(data = df1_2, aes(x=Model, y=IRR)) +
            geom_bar(stat="identity", fill='wheat', width=1.2) +

            geom_errorbar(aes(ymax = IRR_hi, ymin = IRR_low), position='dodge', width = .4) +
    
            ylim(0, ylimall) +
            ggtitle("") +
            theme_minimal() +
            theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.text = element_text(size=fsizeText)) +
            labs(x = "", y = "") +
            theme(axis.text.y=element_blank(),
                axis.ticks.y=element_blank())

p4_comb <- plot_grid(p4, p4_aux, rel_widths=c(0.875, .15))

pdf(paste0("All_Eff_India_T_", target, ".pdf"), height=5, width=7)
print(p1_comb)
dev.off()

df2 <- subset(df, Country=="Indonesia")
df2_1 <- subset(df2, Model!='All')
df2_2 <- subset(df2, Model=='All')

p5 <- ggplot(data = df2_1, aes(x=Model, y=IRR)) +
        geom_bar(stat="identity", fill='lightskyblue', width=0.8) +
        
        geom_errorbar(aes(ymax = IRR_hi, ymin = IRR_low), position='dodge', width = .3) +
        
        ylim(0, ylimall) +
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
        theme_minimal() +
        theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.text = element_text(size=fsizeText), plot.title=element_text(size=fsizeTit)) +
        labs(x = "", y = "") +
        ggtitle("E: Indonesia")

df2_2$Model <- c("Bayesian")
p5_aux <- ggplot(data = df2_2, aes(x=Model, y=IRR)) +
            geom_bar(stat="identity", fill='wheat', width=1.2) +

            geom_errorbar(aes(ymax = IRR_hi, ymin = IRR_low), position='dodge', width = .4) +
    
            ylim(0, ylimall) +
            ggtitle("") +
            theme_minimal() +
            theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.text = element_text(size=fsizeText)) +
            labs(x = "", y = "") +
            theme(axis.text.y=element_blank(),
                axis.ticks.y=element_blank())

p5_comb <- plot_grid(p5, p5_aux, rel_widths=c(0.875, .15))
pdf(paste0("All_Eff_Indonesia_T_", target, ".pdf"), height=5, width=7)
print(p2_comb)
dev.off()

df3 <- subset(df, Country=="Ethiopia")
df3_1 <- subset(df3, Model!='All')
df3_2 <- subset(df3, Model=='All')

p6 <- ggplot(data = df3_1, aes(x=Model, y=IRR)) +
        geom_bar(stat="identity", fill='lightskyblue', width=0.8) +
        
        geom_errorbar(aes(ymax = IRR_hi, ymin = IRR_low), position='dodge', width = .3) +
        
        ylim(0, ylimall) +
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
        theme_minimal() +
        theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.text = element_text(size=fsizeText), plot.title=element_text(size=fsizeTit)) +
        labs(x = "", y = "") +
        ggtitle("F: Ethiopia")

df3_2$Model <- c("Bayesian")
p6_aux <- ggplot(data = df3_2, aes(x=Model, y=IRR)) +
            geom_bar(stat="identity", fill='wheat', width=1.2) +

            geom_errorbar(aes(ymax = IRR_hi, ymin = IRR_low), position='dodge', width = .4) +
    
            ylim(0, ylimall) +
            ggtitle("") +
            theme_minimal() +
            theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.text = element_text(size=fsizeText)) +
            labs(x = "", y = "") +
            theme(axis.text.y=element_blank(),
                axis.ticks.y=element_blank())

p6_comb <- plot_grid(p6, p6_aux, rel_widths=c(0.875, .15))
pdf(paste0("All_Eff_Ethiopia_T_", target, ".pdf"), height=5, width=7)
print(p3_comb)
dev.off()

title0 <- ggdraw() + draw_label(paste0("Newborns"), fontface = 'bold', x = 0.01, hjust = 0, size=fsizeColnames)
p1_comb <- plot_grid(title0, p1_comb, ncol = 1, rel_heights = c(0.1, 1))

title3 <- ggdraw() + draw_label(paste0("Adolescents"), fontface = 'bold', x = 0.01, hjust = 0, size=fsizeColnames)
p4_comb <- plot_grid(title3, p4_comb, ncol = 1, rel_heights = c(0.1, 1))

#Ahora todos juntos con cowplot
p_all <- plot_grid(p1_comb, p4_comb, NULL, p2_comb, p5_comb, NULL, p3_comb, p6_comb, NULL, rel_widths=c(1, 1, .05, 1, 1, .05, 1, 1, .05), ncol=3)

#create common x and y labels
y.grob <- textGrob("IRR (%)", gp=gpar(fontface="bold", col="black", fontsize=fsizeTit), rot=90)
x.grob <- textGrob("Model", gp=gpar(fontface="bold", col="black", fontsize=fsizeTit))
    
#add to plot and save
pdf("../All_Eff_Combined.pdf", width = 15, height = 20)
grid.arrange(arrangeGrob(p_all, left = y.grob, bottom = x.grob))
dev.off()

