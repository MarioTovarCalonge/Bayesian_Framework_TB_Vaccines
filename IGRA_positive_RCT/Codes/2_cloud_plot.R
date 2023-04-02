#!/usr/bin/Rscript


library(ggplot2)
library(grid)
library(gridExtra)
library(viridis)
library(cowplot)

fitG = function(x,y,mu,sig,scale){
    f = function(p){
        d = p[3]*dnorm(x,mean=p[1],sd=p[2])
        sum((d-y)^2)
    }
    optim(c(mu,sig,scale),f)
}

gauss_pdf <- function(x){
    return(dnorm(x, param[1], param[2]))
}

scale_par <- c(.01, .001, 1, .1, .01, .1, .1)

#-------------------------------------------------------------------------------NUBE-----------------------------------------------------------------

for(mode_gen in 1:7){
    dat <- read.table(paste0("Outputs/Model_", sprintf("%d", mode_gen), "/test_sto_A_Mixed_pop.txt"), header=FALSE, skip=1)
    M <- as.matrix(dat)
    eps_axix <- seq(0, 0.995, 0.005)
    AX <- NULL
    AY <- NULL
    for(i in 1:length(eps_axix)){
        AX <- c(AX, rep(eps_axix[i], length(M[, i])))
        AY <- c(AY, M[,i])
    }

    newDF <- data.frame(x=AX, y=AY)

    plot1 <- ggplot(data=newDF, aes(x=x, y=y) ) +
            geom_bin2d(bins = 100) +
            scale_fill_continuous(type = "viridis") +
            coord_cartesian(c(0, 1), c(-50, 101)) +
            theme_bw() +
            labs(x = expression(epsilon), y = expression(paste("V", E["dis"], " (%)")))
    
    legend <- cowplot::get_legend(plot1) 
    plot1<- plot1 + theme(legend.position = "none", axis.title.x=element_text(size=12))
    
    d <- density(newDF$y)
    max_val = d$y[which.max(d$y)]
    print(max_val)
    max_val = round(max_val, 3)
    print(max_val)
    by_val = max_val/2
    bk = round(seq(0, max_val, by_val), 3)
    
    newDF2 <- read.table(paste0("Outputs/Model_", sprintf("%d", mode_gen), "/P(VE_dis|eps)_Mixed_pop.txt"), header=TRUE)
    eps <- as.vector(unlist(newDF2$eps))
    pr <- as.vector(unlist(newDF2$prob))
    
    suppressWarnings(param <- fitG(eps, pr, 0.5, 0.2, scale_par[mode_gen])$par)
    print(param)
    if(param[1] < 0 || param[2] < 0){
        param[1] = 0
        param[2] = 1
    }
    
    const <- (1.0 - integrate(gauss_pdf, 1, Inf)$value)
    
    gauss_pdf_2 <- function(x){
        return(dnorm(x, param[1], param[2])/const)
    }
    
    df_plot <- data.frame(EPS = eps, peps = pr, peps_fit = gauss_pdf_2(eps))
    
    p2 <- ggplot(df_plot, aes(x = EPS, y=peps_fit)) + geom_line(color="grey20") + geom_area(mapping = aes(x = ifelse(EPS>0 & EPS< 1 , EPS, 0)), fill = "grey20") + coord_cartesian(c(0, 1)) + theme_bw() + 
            labs(x = expression(epsilon), y = expression(paste("P(", epsilon, "|V", E["dis"], ")"))) + theme(axis.title.x=element_text(size=12))
    
    p3 <- ggplot(newDF, aes(x = y)) + stat_density() + coord_flip(c(-50, 101))  + theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust=0.7)) + scale_y_continuous(breaks = bk) + labs(x = expression(paste("V", E["dis"], " (%)")))
    
    upper_row <- plot_grid(p2, legend, labels = c('A', ''), label_size = 12,  rel_widths = c(3, 0.8))
    bottom_row <- plot_grid(plot1, p3, labels = c('B', 'C'), label_size = 12,  rel_widths = c(3, 0.8))
    p <- plot_grid(upper_row, bottom_row, ncol = 1, rel_heights=c(0.8, 3))

    ggsave(paste0("Outputs/Model_", sprintf("%d", mode_gen), "/Cloud_2D_hist_", sprintf("%d", mode_gen), "_cowplot.pdf"), p)
    
}


















    
    

