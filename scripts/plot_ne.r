quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(tidyverse))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

filename = args[1]
outf = args[2]
outf2 = args[3]
outf3 = args[4]
outf4 = args[5]
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#filename = "all_stats_per_sample_size.txt"

print("reading file 1...")
d = read_tsv(filename,show_col_types = FALSE)
d$samplesize[d$samplesize == 100] <- "full set"
scaleFUN <- function(x) sprintf("%.2f", x)
d$samplesize = factor(d$samplesize, levels=c(sort(as.numeric(unique(d$samplesize))), "full set"))

d$model_long = d$model
d$model_long[d$model_long == "modelA_1_5000_45ky"] <- "modelA (Ne=5,000)"
d$model_long[d$model_long == "modelA_1_10000_45ky"] <- "modelB (Ne=10,000)"
d$model_long[d$model_long == "modelA_1_2500_45ky"] <- "modelC (Ne=2,500)"

d$ymin = 0; d$ymax = 1
d$stats[d$stats == "overestimation"] <- "Average overestimation"
d$stats[d$stats == "underestimation"] <- "Average underestimation"
d$stats[d$stats == "F1_wi_truepc"] <- "F1_wi"
d$stats[d$stats == "nMCC_wi_truepc"] <- "nMCC_wi"

d2 = d[d$stats %in% c('Average underestimation','Average overestimation'),]

d = d[d$stats != "Average overestimation",]
d$ymax[d$stats == "Average underestimation"]  = 0.00
d$ymin[d$stats == "F1_wi"]  = 0
d$ymin[d$stats == "nMCC_wi"]  = 0.4

d$stats_fixed = factor(d$stats, levels=c("Average overestimation", "F1_wi", "Average underestimation", "nMCC_wi"))


hlines =tibble(stats=unique(d$stats), vals=c(0,NA,0.5))

print("plotting ...")
#mean contains the misestimation (estimated - true)

color_values <- c("modelA (Ne=5,000)" = "red", "modelB (Ne=10,000)" = "gold", "modelC (Ne=2,500)" = "purple")

p1 <- ggplot(d, aes(x=samplesize, y=mean,group=1))+
    geom_hline(data=hlines,aes(yintercept=vals), colour="black")+
    geom_line(data=d[d$model == "modelA_1_5000_45ky",],aes(colour=model_long),linetype=1, size=0.5)+
    geom_line(data=d[d$model == "modelA_1_10000_45ky",],aes(colour=model_long),linetype=1, size=0.5)+ 
    geom_line(data=d[d$model == "modelA_1_2500_45ky",],aes(colour=model_long),linetype=1, size=0.5)+ 
    geom_ribbon(data=d[d$model == "modelA_1_5000_45ky",],aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = model_long), alpha = .2 ) +
    geom_ribbon(data=d[d$model == "modelA_1_10000_45ky",],aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = model_long), alpha = .3  ) +
    geom_ribbon(data=d[d$model == "modelA_1_2500_45ky",],aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = model_long), alpha = .1 ) +
    geom_line(data=d[d$model == "modelA_1_5000_45ky",],aes(colour=model_long,y=median),linetype="dotted", size=0.5)+
    geom_line(data=d[d$model == "modelA_1_10000_45ky",],aes(colour=model_long,y=median),linetype="dotted", size=0.5)+
    geom_line(data=d[d$model == "modelA_1_2500_45ky",],aes(colour=model_long, y=median),linetype="dotted", size=0.5)+
    #ylim(c(-0.01,0.02))+
    labs(x="sample size", title="")+
    scale_colour_manual(values=color_values)+ 
    scale_fill_manual(values=color_values)+ 
    theme_bw()+
    facet_wrap(~stats, ncol=3,scales="free_y")+
    #guides(color = guide_legend(nrow = 2, byrow = TRUE))+
    theme(axis.text.y = element_text(size=11),axis.text.x = element_text(size=11,angle = 90, vjust = 0.5, hjust=1),axis.title = element_text(size=10),axis.title.y= element_blank(),
        strip.background = element_blank(),strip.text=element_text(size=10),
        legend.title=element_blank(), legend.position ="bottom",legend.text=element_text(size=10)) +
    scale_y_continuous(expand = c(0,0), labels=scaleFUN) +
    geom_blank(aes(y = ymin)) + geom_blank(aes(y = ymax))

ggsave(outf, p1,units="in", height=3, width=9)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
d1 = d[d$stats %in% c('F1_wi','nMCC_wi'),]
stats_labs = c('F1_wi'="F1", "nMCC_wi"="nMCC")
hlines =tibble(stats=unique(d1$stats), vals=c(NA,0.5))

p1 <- ggplot(d1, aes(x=samplesize, y=mean,group=1))+
    geom_hline(data=hlines,aes(yintercept=vals), colour="black")+
    geom_line(data=d1[d1$model == "modelA_1_5000_45ky",],aes(colour=model_long),linetype=1, size=0.5)+
    geom_line(data=d1[d1$model == "modelA_1_10000_45ky",],aes(colour=model_long),linetype=1, size=0.5)+ 
    geom_line(data=d1[d1$model == "modelA_1_2500_45ky",],aes(colour=model_long),linetype=1, size=0.5)+ 
    geom_ribbon(data=d1[d1$model == "modelA_1_5000_45ky",],aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = model_long), alpha = .2 ) +
    geom_ribbon(data=d1[d1$model == "modelA_1_10000_45ky",],aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = model_long), alpha = .3  ) +
    geom_ribbon(data=d1[d1$model == "modelA_1_2500_45ky",],aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = model_long), alpha = .1 ) +
    geom_line(data=d1[d1$model == "modelA_1_5000_45ky",],aes(colour=model_long,y=median),linetype="dotted", size=0.5)+
    geom_line(data=d1[d1$model == "modelA_1_10000_45ky",],aes(colour=model_long,y=median),linetype="dotted", size=0.5)+
    geom_line(data=d1[d1$model == "modelA_1_2500_45ky",],aes(colour=model_long, y=median),linetype="dotted", size=0.5)+
    #ylim(c(-0.01,0.02))+
    labs(x="", title="")+
    scale_colour_manual(values=color_values)+ 
    scale_fill_manual(values=color_values)+ 
    theme_bw()+
    facet_wrap(~stats, ncol=2,scales="free_y", labeller =labeller(stats=stats_labs))+
    #guides(color = guide_legend(nrow = 2, byrow = TRUE))+
    theme(axis.text.y = element_text(size=9),axis.text.x = element_text(size=9,angle = 90, vjust = 0.5, hjust=1),axis.title = element_text(size=9),axis.title.y= element_blank(),
        strip.background = element_blank(),strip.text=element_text(size=10),
        legend.title=element_blank(), legend.position ="none") +
    scale_y_continuous(expand = c(0,0)) +
    geom_blank(aes(y = ymin)) + geom_blank(aes(y = ymax))

ggsave(outf2, p1,units="in", height=3, width=6)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

hlines =tibble(stats=unique(d1$stats), vals=c(0,0))

print("plotting 3/4...")

p1 <- ggplot(d2, aes(x=samplesize, y=mean,group=1))+
    geom_hline(data=hlines,aes(yintercept=vals), colour="black")+
    geom_line(data=d2[d2$model == "modelA_1_5000_45ky" ,],aes(colour=model_long, group=stats),linetype=1, size=0.5)+
    geom_line(data=d2[d2$model == "modelA_1_10000_45ky",],aes(colour=model_long, group=stats),linetype=1, size=0.5)+ 
    geom_line(data=d2[d2$model == "modelA_1_2500_45ky",],aes(colour=model_long, group=stats),linetype=1, size=0.5)+ 
    geom_ribbon(data=d2[d2$model == "modelA_1_5000_45ky",],aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = model_long, group=stats), alpha = .2 ) +
    geom_ribbon(data=d2[d2$model == "modelA_1_10000_45ky",],aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = model_long, group=stats), alpha = .3  ) +
    geom_ribbon(data=d2[d2$model == "modelA_1_2500_45ky",],aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = model_long, group=stats), alpha = .1 ) +
    geom_line(data=d2[d2$model == "modelA_1_5000_45ky",],aes(colour=model_long,y=median, group=stats),linetype="dotted", size=0.5)+
    geom_line(data=d2[d2$model == "modelA_1_10000_45ky",],aes(colour=model_long,y=median, group=stats),linetype="dotted", size=0.5)+
    geom_line(data=d2[d2$model == "modelA_1_2500_45ky",],aes(colour=model_long, y=median, group=stats),linetype="dotted", size=0.5)+
    labs(x="", y = 'pct. of genome',title="Average archaic misestimation")+
    scale_colour_manual(values=color_values)+ 
    scale_fill_manual(values=color_values)+ 
    theme_bw()+
    #guides(color = guide_legend(nrow = 2, byrow = TRUE))+
    theme(axis.text.y = element_text(size=9),axis.text.x = element_text(size=9,angle = 90, vjust = 0.5, hjust=1),plot.title = element_text(size=9,hjust = 0.5),
        strip.background = element_blank(),strip.text=element_text(size=8),
        legend.title=element_blank(), legend.position ="none") 

ggsave(outf3, p1,units="in", height=2.7, width=3.1)


hline <- tibble(stats="miss", vals=0)
p1 <- ggplot(d2, aes(x=samplesize, y=mean,group=1))+
    geom_line(data=d2[d2$model == "modelA_1_10000_45ky",],aes(colour=model_long, group=stats),linetype=1, size=0.5)+ 
    geom_line(data=d2[d2$model == "modelA_1_5000_45ky" ,],aes(colour=model_long, group=stats),linetype=1, size=0.5)+
    geom_line(data=d2[d2$model == "modelA_1_2500_45ky",],aes(colour=model_long, group=stats),linetype=1, size=0.5)+ 
    geom_ribbon(data=d2[d2$model == "modelA_1_10000_45ky",],aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = model_long, group=stats), alpha = .3  ) +
    geom_ribbon(data=d2[d2$model == "modelA_1_5000_45ky",],aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = model_long, group=stats), alpha = .2 ) +
    geom_ribbon(data=d2[d2$model == "modelA_1_2500_45ky",],aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = model_long, group=stats), alpha = .1 ) +
    geom_line(data=d2[d2$model == "modelA_1_10000_45ky",],aes(colour=model_long,y=median, group=stats),linetype="dotted", size=0.5)+
    geom_line(data=d2[d2$model == "modelA_1_5000_45ky",],aes(colour=model_long,y=median, group=stats),linetype="dotted", size=0.5)+
    geom_line(data=d2[d2$model == "modelA_1_2500_45ky",],aes(colour=model_long, y=median, group=stats),linetype="dotted", size=0.5)+
#    scale_y_continuous(expand = c(0,0),limits = c(0, 1))+
    labs(x="", y = 'prop. introgressed segments')+
    scale_colour_manual(values=color_values)+ 
    scale_fill_manual(values=color_values)+ 
    facet_wrap(~stats, ncol=1, scales="free")+
    theme_bw()+
    #guides(color = guide_legend(nrow = 2, byrow = TRUE))+
    theme(axis.text.y = element_text(size=9),axis.text.x = element_text(size=9,angle = 90, vjust = 0.5, hjust=1),plot.title = element_text(size=9,hjust = 0.5),
        axis.title.y=element_text(size=9),panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),strip.text=element_text(size=8),
        legend.title=element_blank(), legend.position ="none") 

ggsave(outf4, p1,units="in", height=3.5, width=3.1)

