quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(tidyverse))
quiet(library(grid))
quiet(library(scales))
quiet(library(ggpubr))
#--------------------------------------------------------------------------------
get_summary <- function(tib){
    summ_o = tib %>% group_by(age_groups,minsegl) %>%
        summarize(mean = mean(overestimation_seg), median = median(overestimation_seg), lower = min(overestimation_seg) , upper= max(overestimation_seg)) %>% 
        mutate(stats = "over")
    summ_u <- tib %>%  group_by(age_groups,minsegl) %>%
        summarize(mean = mean(underestimation_seg), median = median(underestimation_seg), lower = min(underestimation_seg) , upper= max(underestimation_seg)) %>% 
        mutate(stats = "under")
    summ <- rbind(summ_o, summ_u)
    return(summ)
}   
#--------------------------------------------------------------------------------
lods = c(1,2,3,4,5,6,10,30,60)
tib=c()
for (i in lods){
    ttmp = read_tsv(paste("modelA_1_45ky/ibdmix_results/power-pop1-allsamples-allparams-", i,".gz", sep=""),show_col_types = FALSE)
    nodes = read_tsv("modelA_1_45ky/results/nodes.tsv", show_col_types = FALSE) %>% select(name, time)
    dtmp_full = merge(ttmp,nodes, by="name", all=TRUE) %>% as_tibble()
    tib = rbind(tib, dtmp_full)
}
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
fullpops <- c("pop1","pop1.modern", "pop1.ancient")
minlods = c(1,2,3,4,5,6,10,30,60)
ntimes = length(minlods)-1
pop2test = "pop1"
output_dir = "."

fullset = tib[tib$popX %in% fullpops & tib$minlod %in% minlods & tib$minsegl >=1000,] %>% distinct()

fullset$n =sapply(fullset$name, function(x) strsplit(x, "_")[[1]][2]) 


tab_age <- c()
for (i in unique(fullset$minlod)){
    partset = fullset[fullset$minlod == i & fullset$popX %in% c("pop1.ancient", "pop1.modern"),]
    partset$age_groups <- cut(partset$time, breaks= c(0,1,15e3,30e3,Inf),
    include.lowest = TRUE,
    labels=c('present-day', '1-15kya', '15-30kya','>30kya'))
    tab_age <- rbind(tab_age,partset)
}

tab_age_mixed <- c()
for (i in unique(fullset$minlod)){
    partset = fullset[fullset$minlod == i & fullset$popX  == "pop1",] 
    partset$age_groups <- cut(partset$time, breaks= c(0,1,15e3,30e3,Inf),
    include.lowest = TRUE,
    labels=c('present-day', '1-15kya', '15-30kya','>30kya'))
    tab_age_mixed <- rbind(tab_age_mixed,partset)
}

alls = rbind(tab_age_mixed, tab_age)

seg_breaks = unique(alls$minsegl)

alls$age_groups = factor(alls$age_groups, levels=c('present-day', '1-15kya', '15-30kya','>30kya'))

plot_summ =  alls[!alls$minlod %in% c(1,2),] %>% 
  group_by(popX,age_groups, minsegl, minlod) %>% 
  summarize(mean_f1 = mean(F1_wi_truepc), mean_nMCC = mean(nMCC_wi_truepc), mean_specificity = mean(TN/(TN+FP)), mean_sensitivity = mean(TP/(TP+FN)),mean_precision=mean(precision_wi_truepc), mean_f1seg = mean(F1_seg))


for (i in unique(plot_summ$popX)){
    plot_i = plot_summ[plot_summ$popX ==i,]
    p1 <- ggplot(plot_i, aes(x=factor(minsegl/1e3), fill=mean_f1, y=factor(minlod)))+ 
    geom_tile(color="white") +  
    facet_wrap(~age_groups, nrow=3)+ 
    labs(y = "minimum LOD score", x="minimum segment length [Kb]")+
    scale_fill_gradientn(colors=c("navy","white","orange","red"),values=rescale(c(0,0.5,0.75,1)),breaks=c(0,0.5,0.75,1),limits=c(0,1), oob=squish)+
    theme_bw()+ theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linetype = "dotted"),
    strip.background = element_blank(),strip.text=element_text(size=9), 
    axis.text.x =element_text(angle=90), legend.position="bottom")+
    guides(fill=guide_colourbar(title.position = "left", title.hjust = -1.5, title.vjust=1, barwidth = 7, barheight = 0.4,ticks.colour= "black", ticks.linewidth = 3,label=T,title = "F1"))

    ggsave(paste(output_dir, "/heatmap_F1_introgression_clusters_", i, ".png", sep=""), p1, units="in", width=7, height=7)
}

# PLOT FOR POP1 VS POP.ANCIENT VS.POP1.MODERN LOOK THE SAME 
i = "pop1"
plot_i = plot_summ[plot_summ$popX ==i,]

    p1 <- ggplot(plot_i, aes(x=factor(minsegl/1e3), fill=mean_specificity, y=factor(minlod)))+ 
    geom_tile(color="white") +  
    facet_wrap(~age_groups, nrow=3)+ 
    labs(y = "minimum LOD score", x="minimum segment length [Kb]")+
    scale_fill_gradientn(colors=c("navy","white","orange"),values=rescale(c(0.99,0.995,1)),breaks=c(0.99,0.995,1),limits=c(0.99,1), oob=squish)+
    theme_bw()+ theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linetype = "dotted"),
    strip.background = element_blank(),strip.text=element_text(size=9), 
    axis.text.x =element_text(angle=90), legend.position="bottom")+
    guides(fill=guide_colourbar(title.position = "left", title.hjust = -1.5, title.vjust=1, barwidth = 15, barheight = 0.4,ticks.colour= "black", ticks.linewidth = 3,label=T,title = "Specificity"))

    ggsave(paste(output_dir, "/heatmap_specificity_introgression_clusters_", i, ".png", sep=""), p1, units="in", width=7, height=7)

    p1 <- ggplot(plot_i, aes(x=factor(minsegl/1e3), fill=mean_sensitivity, y=factor(minlod)))+ 
    geom_tile(color="white") +  
    facet_wrap(~age_groups, nrow=3)+ 
    labs(y = "minimum LOD score", x="minimum segment length [Kb]")+
    scale_fill_gradientn(colors=c("navy","white","orange","red"),values=rescale(c(0,0.8,0.9,1)),breaks=c(0,0.8,0.9,1),limits=c(0,1), oob=squish)+
    theme_bw()+ theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linetype = "dotted"),
    strip.background = element_blank(),strip.text=element_text(size=9), 
    axis.text.x =element_text(angle=90), legend.position="bottom")+
    guides(fill=guide_colourbar(title.position = "left", title.hjust = -1.5, title.vjust=1, barwidth = 13, barheight = 0.4,ticks.colour= "black", ticks.linewidth = 3,label=T,title = "Sensitivity"))

    ggsave(paste(output_dir, "/heatmap_sensitivity_introgression_clusters_", i, ".png", sep=""), p1, units="in", width=7, height=7)

    p1 <- ggplot(plot_i, aes(x=factor(minsegl/1e3), fill=mean_nMCC, y=factor(minlod)))+ 
    geom_tile(color="white") +  
    facet_wrap(~age_groups, nrow=3)+ 
    labs(y = "minimum LOD score", x="minimum segment length [Kb]")+
    scale_fill_gradientn(colors=c("navy","white","red"),values=rescale(c(0.5,0.9,1)),breaks=c(0.5,0.9,1),limits=c(0.5,1), oob=squish)+
    theme_bw()+ theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linetype = "dotted"),
    strip.background = element_blank(),strip.text=element_text(size=9), 
    axis.text.x =element_text(angle=90), legend.position="bottom")+
    guides(fill=guide_colourbar(title.position = "left", title.hjust = -1.5, title.vjust=1, barwidth = 7, barheight = 0.4,ticks.colour= "black", ticks.linewidth = 3,label=T,title = "nMCC"))

    ggsave(paste(output_dir, "/heatmap_nMCC_introgression_clusters_", i, ".png", sep=""), p1, units="in", width=7, height=7)


best_F1_seg = plot_i %>% filter(FP <= 0.1) %>% group_by(age_groups) %>% top_n(1, wt =mean_f1seg)

    p1 <- ggplot(plot_i, aes(x=factor(minsegl/1e3), fill=mean_f1, y=factor(minlod)))+ 
    geom_tile(color="white") +  
    facet_wrap(~age_groups, nrow=3)+ 
    labs(y = "minimum LOD score", x="minimum segment length [Kb]")+
    scale_fill_gradientn(colors=c("navy","white","orange","red"),values=rescale(c(0,0.8,0.9,1)),breaks=c(0,0.8,0.9,1),limits=c(0,1), oob=squish)+
    theme_bw()+ theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linetype = "dotted"),
    strip.background = element_blank(),strip.text=element_text(size=9), 
    axis.text.x =element_text(angle=90), legend.position="bottom")+
    guides(fill=guide_colourbar(title.position = "left", title.hjust = -1.5, title.vjust=1, barwidth = 13, barheight = 0.4,ticks.colour= "black", ticks.linewidth = 3,label=T,title = "F1"))

    ggsave(paste(output_dir, "/heatmap_F1_introgression_clusters_", i, ".png", sep=""), p1, units="in", width=7, height=7)


    p1 <- ggplot(plot_i, aes(x=factor(minsegl/1e3), fill=mean_precision, y=factor(minlod)))+ 
    geom_tile(color="white") +  
    facet_wrap(~age_groups, nrow=3)+ 
    labs(y = "minimum LOD score", x="minimum segment length [Kb]")+
    scale_fill_gradientn(colors=c("navy","white","orange","red"),values=rescale(c(0.9,0.95,0.98,1)),breaks=c(0.9,0.95,0.98,1),limits=c(0.9,1), oob=squish)+
    theme_bw()+ theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linetype = "dotted"),
    strip.background = element_blank(),strip.text=element_text(size=9), 
    axis.text.x =element_text(angle=90), legend.position="bottom")+
    guides(fill=guide_colourbar(title.position = "left", title.hjust = -1.5, title.vjust=1, barwidth = 15, barheight = 0.4,ticks.colour= "black", ticks.linewidth = 3,label=T,title = "Precision"))

    ggsave(paste(output_dir, "/heatmap_precison_introgression_clusters_", i, ".png", sep=""), p1, units="in", width=7, height=7)
