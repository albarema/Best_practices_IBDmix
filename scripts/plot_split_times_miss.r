#--------------------------------------------------------------------------------

quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(tidyverse))
quiet(library("ggpubr"))
quiet(library("glue"))
quiet(library(grid))
quiet(library(scales))
quiet(library(ggplot2))
library(RColorBrewer)
#--------------------------------------------------------------------------------

read_all <- function(mod,lod){
    print(mod)
    infile = paste(mod, "/ibdmix_results/acc-eurasia.modern-allsamples-allparams-", lod,".gz", sep="")
    d = read_tsv(infile,show_col_types = FALSE)
    split_time = as.numeric(str_remove(strsplit(mod, "_")[[1]][2], "k"))*1e3
    d$pop = sapply(d$name, function(x) strsplit(x, "_")[[1]][1])
    d$n =sapply(d$name, function(x) strsplit(x, "_")[[1]][2]) 
    d$class = "modern"
    d$age = 0
    return(d)
}

#--------------------------------------------------------------------------------
outpath="./"
if (!dir.exists(outpath)) dir.create(outpath)

minlods = c(1,2,3,4,5,6,10,30)
mods=c('modelA_1_45ky','modelA_1_20ky','modelA_1_15ky','modelA_1_10ky','modelA_1_5ky','modelA_1_2ky')
df <- c()
for (i in mods){
    for (j in minlods){
        d <- read_all(i, j)
        df = rbind(df, d)
    }
}

df$split_time = sapply(df$mod, function(x) str_remove(strsplit(x, "_")[[1]][3],"ky"))
df$splits <- paste("split at ", df$split_time, " ka")

df$underestimation =df$overlap_introgression - df$true_introgression 
df$overestimation =  df$detected_introgression - df$overlap_introgression

df <- df %>% mutate(underestimation_seg=(overlap_introgression*4e8-true_introgression*4e8)/(true_introgression*4e8), 
    overestimation_seg = (detected_introgression*4e8-overlap_introgression*4e8)/(true_introgression*4e8))

df1 <- df %>% group_by(mod, minsegl, minlod) %>% 
    summarize(mean_under = mean(underestimation_seg), median_under = median(underestimation_seg), max_under= max(underestimation_seg), min_under=min(underestimation_seg) ) %>% 
    mutate(stats = "Underestimation", tot_info=paste(mod, minsegl, minlod, sep="-"))

df2 <- df %>% group_by(mod, minsegl, minlod) %>% 
    summarize(mean_over = mean(overestimation_seg), median_over = median(overestimation_seg), max_over= max(overestimation_seg), min_over=min(overestimation_seg) ) %>% 
    mutate(stats = "Overestimation", tot_info=paste(mod, minsegl, minlod, sep="-"))

df_tot <- merge(df1, df2[,-c(1:3)], by.x=9, by.y=6) %>% as_tibble()
df_tot$split_time = sapply(df_tot$mod, function(x) str_remove(strsplit(x, "_")[[1]][3],"ky"))
df_tot$split_time = factor(df_tot$split_time, levels=unique(sort(as.numeric(df_tot$split_time))))

split_labs = c('2' = "split at 2ka", '20' = "split at 20ka",'5' = "split at 5ka",'15' = "split at 15ka",'45' = "split at 45ka")
lod_labs=c('3'='minium LOD cutoff 3','4'='minium LOD cutoff 4','5'='minium LOD cutoff 5','6'='minium LOD cutoff 6','10'='minium LOD cutoff 10','20'='minium LOD cutoff 20','30'='minium LOD cutoff 30','50'='minium LOD cutoff 50','80'='minium LOD cutoff 80')

df = df[df$popX !="all", ]
df[df$mod == "modelA_2ky",]$popX <- "eurasia.modern"

nea_labs = c(nea1_2 = "nea1_2 (Vindijia-like)", nea2_1="nea2_1 (Altai-like)")
seg_labs = c('50000'="50kb segment len",'1e+05'="100kb segment len")

df$split_time = factor(df$split_time, levels=unique(sort(as.numeric(df$split_time))))
#df$minsegl = factor(df$minsegl, levels=unique(sort(as.numeric(df$minsegl))))

# coloured the bars in the error bars

p <- ggplot(df_tot[!df_tot$minlod %in% c(1,2,6,30,60),], aes(fill=split_time, y=median_under, x=factor(minsegl))) +
    geom_errorbar(aes(ymin=min_under, ymax=max_under), position=position_dodge(width = 1), colour="deeppink")+
    geom_point(aes(y=median_under,colour=split_time, fill=split_time), size=1, shape=21, position=position_dodge(width = 1)) +
    geom_hline(yintercept = 0, linetype = 2, colour='grey') +
    geom_errorbar(aes(ymin=min_over, ymax=max_over),position=position_dodge(width = 1), colour="green")+
    geom_point(aes(y=median_over, colour=split_time, fill=split_time), size=1, shape=21,position=position_dodge(width = 1), show.legend=F) +
    facet_wrap(~minlod, nrow=3,labeller = labeller(split_timnes=split_labs, minlod=lod_labs))+
    scale_y_continuous(n.breaks=10, limits=c(-1,1))+
    scale_colour_manual(values=brewer.pal(9, "Blues")[-c(1:3)])+
    scale_fill_manual(values=brewer.pal(9, "Blues")[-c(1:3)])+
    guides(colour=guide_legend(title="Split time (kya)", override.aes = list(size=3)),fill=guide_legend(title="Split time (kya)"))+
    labs(y = "prop. archaic segment misestimated", x="minimum segment length cutoff")+ theme_bw()+
    theme(panel.grid.major=element_line(linetype = "dotted", size=0.4, colour="grey"),strip.background = element_blank(),
    strip.text=element_text(size=11), axis.text=element_text(size=11), axis.text.x=element_text(angle=90, hjust=.5), legend.position="right", legend.text=element_text(size=11)) 

ggsave(paste(outpath, "dif_introgression_eurasia_lod_20_cols.png", sep=""), p, units="in", width=10.5, height=7.5)

################# --------------------  FIXED -------------------- #################

lod_breaks <- c(2,4,6,10,30,60)
seg_breaks <- unique(as.numeric(df$minsegl))

p <- ggplot(df[df$minlod == 4 ,] , aes(fill=factor(split_time), color = factor(split_time),y=underestimation_seg, group=mod,x=as.numeric(minsegl))) +
    geom_hline(yintercept = 0, linetype = 1, colour='black', linewidth=0.3) +
    geom_point(aes(y=underestimation_seg), size=0.1, alpha=0.2,shape=21, position=position_dodge(width = 1), colour="deeppink", fill="deeppink") +
    #geom_smooth( aes(y=underestimation_seg,colour=split_time, fill=split_time),method = "lm",span = 0.8,formula = y ~ splines::bs(x, 3), se =TRUE,size=0.5) + 
    geom_smooth( aes(y=underestimation_seg),method = "loess",span = 0.4, se =TRUE, alpha=.1) + 
    geom_point(aes(y=overestimation_seg), size=0.1, alpha=0.2,shape=21, position=position_dodge(width = 1), show.legend=F, colour="green", fill="green",) +
    geom_smooth( aes(y=overestimation_seg),method = "loess",span = 0.4, se =TRUE, alpha=0.1) + 
   # geom_smooth(aes(y=overestimation_seg,colour=split_time, fill=split_time),method = "lm",span = 0.8,formula = y ~ splines::bs(x, 3), se =TRUE,size=0.5) + 
    scale_y_continuous(n.breaks=10, limits=c(-1,1))+
    scale_colour_manual(values=brewer.pal(9, "Blues")[-c(1:3)])+
    scale_fill_manual(values=brewer.pal(9, "Blues")[-c(1:3)])+
    scale_x_continuous(breaks=seg_breaks)+
    guides(colour=guide_legend(title="Split time (kya)", override.aes = list(size=2)),fill=guide_legend(title="Split time (kya)"))+
    labs(y = "prop. archaic segment misestimated", x="minimum segment length cutoff")+ 
    theme_bw()+ theme(panel.grid.major=element_line(linetype = 2,color="grey", size=0.4),strip.background = element_blank(),
    strip.text=element_text(size=10), axis.text=element_text(size=9), axis.title=element_text(size=8),legend.position= "none",axis.text.x=element_text(angle=90, hjust=.5, size=8)) 

ggsave(paste(outpath, "dif_introgression_eurasia_LOD4.png", sep=""), p, units="in", width=3.5, height=3)

p <- ggplot(df[df$minsegl == 50000 & df$minlod > 1,], aes(fill=split_time, color = factor(split_time),y=underestimation_seg, group=mod, x=as.numeric(minlod)))+
    geom_hline(yintercept = 0, linetype = 1, colour='black', size=0.3) +
    geom_point(aes(y=underestimation_seg), size=0.1, alpha=0.2,shape=21, position=position_dodge(width = 1), colour="deeppink", fill="deeppink") +
    geom_smooth( aes(y=underestimation_seg,colour=split_time, fill=split_time),method = "loess",span = 0.4, se =TRUE, alpha=.1) +
    geom_point(aes(y=overestimation_seg), size=0.1, alpha=0.2,shape=21,position=position_dodge(width = 1), colour="green", fill="green",show.legend=F) +
    geom_smooth( aes(y=overestimation_seg,colour=split_time, fill=split_time),method = "loess",span = 0.4, se =TRUE, alpha=.1) +
    scale_y_continuous(n.breaks=10, limits=c(-1,1))+
    scale_colour_manual(values=brewer.pal(9, "Blues")[-c(1:3)])+
    scale_fill_manual(values=brewer.pal(9, "Blues")[-c(1:3)])+
    scale_x_continuous(breaks=lod_breaks)+
    guides(colour=guide_legend(title="Split time (kya)", override.aes = list(size=2)),fill=guide_legend(title="Split time (kya)"))+
    labs(y = "prop. archaic segment misestimated", x="minimum LOD cutoff")+ 
    theme_bw()+ theme(panel.grid.major=element_line(linetype = 2,color="grey", size=0.4),strip.background = element_blank(),
    strip.text=element_text(size=10), axis.text=element_text(size=9), axis.title=element_text(size=8),legend.position= "none",axis.text.x=element_text(angle=90, hjust=.5, size=8)) 

ggsave(paste(outpath, "dif_introgression_eurasia_SEG50kb.png", sep=""), p, units="in", width=3.5, height=3)

p <- ggplot(df[df$minsegl == 50000 & df$minlod == 4,], aes(split_time, underestimation_seg)) +
    geom_violin()+
    stat_summary(fun=median, geom="point", size=2, color="deeppink")+
    geom_hline(yintercept = 0, linetype = 2, colour='grey') +
    geom_violin(aes(y=overestimation_seg))+
    stat_summary(aes(y=overestimation_seg),fun=median, geom="point", size=2, color="green")+
    scale_y_continuous(n.breaks=8, limits=c(-1,1))+
    labs(y = "prop. archaic segment misestimated", x= "split time (kya)")+
    theme_bw()+ theme( axis.text=element_text(size=9), axis.title=element_text(size=8.5)) 

ggsave(paste(outpath, "dif_introgression_segs_fixed.png", sep=""), p, units="in", width=3.5, height=3)

p <- ggplot(df[df$minsegl == 50000 & df$minlod == 4,], aes(split_time, true_introgression)) +
    geom_violin()+
    stat_summary(fun=median, geom="point", size=2, color="black")+
    geom_hline(yintercept = 0.0225, linetype = 2, colour='grey') +
    labs(y = "true archaic proportion", x= "split time (kya)")+
    theme_bw()+ theme( axis.text=element_text(size=9), axis.title=element_text(size=8.5)) 

ggsave(paste(outpath, "split_time_true_intro.png", sep=""), p, units="in", width=3.5, height=3)
