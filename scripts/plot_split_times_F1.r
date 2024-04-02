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
read_eurasia <- function(mod, lod){
    print(paste(mod, lod, sep=" "))
    infile = paste(mod, "/ibdmix_results/power-eurasia.modern-allsamples-allparams-", lod,".gz", sep="")
    df = read_tsv(infile,show_col_types = FALSE)
    return(df)
}

#--------------------------------------------------------------------------------
outpath="./"
if (!dir.exists(outpath)) dir.create(outpath)

nea_labs = c(nea1_2 = "nea1_2 (Vindijia-like)", nea2_1="nea2_1 (Altai-like)")
seg_labs = c('50000'="50kb segment len",'1e+05'="100kb segment len")

minlods = c(3,4,5,6,10,30,60)
mods=c('modelA_1_45ky','modelA_1_20ky','modelA_1_15ky','modelA_1_10ky','modelA_1_5ky','modelA_1_2ky')
df <- c()
for (i in mods){
    for (j in minlods){
        d <- read_eurasia(i, j)
        df = rbind(df, d)
    }
}

df$split_time = sapply(df$mod, function(x) str_remove(strsplit(x, "_")[[1]][3],"ky"))
df$splits <- paste("split at ", df$split_time, " ka")

summs_f <- df %>% group_by(mod, minsegl, minlod) %>% 
    summarize(mean = mean(F1_wi_truepc), median = median(F1_wi_truepc), sd = sd(F1_wi_truepc) ) %>% 
    mutate(stats = "Average F1")

summs_f$split_time = sapply(summs_f$mod, function(x) as.numeric(str_remove(strsplit(x, "_")[[1]][3],"ky")))
best_f <- summs_f %>% group_by(mod) %>% filter(abs(median - 1) == min(abs(median - 1))) %>% arrange(desc(minlod), desc(minsegl))
best_f_mean <- summs_f %>% group_by(mod) %>% filter(abs(mean - 1) == min(abs(mean - 1))) %>% arrange(desc(minlod), desc(minsegl))%>% distinct(mod,.keep_all= TRUE)


summs_m <- df %>% group_by(mod, minsegl, minlod) %>% 
    summarize(mean = mean(nMCC_wi_truepc), median = median(nMCC_wi_truepc), sd = sd(nMCC_wi_truepc) ) %>% 
    mutate(stats = "Average nMCC")
summs_m$split_time = sapply(summs_m$mod, function(x) as.numeric(str_remove(strsplit(x, "_")[[1]][3],"ky")))
best_m <- summs_m %>% group_by(mod) %>% filter(abs(median - 1) == min(abs(median - 1))) %>% arrange(desc(minlod), desc(minsegl)) %>% distinct(mod,.keep_all= TRUE)


summs_p <- df %>% group_by(mod, minsegl, minlod) %>% 
    summarize(mean = mean(precision_wi_truepc), median = median(precision_wi_truepc), sd = sd(precision_wi_truepc) ) %>% 
    mutate(stats = "Average precision")
summs_p$split_time = sapply(summs_p$mod, function(x) as.numeric(str_remove(strsplit(x, "_")[[1]][3],"ky")))
best_p <- summs_p %>% group_by(mod) %>% filter(abs(median - 1) == min(abs(median - 1))) %>% arrange(desc(minlod), desc(minsegl)) %>% distinct(mod,.keep_all= TRUE)


summs_r <- df %>% group_by(mod, minsegl, minlod) %>% 
    summarize(mean = mean(recall_wi_truepc), median = median(recall_wi_truepc), sd = sd(recall_wi_truepc) ) %>% 
    mutate(stats = "Average recall")
summs_r$split_time = sapply(summs_r$mod, function(x) as.numeric(str_remove(strsplit(x, "_")[[1]][3],"ky")))
best_r <- summs_r %>% group_by(mod) %>% filter(abs(median - 1) == min(abs(median - 1))) %>% arrange(desc(minlod),desc(minsegl)) %>% distinct(mod,.keep_all= TRUE)


summs_acc <- df %>% group_by(mod, minsegl, minlod) %>% 
    summarize(mean = mean(acc_wi_truepc), median = median(acc_wi_truepc), sd = sd(acc_wi_truepc) ) %>% 
    mutate(stats = "Average accuracy")
summs_acc$split_time = sapply(summs_acc$mod, function(x) as.numeric(str_remove(strsplit(x, "_")[[1]][3],"ky")))
best_acc <- summs_acc %>% group_by(mod) %>% filter(abs(median - 1) == min(abs(median - 1))) %>% arrange(desc(minlod), minsegl) %>% distinct(mod,.keep_all= TRUE)



summs_fs <- df %>% group_by(mod, minsegl, minlod) %>% 
    summarize(mean = mean(F1_seg), median = median(F1_seg), sd = sd(F1_seg) ) %>% 
    mutate(stats = "Average F1 segments")
summs_fs$split_time = sapply(summs_f$mod, function(x) as.numeric(str_remove(strsplit(x, "_")[[1]][3],"ky")))

best_fs_mean <- summs_fs %>% group_by(mod) %>% filter(abs(mean - 1) == min(abs(mean - 1))) %>% arrange(desc(minlod), desc(minsegl)) %>% distinct(mod,.keep_all= TRUE)
best_fs <- summs_fs %>% group_by(mod)  %>% filter(abs(median - 1) == min(abs(median - 1))) %>% arrange(desc(minlod), desc(minsegl)) %>% distinct(mod,.keep_all= TRUE)


################# --------------------  ALL PARAMS -------------------- #################

eu <- ggplot(summs_f, aes(x=factor(minlod), factor(minsegl/1000), fill=mean)) +
    geom_tile(color="white")+
    labs(y="minimum segment length [Kb]", x="minimum LOD score")+
#    geom_tile(data=best_f,colour="black", size =0.3) +
    scale_fill_gradientn(colors=c("navy","white","orange","red"),values=rescale(c(0,0.8,0.9,1)),breaks=c(0,0.8,0.9,1),limits=c(0,1), oob=squish)+
    facet_wrap(~split_time, ncol=3)+
    theme_bw()+ theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linetype = "dotted"),strip.background = element_blank(),
    strip.text=element_text(size=9), axis.text=element_text(size=9),panel.spacing = unit(0.2, "lines")) +
    guides(fill=guide_colourbar(barheight = 13,ticks.colour= "black", ticks.linewidth = 2,label=T,title = "F1"))


ggsave(paste(outpath, "F1_introgression_heatmap.png", sep=""), eu, units="in", width=7, height=5)


eu <- ggplot(summs_fs, aes(x=factor(minlod), factor(minsegl/1000), fill=mean)) +
    geom_tile(color="white")+
#    geom_tile(data=best_fs,colour="black", size =0.3) +
    labs(y="minimum segment length [Kb]", x="minimum LOD score")+
    scale_fill_gradientn(colors=c("navy","white","orange","red"),values=rescale(c(0,0.8,0.9,1)),breaks=c(0,0.8,0.9,1),limits=c(0,1), oob=squish)+
    facet_wrap(~split_time, ncol=3)+
    theme_bw()+ theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linetype = "dotted"),strip.background = element_blank(),
    strip.text=element_text(size=9), axis.text=element_text(size=9),panel.spacing = unit(0.2, "lines")) +
    guides(fill=guide_colourbar(barheight = 13,ticks.colour= "black", ticks.linewidth = 2,label=T,title = "F1"))

ggsave(paste(outpath, "F1_seg_introgression_heatmap.png", sep=""), eu, units="in",width=7, height=5)

################# --------------------  FIXED EITHER -------------------- #################
lod_breaks <- c(2,4,6,10,30)
df$split_time <- factor(df$split_time, levels = unique(sort(as.numeric(df$split_time))))
cols_blue = brewer.pal(9, "Blues")[-c(1:3)]
names(cols_blue) = levels(df$split_time)

eu <- ggplot(df[df$minsegl == 50000 & df$minlod <= 30,], aes(y=F1_wi_truepc, x=minlod,group=mod,color = factor(split_time), fill=factor(split_time))) +
    labs(x="minimum LOD cutoff", y= "F1")+
    geom_point(alpha=0.2, size=0.2, show.legend=F) +
    geom_smooth( method = "lm",span = 0.8,formula = y ~ splines::bs(x, 3), se =TRUE) + 
    scale_colour_manual(name="split time (kya)", values=brewer.pal(9, "Blues")[-c(1:3)])+
    scale_fill_manual(name="split time (kya)", values=brewer.pal(9, "Blues")[-c(1:3)])+
    scale_y_continuous(n.breaks=5, limits=c(0,1))+
    scale_x_continuous(breaks=lod_breaks)+
    theme_bw()+ theme(panel.grid.major=element_line(linetype = 2,color="grey", size=0.4),strip.background = element_blank(),
    strip.text=element_text(size=10), axis.text=element_text(size=9), legend.position= "none",axis.text.x=element_text(angle=90, hjust=.5, size=8)) 

ggsave(paste(outpath, "F1_introgression_eu_50kb.png", sep=""), eu, units="in", width=3.5, height=3)


eu <- ggplot(df[df$minsegl == 50000 & df$minlod <= 30,], aes(y=nMCC_wi_truepc, x=minlod,group=mod,color = split_time, fill= split_time)) +
    labs(x="minimum LOD cutoff", y= "nMCC")+
    geom_point(alpha=0.2, size=0.2, show.legend=F) +
    geom_smooth( method = "lm",span = 0.8,formula = y ~ splines::bs(x, 3), se =TRUE) + 
    scale_colour_manual(name="split time (kya)", values=cols_blue)+
    scale_fill_manual(name="split time (kya)", values=cols_blue)+
    scale_y_continuous(n.breaks=5, limits=c(0.5,1))+
    scale_x_continuous(breaks=lod_breaks)+
    theme_bw()+ theme(panel.grid.major=element_line(linetype = 2,color="grey", size=0.4),strip.background = element_blank(),
    strip.text=element_text(size=10), axis.text=element_text(size=9), legend.position= "none",axis.text.x=element_text(angle=90, hjust=.5, size=8)) 

ggsave(paste(outpath, "MCC_introgression_eu_50kb.png", sep=""), eu, units="in", width=3.5, height=3)

segl_break = unique(df$minsegl)
eu <- ggplot(df[df$minlod == 4,], aes(y=F1_wi_truepc, x=minsegl,group=mod,color = split_time, fill= split_time)) +
    labs(x="min. segment length cutoff", y= "F1")+
    geom_point(alpha=0.2, size=0.2, show.legend=F) +
    geom_smooth( method = "lm",span = 0.8,formula = y ~ splines::bs(x, 3), se =TRUE) + 
    scale_colour_manual(name="split time (kya)", values=brewer.pal(9, "Blues")[-c(1:3)])+
    scale_fill_manual(name="split time (kya)", values=brewer.pal(9, "Blues")[-c(1:3)])+
    scale_y_continuous(n.breaks=5, limits=c(0,1))+
    scale_x_continuous(breaks=segl_break)+
    theme_bw()+ theme(panel.grid.major=element_line(linetype = 2,color="grey", size=0.4),strip.background = element_blank(),
    strip.text=element_text(size=10), axis.text=element_text(size=9), legend.position= "none",axis.text.x=element_text(angle=90, hjust=.5, size=8)) 

ggsave(paste(outpath, "F1_introgression_eu_4lod.png", sep=""), eu, units="in", width=3.5, height=3)



eu <- ggplot(df[df$minlod == 4,], aes(y=nMCC_wi_truepc, x=minsegl,group=mod,color = split_time, fill= split_time)) +
    labs(x="min. segment length cutoff", y= "nMCC")+
    geom_point(alpha=0.2, size=0.2, show.legend=F) +
    geom_smooth( method = "lm",span = 0.8,formula = y ~ splines::bs(x, 3), se =TRUE) + 
    scale_colour_manual(name="split time (kya)", values=cols_blue)+
    scale_fill_manual(name="split time (kya)", values=cols_blue)+
    scale_y_continuous(n.breaks=5, limits=c(0.5,1))+
    scale_x_continuous(breaks=segl_break)+
    theme_bw()+ theme(panel.grid.major=element_line(linetype = 2,color="grey", size=0.4),strip.background = element_blank(),
    strip.text=element_text(size=10), axis.text=element_text(size=9), legend.position= "none",axis.text.x=element_text(angle=90, hjust=.5, size=8)) 

ggsave(paste(outpath, "MCC_introgression_eu_4lod.png", sep=""), eu, units="in", width=3.5, height=3)

################# --------------------  FIXED BOTH -------------------- #################

df$split_time = factor(df$split_time, levels=c("2","5","10","15","20","45"))
eu <- ggplot(df[df$minsegl == 50000 & df$minlod == 4,], aes(y=F1_wi_truepc, x=split_time, fill=split_time)) +
    geom_violin(color="black")+
    scale_fill_manual( values=cols_blue)+
    scale_y_continuous(n.breaks=5, expand = c(0,0), limits=c(0,1))+
    labs(x="split time (kya)", y= "F1")+
    stat_summary(aes(y=F1_wi_truepc),fun=median, geom="point", size=2)+
    theme_bw()+ theme( axis.text=element_text(size=9), axis.title=element_text(size=8.5), legend.position="none") 

ggsave(paste(outpath, "F1_introgression_eu_fixed.png", sep=""), eu, units="in", width=3.5, height=3)