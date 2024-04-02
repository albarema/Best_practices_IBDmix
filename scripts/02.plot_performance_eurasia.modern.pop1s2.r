quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(tidyverse))
quiet(library(grid))
quiet(library(scales))

#--------------------------------------------------------------------------------
path2f = "modelA_1_45ky/ibdmix_results/"
outpath = "summary_plots/"
#--------------------------------------------------------------------------------

ps2 = list.files(path=path2f, pattern="power-eurasia_pop1s2")

dps2 <- c()
for (i in 1:length(ps2)){
    d = read_tsv(paste(path2f, ps2[i], sep=""),show_col_types = FALSE)
    dps2 = rbind(dps2, d)
}

dps2$pop_anc =sapply(dps2$name, function(x) strsplit(x, "_")[[1]][1]) 
col_pops = c("#8c564b", "#c49c94")

get_summary <- function(tib){
    summ_o = tib %>% group_by(pop_anc,minlod,minsegl) %>%
        summarize(mean = mean(nMCC_wi_truepc), median = median(nMCC_wi_truepc), sd = sd(nMCC_wi_truepc) ) %>% 
        mutate(stats = "nMCC")
    summ_u <- tib %>%  group_by(pop_anc,minlod,minsegl) %>%
        summarize(mean = mean(F1_wi_truepc), median = median(F1_wi_truepc), sd = sd(F1_wi_truepc) )  %>% 
        mutate(stats = "F1")
    summ <- rbind(summ_o, summ_u)
    return(summ)
}   

ftab = get_summary(dps2)

#summs_m$min_max = abs(summs_m$mean) - abs(df_max$minimise) 
#best_media <- summs_m %>% group_by(pop_anc) %>% filter(abs(min_max) == max(summs_m$min_max)) %>% arrange(desc(minlod), desc(minsegl))%>% distinct(pop_anc,.keep_all= TRUE)

mod_labs = c(pop1 = "pop1 - 100 present-day individuals", pop2="pop2 - 2 present-day individuals") 

# F1
p <- ggplot(ftab[ftab$stats == "F1",], aes(x=factor(minlod), factor(minsegl), fill=mean)) +
    geom_tile(color="white")+
    labs(y="minimum segment length cutoff", x="minimum LOD cutoff")+
    #geom_tile(data=best_media,colour="black", size =0.3) +
    scale_fill_gradientn(colors=c("blue","white","red"),values=rescale(c(0,0.5,1)),breaks=c(0,0.5,1),limits=c(0,1), oob=squish, name="F1")+
    facet_wrap(~pop_anc, ncol=2,labeller = labeller(pop_anc=mod_labs))+
    geom_text(aes(label = round(mean, 2)),size=2) +
    theme_bw()+ theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linetype = "dotted"),strip.background = element_blank(), legend.text=element_text(size = 7), legend.title=element_text(size = 8),
    strip.text=element_text(size=7.5), axis.text=element_text(size=9),panel.spacing = unit(0.2, "lines")) 

ggsave(paste(outpath,"eurasia.modern.pop1s2_F1_heatmap.png", sep=""), p,units="cm", height=8, width=15)


## ----------------------------------------------------------------
# MISESTIMATION PLOTS
# ## ----------------------------------------------------------------

ps = list.files(path=path2f, pattern="acc-eurasia_pop1s2")

dps <- c()
for (i in 1:length(ps)){
    d = read_tsv(paste(path2f, ps[i], sep=""),show_col_types = FALSE)
    dps = rbind(dps, d)
}


dps$underestimation =dps$overlap_introgression - dps$true_introgression 
dps$overestimation =  dps$detected_introgression - dps$overlap_introgression

dps <- dps %>% mutate(underestimation_seg=(overlap_introgression*4e8-true_introgression*4e8)/(true_introgression*4e8), 
    overestimation_seg = (detected_introgression*4e8-overlap_introgression*4e8)/(true_introgression*4e8))

#dps <- dps %>% mutate(underestimation_seg=(1-mean_prop_true), overestimation_seg = (1-mean_prop_detected))

dps$pop_anc =sapply(dps$name, function(x) strsplit(x, "_")[[1]][1]) 
col_pops = c("#8c564b", "#c49c94")

df1 <- dps %>% group_by(pop_anc, minsegl, minlod) %>% 
    summarize(mean = mean(underestimation_seg), median= median(underestimation_seg), max= max(underestimation_seg), min=min(underestimation_seg) ) %>% 
    mutate(stats = "Sequence underestimation", tot_info=paste(pop_anc, minsegl, minlod, sep="-"))

df2 <- dps %>% group_by(pop_anc, minsegl, minlod) %>% 
    summarize(mean = mean(overestimation_seg), median = median(overestimation_seg), max= max(overestimation_seg), min=min(overestimation_seg) ) %>% 
    mutate(stats = "Sequence overestimation", tot_info=paste(pop_anc, minsegl, minlod, sep="-"))

df_tot = rbind(df1, df2)
df_max <- merge(df1, df2[,-c(1:3)], by.x=9, by.y=6) %>% as_tibble()

# MINIMISE OVER AND UNDERESTIMATIOIN 
df_max <- merge(df1, df2[,-c(1:3)], by.x=9, by.y=6) %>% as_tibble() 
df_max$minimise = abs(df_max$mean.x) + df_max$mean.y*10
best_media <- df_max %>% group_by(pop_anc) %>% filter(abs(minimise) == min(minimise)) %>% arrange(desc(minlod), desc(minsegl))%>% distinct(pop_anc,.keep_all= TRUE)

df_max$minimise_median = abs(df_max$median.x) + df_max$median.y*10
best_median <- df_max %>% group_by(pop_anc) %>% filter(abs(minimise_median) == min(df_max$minimise_median)) %>% arrange(desc(minlod), desc(minsegl))%>% distinct(pop_anc,.keep_all= TRUE)

# plot only values of LOD and segl that makes sense

p <- ggplot(df_tot,aes(x=factor(minlod), factor(minsegl))) +
#    geom_tile(data=best_media,aes(x=factor(minlod), factor(minsegl),fill=mean.x),colour="black", size =0.3,position = "dodge") + #
#    geom_tile(data=best_media,aes(x=factor(minlod), factor(minsegl),fill=mean.y),colour="black", size =0.3,position = "dodge") + #
    geom_tile(color="white",data = df_tot, aes(x=factor(minlod), factor(minsegl), fill=mean))+
    labs(y="minimum segment length cutoff", x="minimum LOD cutoff")+
    scale_fill_gradientn(colors=c("red","deeppink","white","green", "darkgreen"),values=rescale(c(-1,-0.5,0,0.5,1)),breaks=c(-1,-0.5,0,0.5,1),limits=c(-1,1), oob=squish, name="Sequence\nmisestimation")+
    facet_grid(stats~pop_anc, labeller = labeller(pop_anc=mod_labs))+
    geom_text(aes(label = round(mean, 2)), size=2) +
    theme_bw()+ theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linetype = "dotted"),strip.background = element_blank(), legend.text=element_text(size = 7), legend.title=element_text(size = 8),
    strip.text=element_text(size=7.5), axis.text=element_text(size=9),panel.spacing = unit(0.2, "lines")) 

ggsave(paste(outpath,"eurasia_moidern_pop1s2_samplesizes_mis_heatmap_best.png", sep=""), p,units="cm",  height=12, width=15)





p <- ggplot(df_small, aes(x=factor(minlod), factor(minsegl), fill=median)) +
    geom_tile(color="white")+
    labs(y="minimum segment length cutoff", x="minimum LOD cutoff")+
    geom_tile(data=best_median,colour="black", size =0.3,position = "dodge") +
    scale_fill_gradientn(colors=c("red","deeppink","white","green", "darkgreen"),values=rescale(c(-1,-0.5,0,0.5,1)),breaks=c(-1,-0.5,0,0.5,1),limits=c(-1,1), oob=squish,name="Sequence\nmisestimation")+
    facet_grid(stats~pop_anc, labeller = labeller(pop_anc=mod_labs))+
    geom_text(aes(label = round(mean, 2)), size=2) +
    theme_bw()+ theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linetype = "dotted"),strip.background = element_blank(),legend.text=element_text(size = 7), legend.title=element_text(size = 8),
    strip.text=element_text(size=7.5), axis.text=element_text(size=9),panel.spacing = unit(0.2, "lines")) 

ggsave(paste(outpath,"eurasia.modern.pop1s2_samplesizes_mis_heatmap_median_best.png", sep=""), p,units="cm", height=12, width=15)
