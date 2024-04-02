quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(tidyverse))
quiet(library(grid))
quiet(library(scales))

output_dir = "modelA_1_45ky/panels_pops/"
if (!dir.exists(output_dir)) dir.create(output_dir)

#---------------------------------------------------------
lods = c(3,4,5,6,10,30,60)
df= c()
nodes = read_tsv("modelA_1_45ky/results/nodes.tsv", show_col_types = FALSE) %>% select(name, time)

for (i in lods){
    print(i)
    infile = paste("modelA_1_45ky/ibdmix_results/acc-allpops-allsamples-allparams-", i,".gz", sep="")
    ttmp = read_tsv(infile,show_col_types = FALSE)
    dtmp_full = merge(ttmp,nodes, by="name", all=TRUE) %>% as_tibble()
    df = rbind(df, dtmp_full)
}

#---------------------------------------------------------
lodstest = c(3,4,5,6,10)
segtest = c(15000,30000,40000,50000,60000,90000)

popsint = c("eurasia", "eurasia.modern", "eurasia.ancient","pop2", "pop2.modern", "pop2.ancient","pop1", "pop1.modern", "pop1.ancient")

df2 = df %>% filter(popX %in% popsint) %>% mutate(panel = recode(popX,
                            "eurasia" = "superpop",
                            "eurasia.modern" = "superpop.time",
                            "eurasia.ancient" = "superpop.time",
                            "pop2" = "pop",
                            "pop1" = "pop",
                            "pop1.modern" = "pop.time",
                            "pop2.modern" = "pop.time",
                            "pop1.ancient" = "pop.time",
                            "pop2.ancient" = "pop.time"))

df2$n =sapply(df2$name, function(x) strsplit(x, "_")[[1]][2]) 
df2$pop_anc =sapply(df2$name, function(x) strsplit(x, "_")[[1]][1]) 

df2$class = "modern"
df2$class[df2$time != 0] <- "ancient"

#-------------------------------
df_def <- df2 %>% mutate(underestimation_seg=(overlap_introgression-true_introgression)/(true_introgression), 
    overestimation_seg = (detected_introgression-overlap_introgression)/(true_introgression))

labs = c(modern = "present-day individuals", ancient="time-series individuals") 


##### ------------------

get_summary <- function(df_def){
    df1 <- df_def %>% group_by(class,panel, minsegl, minlod) %>% 
    summarize(mean = mean(underestimation_seg), median= median(underestimation_seg), max= max(underestimation_seg), min=min(underestimation_seg) ) %>% 
    mutate(stats = "Sequence underestimation", tot_info=paste(panel, class,minsegl, minlod, sep="-"))

    df2 <- df_def %>% group_by(class,panel, minsegl, minlod) %>% 
    summarize(mean = mean(overestimation_seg), median = median(overestimation_seg), max= max(overestimation_seg), min=min(overestimation_seg) ) %>% 
    mutate(stats = "Sequence overestimation", tot_info=paste(panel, class,minsegl, minlod, sep="-"))

    df_tot = rbind(df1, df2)
    df_small = df_tot[as.numeric(df_tot$minlod) %in% lodstest & as.numeric(df_tot$minsegl) %in% segtest,]
    df_small$minlod = factor(df_small$minlod, levels=sort(unique(df_small$minlod)))

    df_max <- merge(df1, df2[,-c(1:4)], by.x=10, by.y=6) %>% as_tibble() # %>% filter(minlod %in% lodstest, minsegl %in% segtest)
    df_max$minimise_median = abs(df_max$median.x) + df_max$median.y*10
    best_median <- df_max %>% group_by(class, panel) %>% filter(abs(minimise_median) == min(df_max$minimise_median)) %>% arrange(desc(minlod), desc(minsegl))%>% distinct(panel,.keep_all= TRUE)

    return(list(df_small, best_median))
}

only_one = get_summary(df_def)
only_one_df = only_one[[1]]
only_one_best  = only_one[[2]]

p <- ggplot(only_one_df,aes(x=minlod, factor(minsegl))) +
    #geom_tile(data=best_median,aes(x=factor(minlod), factor(minsegl),fill=mean.x),colour="black", size =0.3,position = "dodge") + #
    #geom_tile(data=best_median,aes(x=factor(minlod), factor(minsegl),fill=mean.y),colour="black", size =0.3,position = "dodge") + #
    geom_tile(color="white",data = only_one_df, aes(x=factor(minlod), factor(minsegl), fill=mean))+
    labs(y="minimum segment length cutoff", x="minimum LOD cutoff")+
    scale_fill_gradientn(colors=c("red","deeppink","white","green", "darkgreen"),values=rescale(c(-1,-0.5,0,0.5,1)),breaks=c(-1,-0.5,0,0.5,1),limits=c(-1,1), oob=squish, name="Sequence\nmisestimation")+
    facet_grid(panel~class+stats, labeller = labeller(pop_anc_complete=mod_labs), scales="free")+
    geom_text(aes(label = round(mean, 2)), size=2) +
    theme_bw()+ theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linetype = "dotted"),strip.background = element_blank(), legend.text=element_text(size = 7), legend.title=element_text(size = 8),
    strip.text=element_text(size=7.5), axis.text=element_text(size=9),panel.spacing = unit(0.2, "lines")) 

    ggsave(paste(output_dir,"dif_introgression_clusters_heatmap_best.png", sep=""), p,units="cm",  height=15, width=21)

p <- ggplot(df_def[df_def$minlod == 4 & df_def$minsegl == 50000,],aes(x=panel, y=underestimation_seg)) +    
    geom_violin()+
    stat_summary(fun=median, geom="point", size=2, color="deeppink")+
    geom_hline(yintercept = 0, linetype = 2, colour='grey') +
    geom_violin(aes(y=overestimation_seg))+
    stat_summary(aes(y=overestimation_seg),fun=median, geom="point", size=2, color="green")+
    labs(y = "prop. archaic segment misestimated", x="Target genomes")+
    scale_y_continuous(n.breaks=8, limits=c(-.5,.5))+
    facet_wrap(~class) +
    theme_bw()+ theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linetype = "dotted"),
    strip.background = element_blank(), axis.text.x=element_text(angle=45,hjust=1)) 

ggsave(paste(output_dir, "dif_introgression_default_params.png", sep="/"), p, units="in", width=5, height=5)

