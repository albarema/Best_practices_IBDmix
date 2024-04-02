#plot_misestimacion_group_ages.r

quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(tidyverse))
#--------------------------------------------------------------------------------
df = c()
minlods = c(1,2,3,4,5,6,10,30,60)
output_dir = "length_age"
if (!dir.exists(output_dirpath)) dir.create(output_dir)

pop2test = "pop1"

for (i in minlods){
    dtmp = read_tsv(paste("modelA_1_45ky/ibdmix_results/acc-pop1-allsamples-allparams-", i,".gz", sep=""),show_col_types = FALSE)
    nodes = read_tsv("modelA_1_45ky/results/nodes.tsv", show_col_types = FALSE) %>% select(name, time)
    dtmp_full = merge(dtmp,nodes, by="name", all=TRUE) %>% as_tibble()
    df = rbind(df, dtmp_full)
}

#--------------------------------------------------------------------------------
expected_length <- function(time, m=0.0225, r=1e-8){
  time_difference=55e3-time
  generation_time=29
  t = time_difference/generation_time
  len = (abs(1-m)*r*abs(t-1))^-1
  return(len)
}

get_summary <- function(tib){
    summ_o = tib %>% group_by(age_groups,minsegl, minlod) %>%
        summarize(mean = mean(overestimation_seg), median = median(overestimation_seg), lower = min(overestimation_seg) , upper= max(overestimation_seg)) %>% 
        mutate(stats = "over")
    summ_u <- tib %>%  group_by(age_groups,minsegl, minlod) %>%
        summarize(mean = mean(underestimation_seg), median = median(underestimation_seg), lower = min(underestimation_seg) , upper= max(underestimation_seg)) %>% 
        mutate(stats = "under")
    summ <- rbind(summ_o, summ_u)
    return(summ)
}   

#--------------------------------------------------------------------------------
dlen =  read_tsv("modelA_1_45ky_s50kb_l3_true_tracts_exptected.tsv", show_col_types = FALSE)
dlen$age_groups <- cut(dlen$time, breaks= c(0,1,15e3,30e3,Inf),
    include.lowest = TRUE,
    labels=c('present-day', '1-15kya', '15-30kya','>30kya'))
dlen$estimated_length = expected_length(dlen$time)

expected_l <- dlen %>% 
  group_by(age_groups) %>% 
  summarize(mean = mean(width/1e3), median = median(width/1e3), sd = sd(width/1e3),mean_E = mean(estimated_length/1e3),median_E=median(estimated_length/1e3))

#--------------------------------------------------------------------------------
fullpops = c(paste(pop2test,"modern",sep="."), paste(pop2test,"ancient",sep="."))

df$n =sapply(df$name, function(x) strsplit(x, "_")[[1]][2]) 

df$underestimation =df$overlap_introgression - df$true_introgression 
df$overestimation =  df$detected_introgression - df$overlap_introgression

df <- df %>% mutate(underestimation_seg=(overlap_introgression*4e8-true_introgression*4e8)/(true_introgression*4e8), 
    overestimation_seg = (detected_introgression*4e8-overlap_introgression*4e8)/(true_introgression*4e8))


fullset = df[df$popX %in% fullpops & df$minlod %in% minlods & df$minsegl >=1000,] %>% distinct()

##### bigger groups 
tab_age <- c()
for (i in unique(fullset$minlod)){
    partset = fullset[fullset$minlod == i,]
    partset$age_groups <- cut(partset$time, breaks= c(0,1,15e3,30e3,Inf),
    include.lowest = TRUE,
    labels=c('present-day', '1-15kya', '15-30kya','>30kya'))
    # partset$len_groups <- cut(partset$minsegl, breaks=c(0, seq(20000-1,100000, by=+10000),Inf),
    #     include.lowest = TRUE,
    #     labels=c('0-20kb', '20-30kb', '30-40kb','40-50kb','50-60kb','60-70kb','70-80kb','80-90kb','90-100kb', '100-maxkb'))
    tab_age <- rbind(tab_age,partset)
}
##### bigger groups 
plot_tab <- tab_age %>% group_by(age_groups, minsegl, minlod) %>% 
  summarize(mean_u = mean(underestimation_seg), mean_o = mean(overestimation_seg), median_u = median(underestimation_seg), median_o = median(overestimation_seg))

# weight the overestimation 
df2mim = plot_tab[plot_tab$median_o <= 0.025,]
df2mim$minimise_median = abs(df2mim$median_u) + df2mim$median_o

best_median <- df2mim %>% group_by(age_groups,minlod) %>% filter(abs(minimise_median) == min(minimise_median)) %>% arrange(desc(minsegl)) %>% distinct(age_groups,minlod,.keep_all= TRUE)
best_median_tot <- df2mim %>% group_by(age_groups) %>% filter(abs(minimise_median) == min(minimise_median)) %>% arrange(desc(minsegl), desc(minlod)) 

seg_breaks = unique(tab_age$minsegl)

tab_age$age_groups = factor(tab_age$age_groups, levels=c('present-day', '1-15kya', '15-30kya','>30kya'))


######################## ----------------------------------------------------------------------
# PLOT USED IN PAPER 
######################## ----------------------------------------------------------------------

seg_breaks = unique(tab_age$minsegl/1e3)
dtib = get_summary(tab_age[!tab_age$minlod %in% c(1,2,30,60),])
p = ggplot(dtib ,aes(x=minsegl/1e3, y=median, group=minsegl/1e3)) +
    geom_linerange(data = dtib[dtib$stats == "under",], aes(ymin = lower, ymax = upper), size=.7)+
    geom_point(data = dtib[dtib$stats == "under",], colour="deeppink",alpha=0.5, size=1.7)+ # (width=0.5/length(unique(df1$region))
    geom_line(data = dtib[dtib$stats == "under",], group=1, colour="deeppink")+ # (width=0.5/length(unique(df1$region))
    scale_colour_hue(h = c(0, 90))+
    scale_x_continuous(breaks=seg_breaks)+
    geom_linerange(data = dtib[dtib$stats == "over",], aes(ymin = lower, ymax = upper), colour="black", size=1)+
    geom_point(data = dtib[dtib$stats == "over",], colour="green",alpha=0.5, size=1.7)+ 
    geom_line(data = dtib[dtib$stats == "over",], group=1, colour="green")+
    #geom_vline(data = expected_l, aes(xintercept= median_E*1e3,, color = "Best value"), linetype="solid", size=0.4)+
    geom_vline(data = expected_l, aes(xintercept= mean_E, color = "Expected length"), linetype="solid", size=0.4)+
    geom_vline(data = expected_l, aes(xintercept= median, color = "median"), linetype="dashed", size=0.4)+
    geom_vline(data = expected_l, aes(xintercept= mean, color = "mean"), linetype="dashed", size=0.4)+
    scale_color_manual(name = "", values = c("Expected length"="brown", median = "blue", mean = "black"))+
    facet_grid(minlod ~ factor(age_groups), scales="free_x")+
    geom_hline(yintercept=0, colour='blue', linetype=3)+
    labs(y = "Introgressed sequence misestimation", x="minimum segment length [Kb]")+
    theme_bw()+ theme(axis.text.x =element_text(angle=90), legend.position = "bottom" ) 
ggsave(paste(output_dir, "/segls_dif_introgression_clusters_",pop2test,"line_lods_v2.png", sep=""), p, units="in", width=8, height=10)

# print("plot 2/2...")
seg_breaks = unique(tab_age$minsegl/1e3)
for (i in unique(tab_age$minlod)){
    d1 = tab_age[tab_age$minlod == i,]
    dtib = get_summary(d1)
    p = ggplot(dtib ,aes(x=minsegl/1e3, y=median, group=minsegl/1e3)) +
    geom_linerange(data = dtib[dtib$stats == "under",], aes(ymin = lower, ymax = upper), size=.7)+
    geom_point(data = dtib[dtib$stats == "under",], colour="deeppink",alpha=0.5, size=1.7)+ # (width=0.5/length(unique(df1$region))
    geom_line(data = dtib[dtib$stats == "under",], group=1, colour="deeppink")+ # (width=0.5/length(unique(df1$region))
    scale_colour_hue(h = c(0, 90))+
    scale_x_continuous(breaks=seg_breaks)+
    geom_linerange(data = dtib[dtib$stats == "over",], aes(ymin = lower, ymax = upper), colour="black", size=1)+
    geom_point(data = dtib[dtib$stats == "over",], colour="green",alpha=0.5, size=1.7)+ 
    geom_line(data = dtib[dtib$stats == "over",], group=1, colour="green")+
    #geom_vline(data = expected_l, aes(xintercept= median_E*1e3,, color = "Best value"), linetype="solid", size=0.4)+
    geom_vline(data = expected_l, aes(xintercept= median_E,, color = "Expected length"), linetype="solid", size=0.4)+
    geom_vline(data = expected_l, aes(xintercept= median, color = "median"), linetype="dashed", size=0.4)+
    geom_vline(data = expected_l, aes(xintercept= mean, color = "mean"), linetype="dashed", size=0.4)+
    scale_color_manual(name = "", values = c("Expected length"="brown", median = "blue", mean = "black"))+
    facet_wrap(factor(age_groups)~., scales="free_x", ncol=2)+
    geom_hline(yintercept=0, colour='blue', linetype=3)+
    labs(y = "Introgressed sequence misestimation", x="minimum segment length [Kb]")+
    theme_bw()+ theme(axis.text.x =element_text(angle=90), legend.position = "bottom" ) 
ggsave(paste(output_dir, "/segls_dif_introgression_clusters_",pop2test,"line_lod",i ,"_v2.png", sep=""), p, units="in", width=6, height=7)}

 

######################## ----------------------------------------------------------------------
#PLOT MINIMUM SEGL AND LOD RECORDED

#  p = ggplot(tab_age ,aes(x=factor(minsegl/1e3), y=truelod, group=minsegl)) +
#     geom_boxplot(outlier.alpha=0.3, outlier.size=0.4)+
#     labs(y = "True minimum LOD", x="minimum segment length [Kb]")+
#     scale_y_continuous(n.breaks=6)+
#     stat_summary(fun=median, colour="red", geom="text", show.legend = FALSE, vjust=-0.7, hjust=1.2,size=2.5, aes( label=round(..y.., digits=0)))+
#     stat_summary(fun= function(z) { quantile(z,0.25) }, colour="red", geom="text", show.legend = FALSE, vjust=1.2, hjust=1.2, size=2.5, aes( label=round(..y.., digits=0)))+
#     stat_summary(fun= function(z) { quantile(z,0.75) }, colour="red", geom="text", show.legend = FALSE, vjust=1.2, size=2.5,hjust=1.2, aes( label=round(..y.., digits=0)))+
#     theme_bw()+ theme(axis.text.x =element_text(angle=90) )  

# ggsave(paste(output_dir, "/segls_distr_lod_",pop2test,".v2.png", sep=""), p, units="in", width=5, height=5)


#  p = ggplot(tab_age ,aes(x=factor(minlod), y=truesegl/1e3, group=minlod)) +
#     geom_boxplot(outlier.alpha=0.3, outlier.size=0.4)+
#     labs(y = "True minimum segment length [Kb]", x="minimum LOD score")+
#     scale_y_continuous(n.breaks=6)+
#     stat_summary(fun=median, colour="red", geom="text", show.legend = FALSE, vjust=-0.7, hjust=0.5,size=2, aes( label=round(..y.., digits=0)))+
#     stat_summary(fun= function(z) { quantile(z,0.25) }, colour="red", geom="text", show.legend = FALSE, vjust=-1.2, hjust=.5, size=2, aes( label=round(..y.., digits=0)))+
#     stat_summary(fun= function(z) { quantile(z,0.75) }, colour="red", geom="text", show.legend = FALSE, vjust=1.2, size=2,hjust=.5, aes( label=round(..y.., digits=0)))+
#     theme_bw()+ theme(axis.text.x =element_text(angle=90) )  

# ggsave(paste(output_dir, "/segls_distr_segls_",pop2test,".png", sep=""), p, units="in", width=6, height=5)

# ###### DISTRIBUTION OF THE  ACTUAL MINIMUM SEGMENT IN EACH INDIVIDUAL WHEN WE SET A MINIMUM 

#  p = ggplot(tab_age ,aes(x=factor(minsegl), y=truesegl, group=minsegl)) +
#     geom_boxplot(outlier.alpha=0.3, outlier.size=0.4)+
#     labs(y = "True minimum segment per sampleID", x="Minimum segment length")+
#     scale_y_continuous(n.breaks=6)+
#     stat_summary(fun=median, colour="red", geom="text", show.legend = FALSE, vjust=-0.7, hjust=0.5,size=2.5, aes( label=round(..y.., digits=0)))+
#     theme_bw()+ theme(axis.text.x =element_text(angle=90) )  

# ggsave(paste(output_dir, "/segls_distr_truesegl_withmin_",pop2test,".png", sep=""), p, units="in", width=5, height=5)


#  p = ggplot(tab_age ,aes(x=factor(minlod), y=truelod, group=minlod)) +
#     geom_boxplot(outlier.alpha=0.3, outlier.size=0.4)+
#     labs(y = "True minimum LOD score per sampleID", x="Minimum LOD score")+
#     scale_y_continuous(n.breaks=6)+
#     stat_summary(fun=median, colour="red", geom="text", show.legend = FALSE, vjust=-0.7, hjust=0.5,size=2, aes( label=round(..y.., digits=0)))+
#     theme_bw()+ theme(axis.text.x =element_text(angle=90) )  

# ggsave(paste(output_dir, "/segls_distr_truelod_withmin_",pop2test,".png", sep=""), p, units="in", width=5, height=5)

#  p = ggplot(tab_age ,aes(x=factor(minsegl), y=truelod, group=minsegl)) +
#     geom_boxplot(outlier.alpha=0.3, outlier.size=0.4)+
#     labs(y = "True minimum LOD", x="minimum segment length")+
#     scale_y_continuous(n.breaks=6)+
#     stat_summary(fun=median, colour="red", geom="text", show.legend = FALSE, vjust=-0.5, size=3, aes( label=round(..y.., digits=0)))+
#     facet_wrap(factor(age_groups)~., scales="free_y", ncol=3)+
#     theme_bw()+ theme(axis.text.x =element_text(angle=90) )  

# ggsave(paste(output_dir, "/segls_distr_lod_clusters_",pop2test,".png", sep=""), p, units="in", width=9, height=7)
