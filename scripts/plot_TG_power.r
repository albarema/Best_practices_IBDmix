quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(tidyverse))
quiet(library(scales))


output_dir = "modelA_1_45ky/panels_pops/"
if (!dir.exists(output_dir)) dir.create(output_dir)

lods = c(3,4,5,6,10,30,60)
tib= c()
nodes = read_tsv("modelA_1_45ky/results/nodes.tsv", show_col_types = FALSE) %>% select(name, time)

for (i in lods){
    print(i)
    infile = paste("modelA_1_45ky/ibdmix_results/power-allpops-allsamples-allparams-", i,".gz", sep="")
    ttmp = read_tsv(infile,show_col_types = FALSE)
    dtmp_full = merge(ttmp,nodes, by="name", all=TRUE) %>% as_tibble()
    tib = rbind(tib, dtmp_full)
}

# write_tsv(df, "modelA/final_results/ibdmix-new/power-allsamples-allparams.dist_pops_updated.gz")

# split analyses

popsint = c("eurasia", "eurasia.modern", "eurasia.ancient","pop2", "pop2.modern", "pop2.ancient","pop1", "pop1.modern", "pop1.ancient")

df = tib %>% filter(popX %in% popsint) %>% mutate(panel = recode(popX,
                            "eurasia" = "superpop",
                            "eurasia.modern" = "superpop.time",
                            "eurasia.ancient" = "superpop.time",
                            "pop2" = "pop",
                            "pop1" = "pop",
                            "pop1.modern" = "pop.time",
                            "pop2.modern" = "pop.time",
                            "pop1.ancient" = "pop.time",
                            "pop2.ancient" = "pop.time"))

df$n =sapply(df$name, function(x) strsplit(x, "_")[[1]][2]) 
df$pop_anc =sapply(df$name, function(x) strsplit(x, "_")[[1]][1]) 

df$class = "modern"
df$class[df$time != 0] <- "ancient"

lod_breaks <- c(2,4,6,10,30)
labs = c(modern = "present-day individuals", ancient="time-series individuals") 


summs_f <- df %>% group_by(class, panel, minsegl, minlod) %>% 
    summarize(mean = mean(F1_wi_truepc), median = median(F1_wi_truepc), sd = sd(F1_wi_truepc) ) %>% 
    mutate(stats = "Average F1")
best_f <- summs_f %>% group_by(class) %>% filter(abs(median - 1) == min(abs(median - 1))) %>% arrange(desc(minlod), desc(minsegl)) %>% distinct(panel,.keep_all= TRUE)
best_f_mean <- summs_f %>% group_by(class) %>% filter(abs(mean - 1) == min(abs(mean - 1))) %>% arrange(desc(minlod), desc(minsegl))%>% distinct(panel,.keep_all= TRUE)

summs_fs <- df %>% group_by(class,panel, minsegl, minlod,pop_anc) %>% 
    summarize(mean = mean(F1_seg), median = median(F1_seg), sd = sd(F1_seg) ) %>% 
    mutate(stats = "Average F1 segments")
best_fs_mean <- summs_fs %>% group_by(pop_anc) %>% filter(abs(mean - 1) == min(abs(mean - 1))) %>% arrange(desc(minlod), desc(minsegl)) %>% distinct(panel,.keep_all= TRUE)
best_fs <- summs_fs %>% group_by(pop_anc) %>% filter(abs(median - 1) == min(abs(median - 1))) %>% arrange(desc(minlod), desc(minsegl)) %>% distinct(panel,.keep_all= TRUE)


ap <- ggplot(summs_f, aes(x=factor(minlod), factor(minsegl/1000), fill=mean)) +
    geom_tile(color="white")+
    labs(y="minimum segment length [kb]", x="minimum LOD cutoff")+
    geom_tile(data=best_f,colour="black", size =0.3) +
    scale_fill_gradientn(colors=c("navy","white","orange","red"),values=rescale(c(0,0.8,0.9,1)),breaks=c(0,0.8,0.9,1),limits=c(0,1), oob=squish)+
    facet_grid(panel~class,labeller = labeller(class=labs), )+
    theme_bw()+ theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linetype = "dotted"),strip.background = element_blank(),
    strip.text=element_text(size=9), axis.text=element_text(size=9),panel.spacing = unit(0.2, "lines")) +
    guides(fill=guide_colourbar(barheight = 7,ticks.colour= "black", ticks.linewidth = 1,label=T,title = "F1"))

ggsave(paste(output_dir, "F1_introgression_heatmap_panels.png", sep=""), ap, units="in", width=5, height=7)


#  -------------------------------------------------  F1 ------------------------------------------------------
# PLOT DEFAULT PARAMETERS VIOLIN
#  -------------------------------------------------------------------------------------------------------------
p <- ggplot(df[df$minlod == 4 & df$minsegl == 50000,], aes(panel, y=F1_wi_truepc)) +    
    geom_violin()+
    stat_summary(fun=median, geom="point", size=2)+
    geom_hline(yintercept = 0, linetype = 2, colour='grey') +
    facet_wrap(~class, nrow=1, labeller = labeller(pop_anc=mod_labs), scales="free_x")+
    labs(y = "F1", x="Target genomes")+
    scale_y_continuous(n.breaks=3, limits=c(0,1.05),expand=c(0,0))+
    theme_bw()+ theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linetype = "dotted"),
    strip.background = element_blank(), axis.text.x=element_text(angle=45,hjust=1)) 

ggsave(paste(output_dir, "F1_introgression_default_params.png", sep="/"), p, units="in", width=5, height=3)

