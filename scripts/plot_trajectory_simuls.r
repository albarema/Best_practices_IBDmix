
quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(tidyverse))
quiet(library(grid))
quiet(library(scales))
quiet(library(ggpubr))

#--------------------------------------------------------------------------------
output_dir = "length_age"
minlods=4
minsegls = 5e4

#--------------------------------------------------------------------------------
read_all <- function(mod){
    dataf <- c()
    infile = paste(mod, "/ibdmix_results/acc-allpops-allsamples-allparams-4.gz", sep="")
    d = read_tsv(infile,show_col_types = FALSE)
    nodes = read_tsv( paste(mod, "/results/nodes.tsv", sep=""), show_col_types = FALSE) %>% select(name, time)
    dtmp = merge(d,nodes, by="name", all=TRUE) %>% as_tibble()
    dataf = rbind(dataf, dtmp)
    return(dataf)
}

mods <- paste("modelA_", 1:10, "_45ky", sep = "")

df <- read_all(mods)
fullpops <- c("pop1", "pop2")
fullset = df[df$popX %in% fullpops & df$minsegl==50000 & df$nea == "nea1_2",]

fullset$class = "modern"
fullset$class[fullset$time != 0] <- "ancient"

p <- ggplot(fullset[fullset$class == "ancient",], aes(x=time, y=true_introgression)) + 
  geom_point(size=.6, alpha=0.5) + #  color = clusterAlias, 
  scale_x_continuous(trans = "reverse")+
  ylim(0,0.075)+
  xlab("Years before present") +
  ylab("True Neanderthal proportion") + 
  stat_smooth(method=lm, linewidth=0.55, se=TRUE,linetype = "dashed", fullrange=TRUE,fill="gray", color="blue")+ # fullrange = TRUE,
  stat_cor(method = "pearson" , color="black")+
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), legend.position="bottom",legend.title=element_blank(), strip.background=element_blank())

ggsave(paste(output_dir,"/allreps_true_ancient_modern_trajectory.png", sep=""), p, units="in", width=4, height=3)


p <- ggplot(fullset[fullset$class == "ancient",], aes(x=time, y=detected_introgression)) + # fullset[fullset$class == "ancient",]
  geom_point(size=.6, alpha=0.5) + #  color = clusterAlias, 
  scale_x_continuous(trans = "reverse")+
  ylim(0,0.075)+
  xlab("Years before present") +
  ylab("IBDmix Neanderthal proportion") + 
  stat_smooth(method=lm, linewidth=0.55, se=TRUE,linetype = "dashed", fullrange=TRUE,fill="gray", color="blue")+ # fullrange = TRUE,
  stat_cor(method = "pearson" , color="black")+
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), legend.position="bottom",legend.title=element_blank(), strip.background=element_blank())

ggsave(paste(output_dir,"/allreps_detected_ancient_modern_trajectory.png", sep=""),  p, units="in", width=4, height=3)



########################################## check the pwoer ###################################################
read_power <- function(mod){
    dataf <- c()
    infile = paste(mod, "/ibdmix_results/power-allpops-allsamples-allparams-4.gz", sep="")
    d = read_tsv(infile,show_col_types = FALSE)
    nodes = read_tsv( paste(mod, "/results/nodes.tsv", sep=""), show_col_types = FALSE) %>% select(name, time)
    dtmp = merge(d,nodes, by="name", all=TRUE) %>% as_tibble()
    dataf = rbind(dataf, dtmp)
    return(dataf)
}

mods <- paste("modelA_", 1:10, "_45ky", sep = "")

df <- read_power(mods)
fullpops <- c("pop1", "pop2")
powerset = df[df$popX %in% fullpops & df$minsegl==50000 & df$nea == "nea1_2",]

powerset$class = "modern"
powerset$class[powerset$time != 0] <- "ancient"

p2 <- ggplot(powerset, aes(x=time, y=F1_wi_truepc)) + # fullset[fullset$class == "ancient",]
  geom_point(size=.5, alpha=0.5) + #  color = clusterAlias,
  stat_smooth(method=lm, linewidth=0.55, se=TRUE,linetype = "dashed", fullrange=TRUE,fill="gray", color="black")+ # fullrange = TRUE 
  stat_cor(method = "pearson" , color="black", label.y =.4,size=3)+
  geom_point(data=powerset,aes(y=F1_seg, x=time), colour="red", size=0.5, alpha=.4)+
  stat_smooth(data=powerset, aes(y=F1_seg, x=time),method=lm, linewidth=0.5, se=TRUE,linetype = "dashed", fullrange=TRUE,fill="gray", color="red")+ # fullrange = TRUE,
  stat_cor(data=powerset, mapping=aes(y=F1_seg, x=time),method = "pearson" , color="red", label.y =.2,size=3)+
  scale_x_continuous(trans = "reverse")+
  xlab("Years before present") +
  ylab("F1 score") + 
  ylim(c(0,1))+
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), legend.position="bottom",legend.title=element_blank(), strip.background=element_blank())

ggsave(paste(output_dir,"/F1_ancient_modern_trajectory_reps.png", sep=""), p2, units="in", width=4, height=3)


p2 <- ggplot(powerset[powerset$mod == "modelA_1_45ky",], aes(x=time, y=F1_wi_truepc)) + # fullset[fullset$class == "ancient",]
  geom_point(size=.75, alpha=0.6) + #  color = clusterAlias,
  stat_smooth(method=lm, linewidth=0.55, se=TRUE,linetype = "dashed", fullrange=TRUE,fill="gray", color="black")+ # fullrange = TRUE 
  stat_cor(method = "pearson" , color="black", label.y =.4,size=3)+
  geom_point(data=powerset[powerset$mod == "modelA_1_45ky",],aes(y=F1_seg, x=time), colour="red", size=0.5, alpha=.4)+
  stat_smooth(data=powerset[powerset$mod == "modelA_1_45ky",], aes(y=F1_seg, x=time),method=lm, linewidth=0.55, se=TRUE,linetype = "dashed", fullrange=TRUE,fill="gray", color="red")+ # fullrange = TRUE,
  stat_cor(data=powerset[powerset$mod == "modelA_1_45ky",], mapping=aes(y=F1_seg, x=time),method = "pearson" , color="red", label.y =.2,size=3)+
  scale_x_continuous(trans = "reverse")+
  xlab("Years before present") +
  ylab("F1 score") + 
  ylim(c(0,1))+
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), legend.position="bottom",legend.title=element_blank(), strip.background=element_blank())

ggsave(paste(output_dir,"/F1_ancient_modern_trajectory_1rep.png", sep=""), p2, units="in", width=4, height=3)

