# Can be done both for a step of 5kya for the groupings or 3kya
quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(tidyverse))
quiet(library("ggpubr"))
quiet(library("glue"))
quiet(library(GenomicRanges))

#--------------------------------------------------------------------------------
read_tracks <- function(mod, tracts, nodes){
  # Read true segments
  ttracks <-read_tsv(tracts,show_col_types = FALSE)
  tnodes = read_tsv(nodes,show_col_types = FALSE) %>% select(-c(pop, node_id))
  # merge individual ID
  mtracts = merge(ttracks, tnodes, by=1) %>% as_tibble()
  true_segments <- mtracts %>%
    dplyr:::rename(start = left, end = right) %>%
    mutate(chrom = "chr1") %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  return(true_segments)
}

get_confusionM_segcounts <- function(query_id, subject_id, hits){
  # Get counts
  # We don't have the TN (if they are not predicted how would they be accounted for? which segmetns? the whole 200kb)
  FP = length(query_id[which(!query_id %in% query_id[queryHits(hits)]),])
  FN = length(subject_id[which(!subject_id %in% subject_id[subjectHits(hits)]),])
  FN_info = subject_id[which(!subject_id %in% subject_id[subjectHits(hits)]),]
  TP=length(overlaps)
  #summary(width(subject_id))
  recall_seg = TP/(TP+FN)
  precision_seg =  TP/(TP+FP)
  F1_seg=2*(precision_seg*recall_seg/(precision_seg+recall_seg))
  
  finalvec <- c(recall_seg, precision_seg, F1_seg)
  return(finalvec)
}
get_overlap <- function(query, subject,hits){
  overlaps<- pintersect(query[queryHits(hits)], subject[subjectHits(hits)])
  overlaps$proportion_subject <- width(overlaps) / width(subject[subjectHits(hits)])
  return(overlaps)
}
expected_length <- function(time, m=0.0225, r=1.28e-8){
  generation_time=29
  time_difference = 55e3 - time
  t = time_difference/generation_time
  len = (abs(1-m)*r*abs(t-1))^-1
  return(len)
}


###### ------------------------------------- 45k -------------------------------------#######
model="modelA_1_45ky"
minlod= 3
minsegl = 50000
# ------------------------------------
df <- c()
tracts=paste(model,"/results/",model , "_tracts.tsv.gz",sep="")
nodes=paste(model, "/results/nodes.tsv",sep="")
true_segments <- read_tracks(model, tracts, nodes)
# write_tsv(test, "modelA_1_45ky_s50kb_l3_true_tracts_exptected.tsv")

ibd_estimated = read_tsv(paste(model, "/ibdmix_temp_with_sites/ibd_summary_combined/nea1_2_1_1_1000.gz", sep=""),show_col_types = FALSE) %>%
    mutate(chrom = "chr1") %>%
    filter(slod >= as.numeric(minlod), length  >= as.numeric(minsegl)) %>%
    filter(pop %in% c("pop1", "pop1.ancient", "pop1.modern", "pop2", "pop2.modern", "pop2.ancient")) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
ibd_estimated$ances <- sapply(strsplit(ibd_estimated$ID, "_"),function(x) x[[1]][1])

pop2test = "pop1"
###### ------------------------------------- 45k -------------------------------------#######

get_hist <- function(true_segments, ibd_estimated, pop2test){
  ibd_popX <- ibd_estimated[ibd_estimated$pop == pop2test,]
  query <- true_segments # query now will be the true 
  subject <- ibd_popX # estimated 
  ftab <- c()
  slim_id <- unique(ibd_popX$ID)
  for (m in 1:length(slim_id)){
    query_id <- query[query$name == slim_id[m]]
    subject_id <- subject[subject$ID == slim_id[m]]
          
    hits <- findOverlaps(query_id, subject_id)
    total_pos_hits = length(query_id)
    no_overlaps <- setdiff(c(1:total_pos_hits), queryHits(hits))
    query_id_no_overlaps <- query_id[no_overlaps] %>% as_tibble()  %>% mutate(hit="FALSE", proportion_subject=0)
    overlaps <- get_overlap(query_id, subject_id, hits) %>% as_tibble()   
    allsegs = rbind(overlaps, query_id_no_overlaps)       
    ftab <- rbind(ftab, allsegs)
  }

return(as_tibble(ftab))
}

tib <- get_hist(true_segments, ibd_estimated, pop2test)
df <- rbind(df, tib)

# get a few more columns to group together similar age groups and segment lengths 
df$len_groups <- cut(df$width, breaks=c(0, seq(20000-1,100000, by=+10000),Inf),
    include.lowest = TRUE,
    labels=c('0-20kb', '20-30kb', '30-40kb','40-50kb','50-60kb','60-70kb','70-80kb','80-90kb','90-100kb', '100-maxkb'))

# df$age_groups <- cut(df$time, breaks=c(0,seq(2000-2,41000, by=+5000),Inf),
#   include.lowest = TRUE,
#   labels=c('present-day', '2-7kya', '7-12kya','12-17kya','17-22kya','22-27kya','27-32kya', '32-37kya', '37-44kya'))
df$age_groups <- cut(df$time, breaks=c(0,seq(2000-2,41000, by=+3000),Inf),
  include.lowest = TRUE,
  labels=c('present-day', '2-5kya', '5-8kya','8-11kya','11-14kya','14-17kya','17-20kya','20-23kya','23-26kya', '26-29kya','29-32kya', '32-35kya', '35-38kya', '38-41kya','41-44kya'))

df$estimated_length = expected_length(df$time)

df$proportion_missing = 1- df$proportion_subject

#write_tsv(p,"expected_length_age_sampled.tsv")
#df = read_tsv("expected_vs_true_lengths_per_group.txt")

summ <- df %>% 
  group_by(age_groups) %>% 
  summarize(mean = mean(width/1e3), median = median(width/1e3), sd = sd(width/1e3))
expected_length = df %>% group_by(age_groups) %>% summarize(mean_E = mean(estimated_length/1e3), 
                                                            median_E=median(estimated_length/1e3))
summ = merge(summ, expected_length, by=1)
summ$age_groups = factor(summ$age_groups,levels=c('present-day', '2-5kya', '5-8kya','8-11kya','11-14kya','14-17kya','17-20kya','20-23kya','23-26kya', '26-29kya','29-32kya', '32-35kya', '35-38kya', '38-41kya','41-44kya'))
legend_labs = c("mean" = "mean(true length)","median"="median(true length)", "mean_E" = "mean(expected length)")

#   plot mean/median of the fragments simulated vs. true 
new_lines = ggplot(summ) +
  geom_line(aes(x =age_groups,y=mean,group=1,colour="mean"),, linetype=2, size=0.75)+
  geom_line(aes(x =age_groups,y=median,group=1,colour="median"), linetype=2, size=0.75)+
  geom_line(aes(x =age_groups,y=mean_E, group=1,colour="mean_E"), linetype=1, size=1)+
  #geom_ribbon(aes(x =age_groups,y=mean,group=1, ymin = mean - sd, ymax = mean + sd), fill = "blue",alpha = .2,show.legend =F) +
  labs(y="length of segments [kb]")+
  #coord_cartesian(ylim = c(0, 500)) + 
  scale_color_manual(name = "", 
          values= c("mean" = "blue", "median" = "dodgerblue", "mean_E" = "black"),
          labels =  c("mean" = "mean(true simulated)","median"="median(true simulated)", "mean_E" = "mean(expected)"))+
  geom_hline(yintercept=50, colour="red", linetype="dotted", size=0.5)+
  scale_y_continuous(breaks=seq(0, 150, by = 50))+
  theme_bw()+
  #guides(color = guide_legend(nrow = 2, byrow = TRUE))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linetype = 2),
  strip.background = element_blank(),strip.text=element_text(size=12),
  legend.position ="right", axis.title.x=element_blank(),
  axis.text.x = element_text(size=11,angle=90,hjust=1),axis.text = element_text(size=12))

ggsave("expected_vs_true_3kya_step.png", new_lines, units="in", width=7, height=7)

#([1−m]r[t−1])−1, where t is the number of generations since a proportion m of the population was replaced by migrants from another population, and r is the recombination rate (in units of Morgans) per bp
# 1/(r*t)

new_lines = ggplot(df,aes(x =age_groups,y=width/1e3)) +
  geom_boxplot(colour="gray")+
  geom_line(data = expected_length,aes(x =age_groups,y=mean_E, group=1), size=1, alpha=.6,colour="black")+
  #geom_ribbon(aes(x =age_groups,y=mean,group=1, ymin = mean - sd, ymax = mean + sd), fill = "blue",alpha = .2,show.legend =F) +
  labs(y="true length simulated tracts [kb]")+
  #coord_cartesian(ylim = c(0, 500)) + 
  geom_hline(yintercept=50, colour="red", linetype="dotted", size=0.5)+
  scale_y_continuous(breaks=seq(0, 1600, by = 250))+
  theme_bw()+
  #guides(color = guide_legend(nrow = 2, byrow = TRUE))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(linetype = 2),
  strip.background = element_blank(),strip.text=element_text(size=12),
  legend.position ="none", axis.title.x=element_blank(),
  axis.text.x = element_text(size=11,angle=90,hjust=1),axis.text = element_text(size=12))
ggsave("boxplot_expected_vs_true_2_5kya_step.png", new_lines, units="in", width=5, height=7)
