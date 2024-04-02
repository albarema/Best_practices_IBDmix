# new figures sampling in ancient inds 
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(tidyverse))
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# setwd("reps_samplesize")
# sampling by ages
podf ="modelA_1_45ky/ibdmix_results/power-ancient_byages-allsamples-allparams-4.gz"
acdf="modelA_1_45ky/ibdmix_results/acc-ancient_byages-allsamples-allparams-4.gz"
# random sampling 
pod = "modelA_1_45ky/ibdmix_results/power-ancient_random-allsamples-allparams-4.gz"
acd = "modelA_1_45ky/ibdmix_results/acc-ancient_random-allsamples-allparams-4.gz"
outpath = "reps_ancient/"
if (!dir.exists(outpath)) dir.create(outpath)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

extract_digit_after_dot_s <- function(input_string) {
  # Use regular expression to match the digit after ".s"
  match_result <- regmatches(input_string, regexpr("\\.s(\\d+)", input_string))

  # Check if a match is found
  if (!is.na(match_result) && length(match_result) >= 1) {
    # Extract the matched digit
    digit_after_dot_s <-  gsub(".s", "", match_result)

    return(as.numeric(digit_after_dot_s))
  } else {
    # Return a default value if no match is found
    return(100)
  }
}
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

smalls <- c(1:15)
d1 = read_tsv(podf,show_col_types = FALSE )%>% filter(minsegl == 50000) %>% mutate(type = "sampling_age_based")
d2 = read_tsv(acdf,show_col_types = FALSE) %>% filter(minsegl == 50000) %>% mutate(type = "sampling_age_based")
ac = read_tsv(acd,show_col_types = FALSE) %>% filter(minsegl == 50000) %>% mutate(type = "randomise_w_replacement")
po = read_tsv(pod,show_col_types = FALSE) %>% filter(minsegl == 50000) %>% mutate(type = "randomise_w_replacement")

dpower = rbind(d1,po)
dacc = rbind(d2,ac)

dpower$samplesize <- lapply(dpower$popX, extract_digit_after_dot_s) %>% as.numeric()
dacc$samplesize <- lapply(dacc$popX, extract_digit_after_dot_s) %>% as.numeric()

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
get_summary <- function(smalltab){
    summ_o = smalltab %>% group_by(type,samplesize,minlod,minsegl) %>%
        summarize(mean = mean(nMCC_wi_truepc), median = median(nMCC_wi_truepc), sd = sd(nMCC_wi_truepc) ) %>% 
        mutate(stats = "nMCC")
    summ_u <- smalltab %>%  group_by(type,samplesize,minlod,minsegl) %>%
        summarize(mean = mean(F1_wi_truepc), median = median(F1_wi_truepc), sd = sd(F1_wi_truepc) )  %>% 
        mutate(stats = "F1")
    summ <- rbind(summ_o, summ_u)
    return(summ)
} 

type_labs = c('sampling_age_based'='sampling grouped by age', 'randomise_w_replacement'= "random sampling with replacement")

ftab_all = get_summary(dpower)
ftab_sf1 = ftab_all[ftab_all$stats == "F1" & ftab_all$samplesize <=15, ]

p1 <- ggplot(ftab_sf1, aes(x=samplesize, y=mean,group=1))+
    geom_line(data=ftab_sf1[ftab_sf1$type != "randomise_w_replacement",],aes(colour=type, group=stats),linetype=1, size=0.5, alpha=.4)+ 
    geom_line(data=ftab_sf1[ftab_sf1$type == "randomise_w_replacement" ,],aes(colour=type, group=stats),linetype=1, size=0.5, alpha=.4)+
    geom_ribbon(data=ftab_sf1[ftab_sf1$type != "randomise_w_replacement",],aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type, group=stats), alpha = .2  ) +
    geom_ribbon(data=ftab_sf1[ftab_sf1$type == "randomise_w_replacement",],aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type, group=stats), alpha = .2 ) +
    geom_line(data=ftab_sf1[ftab_sf1$type != "randomise_w_replacement",],aes(colour=type,y=median, group=stats),linetype="dotted", size=0.5, alpha=.4)+
    geom_line(data=ftab_sf1[ftab_sf1$type == "randomise_w_replacement",],aes(colour=type,y=median, group=stats),linetype="dotted", size=0.5, alpha=.4)+
    scale_y_continuous(expand = c(0,0),limits = c(0, 1))+
    scale_x_continuous(n.breaks = 8)+
    labs(x="sample size", y="", title="F1")+
    scale_colour_manual(values=c("red", "blue"), labels=type_labs)+ 
    scale_fill_manual(values=c("red", "blue"), labels=type_labs)+ 
    theme_bw()+
    #guides(color = guide_legend(nrow = 2, byrow = TRUE))+
    theme(axis.text.y = element_text(size=8),axis.text.x = element_text(size=8,angle = 90, vjust = 0.5, hjust=1),axis.title = element_text(size=9,hjust = 0.5),
        legend.title=element_blank(), legend.position ="bottom", legend.text=element_text(size=8),plot.title = element_text(size=10,hjust = 0.5)) +
    guides(color = guide_legend(nrow = 2))      


ggsave(paste(outpath, "pops.ancient_sampling_types_f1.png", sep=""), p1,units="in", height=3, width=3.5)

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
print("reading file 2...")
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------


get_summary_mis <- function(tib){
    smalltab = tib %>% mutate(overestimation=(as.numeric(detected_introgression) - as.numeric(overlap_introgression)), 
    underestimation = (as.numeric(overlap_introgression) - as.numeric(true_introgression)))
    summ_o <- smalltab %>%
        group_by(type, samplesize) %>% 
        summarize(mean = mean(overestimation), median = median(overestimation), sd = sd(overestimation) ) %>% 
        mutate(stats = "Average overestimation")
    summ_u <- smalltab %>%
        group_by(type, samplesize) %>% 
        summarize(mean = mean(underestimation), median = median(underestimation), sd = sd(underestimation) )  %>% 
        mutate(stats = "Average underestimation")
    summ <- rbind(summ_o, summ_u)
    return(summ)
} 


d2 = get_summary_mis(dacc) %>% filter(samplesize <=15)
type_labs = c('sampling_age_based'='sampling grouped by age', 'randomise_w_replacement'= "random sampling with replacement")

hlines =tibble(stats=unique(d2$stats), vals=c(0,0))


p3 <- ggplot(d2, aes(x=as.numeric(samplesize), y=mean,group=1))+
    geom_hline(data=hlines,aes(yintercept=vals), colour="black")+
    geom_line(data=d2[d2$type != "randomise_w_replacement",],aes(colour=type, group=stats),linetype=1, size=0.5, alpha=.4)+ 
    geom_line(data=d2[d2$type == "randomise_w_replacement" ,],aes(colour=type, group=stats),linetype=1, size=0.5, alpha=.4)+
    geom_ribbon(data=d2[d2$type != "randomise_w_replacement",],aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type, group=stats), alpha = .2  ) +
    geom_ribbon(data=d2[d2$type == "randomise_w_replacement",],aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = type, group=stats), alpha = .2 ) +
    geom_line(data=d2[d2$type != "randomise_w_replacement",],aes(colour=type,y=median, group=stats),linetype="dotted", size=0.5, alpha=.4)+
    geom_line(data=d2[d2$type == "randomise_w_replacement",],aes(colour=type,y=median, group=stats),linetype="dotted", size=0.5, alpha=.4)+
    labs(x="sample size", y = 'pct. of genome',title="Average archaic misestimation")+
    scale_colour_manual(values=c("red", "blue"), labels=type_labs)+ 
    scale_fill_manual(values=c("red", "blue"), labels=type_labs)+ 
    scale_x_continuous(n.breaks = 10)+
    theme_bw()+
    theme(axis.text.y = element_text(size=8),axis.text.x = element_text(size=8,angle = 90, vjust = 0.5, hjust=1),axis.title = element_text(size=9,hjust = 0.5),
        legend.title=element_blank(), legend.position ="bottom", legend.text=element_text(size=8),plot.title = element_text(size=10,hjust = 0.5)) +
    guides(color = guide_legend(nrow = 2))    

ggsave(paste(outpath, "pop2.ancient_sampling_types_miss.png", sep="" ),  p3,units="in", height=3, width=3.5)