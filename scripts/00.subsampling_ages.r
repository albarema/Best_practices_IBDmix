# get multiple replicates for different sample sizes

quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(tidyverse))

## ----------------------------------------------------------------------------
POPS=c("pop1", "pop2")

smalls = seq(2,15)
samsizes=c(smalls,20,30,40,50,80)
reps = seq(1,25)

df = read_tsv("age_groups/modelA_1_45ky/results/nodes.tsv", show_col_types = FALSE) %>% 
                select(name, time) %>% 
                filter(time != 0) %>%
                unique() %>% 
                arrange(desc(time))


check = c()
for (i in POPS){
    dpop = df[startsWith(df$name, i),]
    for (s in smalls){
        r = 0
        for (j in seq(1,length(dpop$name), s)){
            r = r + 1 
            subsam = dpop$name[j:(j+s-1)] %>% as_tibble()
            if (length(dpop$name) < j+s-1 ){next}
            else{
                write_tsv(subsam, paste("reps_samplesize/samples_ages/", i,".s",s,".rep",r,".txt", sep=""), col_names=F)
                check <- rbind(check,  paste(i,".s",s,".rep",r, sep=""))
            }
        }
    }
}

# for (t in 1:length(unique(d$times))){
#     ti = unique(d$times))[t]
#     subsam = d[d$times == ti,]$X1 %>% as_tibble()
#     #write_tsv(subsam, paste("introgression-sims/samples_ages/", i,".s",2,".rep",t,".txt", sep=""), col_names=F)
# }

# d=c()
# for (i in POPS){ d <- rbind(d, i)}
# df = rbind(check, d)

colnames(check) <- c("pop_name")
write_tsv(as_tibble(check), "reps_samplesize/samples_ages/pops_reps_perages_list.tsv")

