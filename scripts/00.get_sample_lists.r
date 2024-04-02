# get file in samples/{population}.txt
quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(tidyverse))

args = commandArgs(trailingOnly=TRUE)

sampl = args[1]
sample_name = args[2]
outpath = args[3]

# d = read_tsv("modelA_10k/samples-wcluster.txt",col_names=F, show_col_types = FALSE)
d = read_tsv(sampl, col_names=F, show_col_types = FALSE)
split_time = as.numeric(str_remove(strsplit(sample_name, "_")[[1]][2], "k"))*1e3

by_sampling=-500
time_age=seq(split_time + by_sampling, 2000, by = by_sampling)
time_anc_size =length(time_age)*2
anc = c(1:time_anc_size)

pops <- c('pop1', 'pop2')

all = d[!startsWith(d$X2, "nea"),1] %>% as_tibble()
write_tsv(all,paste(outpath, "all.txt", sep="/"), col_names=F)

eur = d[!startsWith(d$X2, "nea") & !startsWith(d$X2, "AFR"),1] %>% as_tibble()
write_tsv(eur,paste(outpath, "eurasia.txt", sep="/"), col_names=F)

for (i in pops){
    newd = d[d$X2 == i,]$X1 %>% as_tibble()
    write_tsv(newd, paste(outpath,  "/",i, ".txt", sep=""), col_names=F)
    newd$class = sapply(newd$value, function(x) strsplit(x, "_")[[1]][2])
    moder = c(time_anc_size+1:dim(newd)[1])
    anc_d = newd[newd$class %in% anc,1] %>% as_tibble()
    mod_d = newd[newd$class %in% moder,1] %>% as_tibble()
    write_tsv(anc_d, paste(outpath, "/", i, ".ancient.txt", sep=""), col_names=F)
    write_tsv(mod_d, paste(outpath, "/",i, ".modern.txt", sep=""), col_names=F)

}

