quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(tidyverse))

setwd("/projects/racimolab/people/gsd818/arcIntro/introgression-sims/msprime/seeds/")

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
fullset$n =sapply(fullset$name, function(x) strsplit(x, "_")[[1]][2]) 

anc = c(1:172)
mod = c(173:272)
fullset$class <- NA
for (i in 1:dim(fullset)[1]){
    if (fullset$n[i] %in% anc){fullset$class[i] <- "ancient"}
    if (fullset$n[i] %in% mod){fullset$class[i] <- "modern"}
}
fullset$population = paste(fullset$popX, fullset$class, sep=" ")
fullset = arrange(fullset, as.numeric(n))
#write_tsv(fullset, "acc-allsamples-allseeds-lod3-nea1_2-5e4bp.gz")

ymax <- max(fullset$true_introgression)+0.005
p <- ggplot(fullset, aes(popX, true_introgression)) +
    geom_boxplot(outlier.shape = NA, show.legend=F) +
    geom_hline(yintercept = 0.0225, linetype = 2, colour="goldenrod3",size=1.5) +
    geom_jitter(size=1.1,alpha = 0.25, aes(shape=population,fill=population)) +
    ylim(c(0,ymax))+
    scale_fill_manual(name = "Population",
                      labels = c("pop1 ancient", "pop1 modern", "pop2 ancient", "pop2 modern"),
                      values = c("red", "red", "blue", "blue")) +   
    scale_colour_manual(name = "Population",
                      labels = c("pop1 ancient", "pop1 modern", "pop2 ancient", "pop2 modern"),
                      values = c("red", "red", "blue", "blue")) +   
    scale_shape_manual(name = "Population",
                      labels = c("pop1 ancient", "pop1 modern", "pop2 ancient", "pop2 modern"),
                     values = c(21,24,21,24),guide = guide_legend(override.aes = list(size = 2.5,alpha = 1)))+
    theme_bw()+
    theme(axis.title.x = element_blank()) +
    ggtitle("True proportion based on tree sequences") +
    labs(y = "Neanderhal ancestry proportion")

ggsave("true_introgression_treeseq.png", p)

fullset$seed <- sapply(fullset$mod, function(x) strsplit(x, "_")[[1]][2])
fullset$seed = factor(as.numeric(fullset$seed), levels=1:10)

p <- ggplot(fullset, aes(popX, true_introgression)) +
    geom_boxplot(outlier.shape = NA, show.legend=F) +
    geom_hline(yintercept = 0.0225, linetype = 2, colour="goldenrod3",size=1.5) +
    geom_jitter(size=1.1,alpha = 0.25, aes(shape=population,fill=population)) +
    ylim(c(0,ymax))+
    scale_fill_manual(name = "Population",
                      labels = c("pop1 ancient", "pop1 modern", "pop2 ancient", "pop2 modern"),
                      values = c("red", "red", "blue", "blue")) +   
    scale_colour_manual(name = "Population",
                      labels = c("pop1 ancient", "pop1 modern", "pop2 ancient", "pop2 modern"),
                      values = c("red", "red", "blue", "blue")) +   
    scale_shape_manual(name = "Population",
                      labels = c("pop1 ancient", "pop1 modern", "pop2 ancient", "pop2 modern"),
                     values = c(21,24,21,24),guide = guide_legend(override.aes = list(size = 2.5,alpha = 1)))+
    theme_bw()+
    theme(axis.title.x = element_blank(),legend.position = "bottom", legend.direction = "horizontal") +
    ggtitle("True proportion based on tree sequences") +
    facet_wrap(~seed, nrow=2)+
    labs(y = "Neanderhal ancestry proportion")

ggsave("true_introgression_treeseq_faceted.png", p)

################################################################################################################ 
# MODEL A
################################################################################################################

p <- ggplot(fullset[fullset$mod == "modelA_1_45ky",], aes(popX, true_introgression)) +
    geom_boxplot(outlier.shape = NA, show.legend=F) +
    geom_hline(yintercept = 0.0225, linetype = 2, colour="goldenrod3",size=1.5) +
    geom_jitter(size=1.1,alpha = 0.25, aes(shape=population,fill=population)) +
    ylim(c(0,ymax))+
    scale_fill_manual(name = "Population",
                      labels = c("pop1 ancient", "pop1 modern", "pop2 ancient", "pop2 modern"),
                      values = c("red", "red", "blue", "blue")) +   
    scale_colour_manual(name = "Population",
                      labels = c("pop1 ancient", "pop1 modern", "pop2 ancient", "pop2 modern"),
                      values = c("red", "red", "blue", "blue")) +   
    scale_shape_manual(name = "Population",
                      labels = c("pop1 ancient", "pop1 modern", "pop2 ancient", "pop2 modern"),
                     values = c(21,24,21,24),guide = guide_legend(override.aes = list(size = 2.5,alpha = 1)))+
    theme_bw()+
    theme(axis.title.x = element_blank()) +
    ggtitle("True proportion based on tree sequences") +
    labs(y = "Neanderhal ancestry proportion")

ggsave("true_introgression_treeseq_seed1_ylims.png", p)