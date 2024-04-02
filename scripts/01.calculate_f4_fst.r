# calculate f4 ratio 
quiet <- function(x) { suppressMessages(suppressWarnings(x)) }

quiet(library(slendr))
quiet(library(dplyr))
quiet(library(admixr))
quiet(library(ggplot2))
quiet(library(purrr))
quiet(library(readr))
quiet(library(optparse))

# init env to make sure we load the right version of packages 
init_env()

option_list = list(
    make_option(c("-m", "--model"), type="character", default=NULL, help="Model name"),
    make_option(c("-t", "--trees"), type="character", default=5000, help="Tree sequences"),
    make_option(c("-s", "--seed"), type="character", default=314159265, help="seed"),
    make_option(c("-i", "--time"), type="character", default=45, help="split"),
    make_option(c("-e", "--eigen"), type="character", default=NULL, help="eigen data"),
    make_option(c("-o", "--output"), type="character", default=NULL, help="output name")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

mod <- opt$model
treefile <- opt$trees
outfile <- opt$output
seedN = opt$seed
split_time <- as.numeric(opt$time)*1e3

model_dir <- paste("./", mod, "/model",sep="")
output_dir <-  paste("./", mod, "/results",sep="")

if (!dir.exists(model_dir)) dir.create(model_dir)
if (!dir.exists(output_dir)) dir.create(output_dir)

print("reading files")
eigen_data = eigenstrat(paste(output_dir,"/eigenstrat", sep=""))

model <- read_model(model_dir)
ts <- ts_load(model = model, file = treefile) %>%
      ts_mutate(mutation_rate = 1e-8, random_seed = seedN)

samples <- ts_samples(ts)

print("plotting")
# f4 statistics 
# nea1_2 

f4_result <-
  filter(samples, time == 0, pop %in% c("AFR", "pop1", "pop2")) %>%
  group_by(pop) %>%
  sample_n(25) %>%
  pull(name) %>%
  f4(W = "AFR_1", X = ., Y = "nea1_2", Z = "chimp", data = eigen_data)


inner_join(f4_result, samples, by = c("X" = "name")) %>%
  ggplot(aes(pop, f4, color = pop)) +
    geom_pointrange(aes(ymin = f4 - 2 * stderr, ymax = f4 + 2 * stderr, group = X),
                    position = position_dodge(width = 1)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    geom_hline(yintercept = 0, linetype = 2) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
    ggtitle("f4(AFR1, <African or European>; Neanderthal, Chimp)",
            "Statistical test for the presence of Neanderthal ancestry")

 ggsave(paste(output_dir, "f4_result_stats.png", sep="/"))


 ### f4 ratio

 f4ratio_result <-
  filter(samples, pop %in% c("pop1", "pop2")) %>%
  pull(name) %>%
  f4ratio(X = ., A = "nea1_2", B = "nea2_1", C = "AFR_1", O = "chimp", data = eigen_data)
#B=Altai=nea2_1 and A=Vindijia=nea1_2

inner_join(f4ratio_result, samples, by = c("X" = "name")) %>%
  filter(pop %in% c("pop1", "pop2")) %>%
  ggplot(aes(time, alpha)) +
    geom_pointrange(aes(ymin = alpha - 2 * stderr, ymax = alpha + 2 * stderr)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_smooth(method = "lm") +
    coord_cartesian(ylim = c(0, 0.1)) +
    xlim(split_time, 0) +
    facet_wrap(~ pop) +
    ggtitle("f4-ratio estimate of Neanderthal ancestry in Europeans over time") +
    labs(x = "time [years before present]", y = "Neanderhal ancestry proportion")

 ggsave(paste(output_dir, "f4_result_ratio_time.png", sep="/"))

inner_join(f4ratio_result, samples, by = c("X" = "name")) %>%
filter(time == 0) %>%
  ggplot(aes(pop, alpha, color = pop)) +
    geom_pointrange(aes(ymin = alpha - 2 * stderr, ymax = alpha + 2 * stderr, group = X),
                    position = position_dodge(width = 1)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    geom_hline(yintercept = 0, linetype = 2) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
    ggtitle("f4-ratio estimate of Neanderthal ancestry in between Europe and Asia") +
    labs(y = "Neanderhal ancestry proportion")

 ggsave(paste(output_dir, "f4_result_ratio.png", sep="/"))

print("calculating Fst")

# sets from slendr package
sample_sets <- ts_samples(ts) %>%
  split(., .$pop) %>%
  lapply(function(pop) pop$name)

# create new sets to check FST based on ages
emh_samples <- c()
modern <- c()
modern$time <- 0
modern$class <- "modern"
emh_samples$class <- "ancient"

if (split_time == 45e3){
  time_age=seq(split_time + (-500), 2000, by = (-500))}
if (split_time < 45e3){
  time_age=seq(split_time + (-500), 2000, by = (-500))}

emh_samples$time <- sort(rep(time_age,2), decreasing=TRUE)
emh_samples$n <- c(1:length(emh_samples$time))
emh_samples$pop1 <- paste("pop1", emh_samples$n, sep="_")
emh_samples$pop2 <- paste("pop2", emh_samples$n, sep="_")
modern$n <- c((length(emh_samples$time)+1): (length(emh_samples$time) + 100))
modern$pop1 <-  paste("pop1", modern$n, sep="_")
modern$pop2 <-  paste("pop2", modern$n, sep="_")
tot = rbind(as_tibble(emh_samples), as_tibble(modern))

if (split_time  == 45e3){
  tot$age_groups <- cut(tot$time, breaks=c(0,seq(2000-2,41000, by=+3000),Inf),
  include.lowest = TRUE,
  labels=c('present-day', '2-5kya', '5-8kya','8-11kya','11-14kya','14-17kya','17-20kya','20-23kya','23-26kya', '26-29kya','29-32kya', '32-35kya', '35-38kya', '38-41kya','41-44kya'))

}
sample_set <- c()
sample_set$pop1.ancient <- tot[tot$class == "ancient",]$pop1
sample_set$pop2.ancient <- tot[tot$class == "ancient",]$pop2
sample_set$pop2.modern <- tot[tot$class == "modern",]$pop2
sample_set$pop1.modern <- tot[tot$class == "modern",]$pop1

all_sets = c(sample_set,sample_sets)
fst <- ts_fst(ts, sample_sets = all_sets) %>% as.matrix()

write_tsv(as_tibble(fst),paste(output_dir, "fst_pairwise.tsv", sep="/"))