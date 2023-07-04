# calculate f4 ratio 
quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(slendr))
quiet(library(dplyr))
quiet(library(admixr))
quiet(library(dplyr))
quiet(library(ggplot2))
quiet(library(purrr))
quiet(library(readr))
quiet(library(optparse))


reticulate::use_condaenv("slendr", required = TRUE)

option_list = list(
    make_option(c("-m", "--model"), type="character", default=NULL, help="Model name"),
    make_option(c("-t", "--trees"), type="character", default=5000, help="Tree sequences"),
    make_option(c("-s", "--seed"), type="character", default=314159265, help="seed"),
    make_option(c("-o", "--output"), type="character", default=NULL, help="output name")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

mod <- opt$model
treefile <- opt$trees
outfile <- opt$output

model_dir <- paste("./", mod, "/model",sep="")
output_dir <-  paste("./", mod, "/results",sep="")

if (!dir.exists(model_dir)) dir.create(model_dir)
if (!dir.exists(output_dir)) dir.create(output_dir)

print("reading files")
model <- read(model_dir)
root_ids <- model$splits %>%
    dplyr::filter(pop %in% c("archaic", "AMH")) %>%
    .$pop_id %>%
    { . + 1 }

n_pops <- length(model$populations)
migration_matrix <- matrix(0, nrow = n_pops, ncol = n_pops)
migration_matrix[root_ids[1], root_ids[2]] <- 1
migration_matrix[root_ids[2], root_ids[1]] <- 1

ts <- ts_load(model, file = treefile,
              recapitate = TRUE, Ne = 10000, recombination_rate = 1e-8,
              simplify = TRUE,
              mutate = TRUE, mutation_rate = 1e-8, random_seed = SEED,
              migration_matrix = migration_matrix)

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