# Run from /projects/racimolab/people/gsd818/arcIntro/introgression-sims
# Rscript introgression.R -n 5000 -m modelA --time
# Simulating Neanderthal introgression data using *slendr*
# devtools::install_github("bodkan/slendr")
quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(slendr))
quiet(library(dplyr))
quiet(library(admixr))
quiet(library(ggplot2))
quiet(library(purrr))
quiet(library(readr))
quiet(library(optparse))


init_env()

option_list = list(
    make_option(c("-n", "--ne"), type="character", default=5000, help="Ne pop1"),
    make_option(c("-m", "--model"), type="character", default=NULL, help="Model name"),
    make_option(c("-s", "--seed"), type="character", default=314159265, help="seed"),
    make_option(c("-t", "--time"), type="character", default=45e3, help="split time"),
    make_option(c("-l", "--length"), type="character", default=200e6, help="segment length")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$ne)){
    print <- help(opt_parser)
    stop("ne for pop1 name must be supplied.n", call.=FALSE)
}

mod <- opt$model
ne_pop1 <- opt$ne
seglen <- as.numeric(opt$length)
SEED <- opt$seed
split_time <-  as.numeric(opt$time)

set.seed(SEED)

#  ---------------- create dirs ----------------
#mod = "model_A_seed"
#seed for model_A_seed_3687056
sim_m <- paste("./", mod, sep="")
if (!dir.exists(sim_m)) dir.create(sim_m)

model_dir <- paste("./", mod, "/model",sep="")
output_dir <-  paste("./", mod, "/results",sep="")

if (!dir.exists(model_dir)) dir.create(model_dir)
if (!dir.exists(output_dir)) dir.create(output_dir)

# Time difference of 1.070544 hours - 10k
# ---------------- populations ----------------
ancestor <- population("ancestor", time = 600e3, N = 10000, remove = 549e3)
archaic <- population("archaic", time=550e3, N=3000, parent = ancestor)
neaA <- population("neaA", time = 350e3, N = 1000, parent =archaic)
nea2  <- population("nea2",N=1000,time=130e3,parent=neaA, remove=35e3)
nea1A <- population("nea1A",N=1000, time=130e3, parent=neaA)
nea1 <- population("nea1", N= 1000, time=90e3,parent= nea1A,remove=35e3)
intN <- population("intN", N=1000, time=90e3,parent= nea1A,remove=35e3)
AMH <- population("AMH",N=10000,time=550e3, parent = ancestor)
OOA <- population("OOA",N=2200, time=65e3,parent=AMH)
Eurasia <- population("Eurasia",N=5000,time=56e3,parent=OOA)
pop1 <- population("pop1",N=ne_pop1,time=split_time,parent=Eurasia) # testing 5k and 10k
pop2 <- population("pop2",N=5000,time=split_time,parent=Eurasia)
AFR <- population("AFR",N=10000,time=65e3, parent=AMH)

#1.75
# ---------------- gene flow arc to eurasia ----------------
prop1 <- 0.0225
gf <- gene_flow(from = intN, to = Eurasia, rate = prop1, start = 55000, end = 50000)

pops <- list(ancestor, archaic, neaA, nea2, nea1A, nea1, intN, AMH, OOA, Eurasia, pop1, pop2, AFR)

# ---------------- compile model ----------------
model <- compile_model(
  populations = pops, gene_flow = gf,
  generation_time = 29,
  path = model_dir, overwrite = TRUE
)

plot_model(model)

# ---------------- sampling ----------------
nea_samples <- schedule_sampling(model, times = c(70000,40000), list(nea1, 1),list(nea2, 1))

present_samples <- schedule_sampling(model, times = 0, list(AFR, 300), list(pop1, 100), list(pop2, 100))
emh_samples <- schedule_sampling(model, times = seq(split_time, 2000, by = -500), list(pop1, 2), list(pop2, 2))

samples <- rbind(nea_samples, present_samples, emh_samples)


# ---------------- slim ----------------
start_time <- Sys.time()
ts <- msprime(
  model, sequence_length = seglen, recombination_rate = 1e-8,
  samples = samples,
  verbose = TRUE, random_seed = SEED
) %>% ts_mutate(1e-8)

ts_save(ts, file = file.path(output_dir, paste(mod, "_output", sep="")))

end_time <- Sys.time()
end_time - start_time

#print("If eigenstrat files exists - only to read")
#prefix <- file.path(output_dir, "eigenstrat")
#eigen_data <- ts_eigenstrat(ts, prefix = prefix, outgroup = "chimp")


print("Creating vcf files")

ts_vcf(ts, chrom="1", path = file.path(output_dir, "output.vcf.gz"))

print("Creating eigenstrat files")
prefix <- file.path(output_dir, "eigenstrat")
eigen_data <- ts_eigenstrat(ts, prefix = prefix, outgroup = "chimp")

print("nodes...")

node_table <-
  ts_nodes(ts) %>%
  filter(sampled) %>%
  filter(pop %in% c("pop1", "pop2")) %>%
  as_tibble() %>%
  select(name, pop, time, node_id)

write_tsv(node_table, file.path(output_dir,"nodes.tsv"))




