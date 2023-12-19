# Step 1: GET VCF for simulations
# snakemake --snakefile slim_vcf_intro_seed.smk --cores 50 --use-conda

import numpy as np
import pandas as pd

SEED= config['model']['seed']
MODEL = config['model']
TIMES= config['model']['times']


rule all:
    input:
        #expand("{model}_{seed}_{time}ky/results/output.vcf.gz", model=MODEL, seed=SEED, time=TIMES),
        expand("{model}_{seed}_{time}ky/results/{model}_{seed}_{time}ky_tracts.tsv.gz",model=MODEL, seed=SEED, time=TIMES)

rule slim:
    output:
        vcf="{model}_{seed}_{time}ky/results/output.vcf.gz",
        trees="{model}_{seed}_{time}ky/results/{model}_{seed}_{time}ky_output_ts.trees",
        node="{model}_{seed}_{time}ky/results/nodes.tsv",
        tracts="{model}_{seed}_{time}ky/results/{model}_{seed}_{time}ky_tracts.tsv.gz"
    conda: config['envs']
    params:
        genome = config['model']['genome_length'],
        pop_size = config['model']['ne']
    log: '{model}_{seed}/logs/slim.log'
    shell:
        "Rscript scripts/00.introgression.R "
        "--ne {params.pop_size} "
        "--length {params.genome} "
        "--time {wildcards.time}e3 "
        "--model {wildcards.model}_{wildcards.seed}_{wildcards.time}ky "
        "--seed {wildcards.seed} 2> {log}"

rule f4_ratio:
    input:
        trees="{model}_{seed}_{time}ky/results/{model}_{seed}_{time}ky_output_ts.trees",
    output:
        "{model}_{seed}_{time}ky/results/f4_results_ratio.png"
    log: '{model}_{seed}_{time}ky/logs/f4.log'
    shell:
        "Rscript scripts/00.calculate_f4_ratio.r "
        "--trees {input.trees} "
        "--model {wildcards.model}_{wildcards.seed}_{wildcards.time}ky "
        "-o {output} 2> {log}"
