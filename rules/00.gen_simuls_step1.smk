# Step 1: GET VCF from simulations 
# snakemake --snakefile rules/00.gen_simuls_step1.smk --cores 50 --use-conda

configfile: "/maps/projects/racimolab/people/gsd818/arcIntro/introgression-sims/msprime/config_simul.yaml"

import numpy as np
import pandas as pd

SEED= config['model']['seeds']
MODEL = config['model']['name']
TIMES= config['model']['times']


rule all:
    input:
        expand("{model}_{seed}_{time}ky/results/f4_result_ratio.png",model=MODEL, seed=SEED, time=TIMES)

rule msprime:
    # simulate genomes
    output:
        vcf="{model}_{seed}_{time}ky/results/output.vcf.gz",
        trees="{model}_{seed}_{time}ky/results/{model}_{seed}_{time}ky_output_ts.trees",
        node="{model}_{seed}_{time}ky/results/nodes.tsv",
        tracts="{model}_{seed}_{time}ky/results/{model}_{seed}_{time}ky_tracts.tsv.gz"
#    conda: config['envs']
    params:
        genome = config['model']['genome_length'],
        pop_size = config['model']['ne'],
        spath = config['spath']
    log: '{model}_{seed}_{time}ky/logs/msprime.log'
    shell:
        "Rscript {params.spath}/scripts/00.introgression.R "
        "--ne {params.pop_size} "
        "--length {params.genome} "
        "--time {wildcards.time}e3 "
        "--model {wildcards.model}_{wildcards.seed}_{wildcards.time}ky "
        "--seed {wildcards.seed} 2> {log}"

rule f4_ratio:
    # get f4, f4 ratio and Fst
    input:
        trees="{model}_{seed}_{time}ky/results/{model}_{seed}_{time}ky_output_ts.trees",
    output:
        "{model}_{seed}_{time}ky/results/f4_result_ratio.png"
    params:
        spath = config['spath']
    log: '{model}_{seed}_{time}ky/logs/f4.log'
    shell:
        "Rscript {params.spath}/scripts/01.calculate_f4_fst.r "
        "--trees {input.trees} "
        "--time {wildcards.time} "
        "--model {wildcards.model}_{wildcards.seed}_{wildcards.time}ky "
        "-o {output} 2> {log}"
