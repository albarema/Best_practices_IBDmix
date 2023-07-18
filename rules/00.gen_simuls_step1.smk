# Step 1: GET VCF for simulations
# snakemake --snakefile slim_vcf_intro_seed.smk --cores 50 --use-conda

import numpy as np
import pandas as pd

SEED= config['seed'] 
MODEL = config['models']
TIMES= config['times']


rule all:
    input:
        #expand("{model}_{seed}_{time}ky/results/output.vcf.gz", model=MODEL, seed=SEED, time=TIMES),
        expand("{model}_{seed}_{time}ky/results/{model}_{seed}_{time}ky_tracts.tsv.gz",model=MODEL, seed=SEED, time=TIMES)

rule slim:
    output:
        vcf=protected("{model}_{seed}_{time}ky/results/output.vcf.gz"),
        trees=protected("{model}_{seed}_{time}ky/results/{model}_{seed}_{time}ky_output_ts.trees"),
        node=protected("{model}_{seed}_{time}ky/results/nodes.tsv")
    conda: config['envs']
    log: '{model}_{seed}/logs/slim.log'
    shell:
        "Rscript scripts/00.introgression.R "
        "--ne 5000 "
        "--length 200e6 "
        "--time {wildcards.time}e3 "
        "--model {wildcards.model}_{wildcards.seed}_{wildcards.time}ky "
        "--seed {wildcards.seed} 2> {log}"

rule track:
    input:
        vcf="{model}_{seed}_{time}ky/results/output.vcf.gz",
        trees="{model}_{seed}_{time}ky/results/{model}_{seed}_{time}ky_output_ts.trees",
        node="{model}_{seed}_{time}ky/results/nodes.tsv"
    output:
        "{model}_{seed}_{time}ky/results/{model}_{seed}_{time}ky_tracts.tsv.gz"
    log: '{model}_{seed}_{time}ky/logs/trees_track.log'
    params:
        model="{model}_{seed}_{time}ky/model/",
    shell:
        "python scripts/00.detect_tracts.py "
        "--slendr {params.model} "
        "--trees  {input.trees} "
        "--output {output} 2> {log}"

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
