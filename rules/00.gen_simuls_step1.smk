# GET VCF for simulations
#snakemake --snakefile slim_vcf_intro_seed.smk --cores 50 --use-conda

import numpy as np
import pandas as pd
from numpy.random import randint

# conda activate aligns.

SEED=range(1, 12)
MODEL = config['models']
EPATH = onfig['envs']

#For randome seeds: randint(0, 999999999, 10)


rule all:
    input:
        #expand("{model}_{seed}/results/output.vcf.gz", model=MODEL, seed=SEED),
        expand("{model}_{seed}/results/{model}_{seed}_tracts.tsv.gz",model=MODEL, seed=SEED)

rule slim:
    output:
        vcf=protected("{model}_{seed}/results/output.vcf.gz"),
        trees=protected("{model}_{seed}/results/{model}_{seed}_output_ts.trees"),
        node=protected("{model}_{seed}/results/nodes.tsv")
    conda: EPATH
    log: '{model}_{seed}/logs/slim.log'
    shell:
        "Rscript scripts/00.introgression.R "
        "--ne 5000 "
        "--model {wildcards.model}_{wildcards.seed} "
        "--seed {wildcards.seed} 2> {log}"

rule track:
    input:
        vcf="{model}_{seed}/results/output.vcf.gz",
        trees="{model}_{seed}/results/{model}_{seed}_output_ts.trees",
        node="{model}_{seed}/results/nodes.tsv"
    output:
        "{model}_{seed}/results/{model}_{seed}_tracts.tsv.gz"
    log: '{model}_{seed}/logs/trees_track.log'
    params:
        model="{model}_{seed}/model/",
    shell:
        "python scripts/00.detect_tracts.py "
        "--slendr {params.model} "
        "--trees  {input.trees} "
        "--output {output} 2> {log}"

rule f4_ratio:
    input:
        trees="{model}_{seed}/results/{model}_{seed}_output_ts.trees",
    output:
        "{model}_{seed}/results/f4_results_ratio.png"
    log: '{model}_{seed}/logs/f4.log'
    shell:
        "Rscript scripts/00.calculate_f4_ratio.r "
        "--trees {input.trees} "
        "--model {wildcards.model}_{wildcards.seed} "
        "-o {output} 2> {log}"
