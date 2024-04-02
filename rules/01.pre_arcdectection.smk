# THIS FILES GENERATES THE VCF FILES THAT ARE NEEDED FOR IBDMIX AND IN THE RIGHT FORMAT
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
configfile: '/projects/racimolab/people/gsd818/arcIntro/introgression-sims/msprime/config_pre_arcdetection.yaml'
# TO BE RUN FROM THE DIR WHERE SPECIFIC SIMS VERSION IS (where the config.yaml file should be  he right paths)

## --------------------------------------------------------------------------------
import pandas as pd
import csv
from snakemake.io import expand
import itertools

## --------------------------------------------------------------------------------
paths = config['paths']
VERSIOND=config['versiond']
CHROM = "1"
ARCS=[ "nea1_2"] # "nea2_1",

sample_df = pd.read_table("/maps/projects/racimolab/people/gsd818/arcIntro/introgression-sims/list_allpops.txt") 
#sample_df = pd.read_table("/maps/projects/racimolab/people/gsd818/arcIntro/introgression-sims/msprime/seeds/samples/list_pops.tsv")
SUBSETS=pd.unique(sample_df['pop_name'])

## --------------------------------------------------------------------------------
# conda activate neon
wildcard_constraints:
    pops="[^-]+",
    arc1="[^-]+",
    arc2="[^-]+",
    subset="[^-]+",


## ----------------------------------------------------------------

# This is commented out if we are running the pipeline from Snakefile 

rule all:
    input:
        # files needed for IBDMIX, file for archaic, file for the rest of the samples, sample lists
        expand("{vd}/vcf/arc-{arc}.vcf.gz", vd=VERSIOND, arc=ARCS),
        expand("{vd}/vcf/afro-eurasia.vcf.gz",vd=VERSIOND),
        # expand("{vd}/samples/{subset}.txt",vd=VERSIOND, subset=['pop1.modern', 'pop2.modern', 'eurasia', 'all'])


## ----------------------------------------------------------------

rule final_vcf:
    # remove archaic samples from vcf file needed for IBDmix and save file as bgzp compressed (NEEDED for IBDmix)
    input:
        vcf=os.path.join(paths['sim_dir'], "{vd}/results/output.vcf.gz"),
        exsamples="/projects/racimolab/people/gsd818/arcIntro/introgression-sims/msprime/arc.samples.txt",
    output:
        "{vd}/vcf/afro-eurasia.vcf.gz"
    shell:
        "bcftools view -S ^{input.exsamples} {input.vcf} -Oz > {output} ; "
        "bcftools index {output}"

rule get_main_pops_split_times:
    # get samples' list files for each group/cluster that has to be analyse using IBDmix
    input:
        "{vd}/samples-wcluster.txt",
    output:
        "{vd}/samples/{subset}.txt"
    params:
        tools="/projects/racimolab/people/gsd818/arcIntro/scripts",
        outpath="{vd}/samples"
    shell:
        "Rscript {params.tools}/scripts/00.get_sample_lists.r {input} {wildcards.vd} {params.outpath}"

rule arc_vcf:
    # generate archaic vcf as a separate file for IBDmix
    input:
        vcf=os.path.join(paths['sim_dir'], "{vd}/results/output.vcf.gz"),
    output:
        "{vd}/vcf/arc-{arc}.vcf.gz"
    shell:
        "bcftools view -s {wildcards.arc} {input.vcf} -Oz > {output} ; "
        "bcftools index {output}"