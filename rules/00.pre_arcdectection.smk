
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
configfile: 'config_pre_arcdetection.yaml'
# TO BE RUN FROM THE DIR WHERE SPECIFIC SIMS VERSION IS (where the config.yaml file should be witht he right paths)

## --------------------------------------------------------------------------------
import pandas as pd
import csv
from snakemake.io import expand
import itertools

## --------------------------------------------------------------------------------
paths = config['paths']
VERSIOND=config['versiond']
CHROM = "1"
ARCS=["nea2_1", "nea1_2"]

sample_df = pd.read_table("/maps/projects/racimolab/people/gsd818/arcIntro/introgression-sims/list_allpops.txt") 
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
        # files needed for IBDMIX
        expand("{vd}/vcf/arc-{arc}.vcf.gz", vd=VERSIOND, arc=ARCS),
        expand("{vd}/vcf/afro-EUR.vcf.gz",vd=VERSIOND, arc=ARCS),
        expand("{vd}/samples/{subset}.txt",vd=VERSIOND, subset=['pop1.modern', 'pop2.modern', 'eurasia', 'all'])
        # files needed for SPRIME
        #expand("{vd}/vcf/afro-sub.{subset}.vcf.gz", vd=VERSIOND, subset=SUBSETS), # changed sub-{subset} for sub.{subset}
        #expand("{vd}/vcf/nean-mask-{pops}.bed",vd=VERSIOND, pops='EUR')

## ----------------------------------------------------------------

rule get_pops:
    input:
        vcf=os.path.join(paths['sim_dir'], "{vd}/results/output.vcf.gz"),
    output:
        sam="{vd}/samples.txt",
        arc="{vd}/arc.samples.txt",
        all="{vd}/samples-wcluster.txt",
        afr="{vd}/samples-outgroup-AFR.txt",
        merged="{vd}/samples-wcluster-EUR-merged.txt"
    params:
        tools="/projects/racimolab/people/gsd818/arcIntro/scripts"
    shell:
        "bcftools query -l {input} > {output.sam} ; "
        "Rscript {params.tools}/pops_sims.r {output.sam} {output.arc} {output.afr} {output.all} {output.merged} "

rule final_vcf:
    input:
        vcf=os.path.join(paths['sim_dir'], "{vd}/results/output.vcf.gz"),
        exsamples="{vd}/arc.samples.txt",
    output:
        "{vd}/vcf/afro-EUR.vcf.gz"
    shell:
        "bcftools view -S ^{input.exsamples} {input.vcf} -Oz > {output} ; "
        "bcftools index {output}"

rule get_main_pops_split_times:
    input:
        "{vd}/samples-wcluster.txt",
    output:
        "{vd}/samples/{subset}.txt"
    params:
        tools="/projects/racimolab/people/gsd818/arcIntro/scripts",
        outpath="{vd}/samples"
    shell:
        "Rscript {params.tools}/dist_pops_sims_files.r {input} {wildcards.vd} {params.outpath}"

rule arc_vcf:
    input:
        vcf=os.path.join(paths['sim_dir'], "{vd}/results/output.vcf.gz"),
    output:
        "{vd}/vcf/arc-{arc}.vcf.gz"
    shell:
        "bcftools view -s {wildcards.arc} {input.vcf} -Oz > {output} ; "
        "bcftools index {output}"

rule subsets_vcf:
    input:
        vcf=os.path.join(paths['sim_dir'], "{vd}/results/output.vcf.gz"),
        exsamples="/projects/racimolab/people/gsd818/arcIntro/introgression-sims/samples/{subset}.txt", # the sampling does not change
    output:
        "{vd}/vcf/afro-sub.{subset}.vcf.gz"
    shell:
        "bcftools view -S {input.exsamples} {input.vcf} -Oz > {output} ; "
        "bcftools index {output}"
rule get_map: 
    """
    Create genetic map, all distances set to 0
    """
    input:
        "{vd}/vcf/afro-{pops}.vcf.gz"
    output:
        temp("{vd}/vcf/afro-{pops}-temps.map")
    shell:
        """ bcftools query -f '%CHROM\t%POS\n' {input} | awk -vOFS="\t" '{{print 1, ".", 0, $2}}' > {output} """
rule rmv_dp: 
    input:
        "{vd}/vcf/afro-{pops}-temps.map"
    output:
        temp("{vd}/vcf/afro-{pops}-temps2.map")
    shell:
        """
        LC_ALL=C sort -t $'\t' -k1,1 -k4,4n {input} |\
        awk -F '\t' '/^#/ {{print;prev="";next;}} {{key=sprintf("%s\t%s\t%s",$1,$2,$4);if(key==prev) next;print;prev=key;}}' > {output}
        """

rule inter_map:
    input:
        "{vd}/vcf/afro-{pops}-temps2.map"
    output:
        "{vd}/vcf/afro-{pops}-interpolated.map"
    params:
        tools="/projects/racimolab/people/gsd818/arcIntro/scripts",
        len=2e8
    shell:
        "python {params.tools}/map-interpolation.py {input} {params.len} {output}"

rule bed_file:
    input:
        "{vd}/vcf/afro-{pops}-interpolated.map"
    output:
        "{vd}/vcf/nean-mask-{pops}.bed"
    params:
        tools="/projects/racimolab/people/gsd818/arcIntro/scripts",
    shell:
        "Rscript {params.tools}/get_neanbed.R {input} {output}"