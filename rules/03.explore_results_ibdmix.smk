# calculate and plot power and MCC curves- Run from segment_length_200kb 
## ---------------------------------------------------------------
# Rscript calc_precision_200Mb.R -m {input.model} -l 100000 --nean nea1_2 --popn 1,2 --age {params.age} --outf {output}

import pandas as pd
## ----------------------------------------------------------------

# confi for split_times 
# configfile: '/projects/racimolab/people/gsd818/arcIntro/introgression-sims/msprime/split_times/config_explore_all.yaml'
# for pop_size
# configfile: '/projects/racimolab/people/gsd818/arcIntro/introgression-sims/msprime/pop_size/config_explore.yaml'
# rep_samplesize, age_groups
configfile: '/projects/racimolab/people/gsd818/arcIntro/introgression-sims/msprime/config_explore_all.yaml'

# modelA for age_groups
# MOD = config['model']['seeds'] 
MOD = config['model']['modelA']
#LEN=config['length_seg']
NEA=config['neas']

# uncomment for simuls with pop list files for example reps_samplesize, ne, age_groups
pops_list =  pd.read_table(config['pops_list'])
POPS = pops_list['pop_name'].tolist()

MINSEGS=['1000', '15000', '30000', '40000', '50000', '60000', '90000', '120000']
MINLODS=[1,2,3,4,5,6,10,30,60]

# MINLODS = 3, 
# Uncomment for reps_samplesize 
# MINSEGS = 50000
# uncomment for seeds
# POPS = config['pops_groups']['pops']
# uncomment for split_time
# POPS = ['eurasia.modern']

## ----------------------------------------------------------------

rule all:
    input:
        # for seeds, age_groups
        expand("{model}/ibdmix_results/power-{super}-{pops}-allparams-{lod}.gz", model=MOD, lod=[30,60],super="allpops", pops="allsamples"),
        expand("{model}/ibdmix_results/acc-{super}-{pops}-allparams-{lod}.gz", model=MOD, lod=[30,60],super="allpops", pops="allsamples")
        # cal detection rate and accuracy
        #expand("{model}/ibdmix_results/power-{super}-{pops}-allparams-{lod}.gz", model=MOD, lod=MINLODS,super="eurasia.modern", pops="allsamples"),
        #expand("{model}/ibdmix_results/acc-{super}-{pops}-allparams-{lod}.gz", model=MOD, lod=MINLODS,super="eurasia.modern", pops="allsamples"), 
        
        # allsamples if not pops, -reps.txt

## ----------------------------------------------------------------

rule power:
    input:
        tracts="{model}/results/{model}_tracts.tsv.gz",
        # uncomment for most models with 1_1000: split time, age_groups       
        res="{model}/ibdmix_temp_with_sites/ibd_summary/{nean}_{pops}_1_1_1000.gz"
        # uncomment for most models with less params tested: reps, ne       
        # res="{model}/ibdmix_temp_with_sites/ibd_summary/{nean}_{pops}_1_3_20000.gz"
    output:
        out1 = temp("{model}/ibdmix_results/power-samples-{super}-{nean}-{len}-{lod}-{pops}.txt"),
        out2 = temp("{model}/ibdmix_results/acc-samples-{super}-{nean}-{len}-{lod}-{pops}.txt"),
    params:
        scriptsP="/projects/racimolab/people/gsd818/arcIntro/introgression-sims/msprime/scripts",
    shell:
        "Rscript {params.scriptsP}/02.calc_performance_scores_200Mb.r "
        "--tracts {input.tracts} "
        "--infile {input.res} "
        "--model {wildcards.model} "
        "--segl 1000 "
        "--minsegl {wildcards.len} "
        "--lod 1 "
        "--minlod {wildcards.lod} "
        "--nean {wildcards.nean} "
        "--popn {wildcards.pops} "
        "--outf {output.out1} "
        "--acc {output.out2}"


rule combine_all:
    input:
        alls=expand("{model}/ibdmix_results/power-samples-{super}-{nean}-{len}-{lod}-{pops}.txt",nean = NEA, pops=POPS, len=MINSEGS, allow_missing=True), #  pops=POPS,
        one=expand("{model}/ibdmix_results/power-samples-{super}-{nean}-{len}-{lod}-{pops}.txt",nean=NEA[0], pops=POPS[0],len=50000, allow_missing=True) #  pops=POPS[1]
    output:
        "{model}/ibdmix_results/power-{super}-allsamples-allparams-{lod}.gz"
    shell:
        "head -n1 {input.one} | gzip -c > {output};"
        " tail -n+2 -q {input.alls} | gzip -c >> {output}"

rule combine_acc:
    input:
        alls=expand("{model}/ibdmix_results/acc-samples-{super}-{nean}-{len}-{lod}-{pops}.txt",nean = NEA ,len=MINSEGS, pops=POPS, allow_missing=True),
        one=expand("{model}/ibdmix_results/acc-samples-{super}-{nean}-{len}-{lod}-{pops}.txt",nean=NEA[0], len=50000,pops=POPS[0], allow_missing=True) # eurasia.modern
    output:
        "{model}/ibdmix_results/acc-{super}-allsamples-allparams-{lod}.gz"
    shell:
        "head -n1 {input.one} | gzip -c > {output};"
        " tail -n+2 -q {input.alls} | gzip -c >> {output}"

rule preplot:
    input:
        acc="{model}/ibdmix_results/acc-allsamples-allparams-reps.txt",
        po="{model}/ibdmix_results/power-allsamples-allparams-reps.txt"
    output:
        acc="{model}/ibdmix_results/acc-allsamples-lod3-complete.txt",
        po="{model}/ibdmix_results/power-allsamples-lod3-complete.txt"
    params:
        scriptsP="/projects/racimolab/people/gsd818/arcIntro/scripts"
    shell:
        "Rscript {params.scriptsP}/02.pre_plotting_filters_samplesize.r "
        "{input.po} "
        "{input.acc} "
        "{output.po} "
        "{output.acc} "

rule plot2: # needs summary stats with all lods (to expand in previous rules)
    input:
        acc="{model}/ibdmix_results/acc-allsamples-lod.txt", 
        po="{model}/ibdmix_results/power-allsamples-lod.txt"
    output:
        roc="{model}/figures/test-all_lods.png",
    params:
        scriptsP="/projects/racimolab/people/gsd818/arcIntro/scripts"
    shell:
        "Rscript {params.scriptsP}/power_plotting_lod.r "
        "{input.acc} "
        "{input.po} "
        "{wildcards.pops} "
