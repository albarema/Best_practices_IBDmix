## ---------------------------------------------------------------
import pandas as pd
## ----------------------------------------------------------------
configfile: 'config_explore.yaml'


MOD = config['model']
LEN=config['length_seg']
NEA=config['neas']
pops_list =  pd.read_table(config['pops_list'])
POPS = pops_list['pop_name'].tolist()


rule all:
    input:
        # cal detection rate and accuracy - calculates the propotion of misestimation according to the whole simulated genome (not only introgressed)
        expand("{model}/final_results/ibdmix/{model}-{len}-{pops}-stats_per_sample_size.txt", model=['modelC', 'modelA', 'modelB'], len='50000', pops='pop1'),
        # cal proportion of of the introgressed genome misestimation 
        expand("{model}/final_results/ibdmix/{model}-{len}-{pops}-stats_per_sample_size_v2.txt", model=['modelC', 'modelA', 'modelB'], len='50000', pops='pop1'),
        #
        expand("compare_Ne_models/all_models_miss_per_sample_size-{len}-{pops}-v2.txt", len='50000', pops=['pop1']),
        expand("compare_Ne_models/all_models_compare_under_line_all_facet_i-{len}-{pops}.png", len='50000', pops=['pop1', 'pop2']),
        expand("{model}/new_figures/misestimation_{pops}_box_{segl}.png",model=['modelA', 'modelB', 'modelC'], pops=['pop1', 'pop2'], segl='50000')



rule plot1_segs:
    input:
        acc="{model}/final_results/ibdmix/acc-allsamples-lod3-complete.txt",
        po="{model}/final_results/ibdmix/power-allsamples-lod3-complete.txt"
    output:
        roc="{model}/new_figures/test_roc_{pops}_points_{segl}.png",
        acc="{model}/new_figures/misestimation_{pops}_box_{segl}.png",
    params:
        scriptsP="/projects/racimolab/people/gsd818/arcIntro/scripts"
    shell:
        "Rscript {params.scriptsP}/02.plot_stats_power_per_segl.r "
        "{input.po} " 
        "{input.acc} "
        "{wildcards.pops} "
        "{wildcards.model}/new_figures/ "
        "{wildcards.segl} "

rule power_compare: 
    input:
        acc="{model}/final_results/ibdmix/acc-allsamples-lod3-complete.txt",
        po="{model}/final_results/ibdmix/power-allsamples-lod3-complete.txt"
    output:
        txt="{model}/final_results/ibdmix/{model}-{len}-{pops}-stats_per_sample_size.txt"
    params:
        scriptsP="/projects/racimolab/people/gsd818/arcIntro/scripts"
    shell:
        "Rscript {params.scriptsP}/02.01.get_summary_performance_models.r"
        " {input.po}"
        " {input.acc}"
        " {wildcards.pops}"
        " {output.txt}"
        " {wildcards.len}"
        " {wildcards.model}"

rule concat: 
    input:
        alls=expand("{model}/final_results/ibdmix/{model}-{len}-{pops}-stats_per_sample_size.txt", model=['modelC', 'modelA', 'modelB'], allow_missing=True),
        one=expand("{model}/final_results/ibdmix/{model}-{len}-{pops}-stats_per_sample_size.txt", model='modelA', allow_missing=True)
    output:
        "compare_Ne_models/all_models_stats_per_sample_size-{len}-{pops}.txt"
    shell:
        "head -n1 {input.one} > {output};"
        " tail -n+2 -q {input.alls} >> {output}"

rule miss_v2: 
    input:
        acc="{model}/final_results/ibdmix/acc-allsamples-lod3-complete.txt",
    output:
        txt="{model}/final_results/ibdmix/{model}-{len}-{pops}-stats_per_sample_size_v2.txt"
    params:
        scriptsP="/projects/racimolab/people/gsd818/arcIntro/scripts"
    shell:
        "Rscript {params.scriptsP}/02.01.b.get_summary_proportion_models.r"
        " {input.acc}"
        " {wildcards.pops}"
        " {output.txt}"
        " {wildcards.len}"
        " {wildcards.model}"

rule concat_v2: 
    input:
        alls=expand("{model}/final_results/ibdmix/{model}-{len}-{pops}-stats_per_sample_size_v2.txt", model=['modelC', 'modelA', 'modelB'], allow_missing=True),
        one=expand("{model}/final_results/ibdmix/{model}-{len}-{pops}-stats_per_sample_size_v2.txt", model='modelA', allow_missing=True)
    output:
        "compare_Ne_models/all_models_miss_per_sample_size-{len}-{pops}-v2.txt"
    shell:
        "head -n1 {input.one} > {output};"
        " tail -n+2 -q {input.alls} >> {output}"


rule plot_compare:
    input:
        f1="compare_Ne_models/all_models_stats_per_sample_size-{len}-{pops}.txt",
        f2="compare_Ne_models/all_models_miss_per_sample_size-{len}-{pops}-v2.txt"
    output:
        p1="compare_Ne_models/all_models_compare_under_line_all_facet_i-{len}-{pops}.png",
        p2="compare_Ne_models/all_models_compare_stats_line_all_facet_i-{len}-{pops}.png",
        p3="compare_Ne_models/all_models_compare_only_miss_line_all_facet_i-{len}-{pops}.png",
        p4="compare_Ne_models/all_models_compare_only_miss_inintro_line_all_facet_i-{len}-{pops}.png"
    params:
        scriptsP="/projects/racimolab/people/gsd818/arcIntro/scripts"
    shell:
        "Rscript {params.scriptsP}/02.02.compare_Ne_models_missestimation.r "
        "{input.f1} "
        "{input.f2} "
        "{output.p1} "
        "{output.p2} "
        "{output.p3} "
        "{output.p4} "

