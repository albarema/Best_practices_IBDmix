# conda activate neon
# RUN IBD for simulations using this file. It will use the snakemake file from the common IBDmix dir.Snakefile from this dir
configfile: '/projects/racimolab/people/gsd818/arcIntro/introgression-sims/msprime/pop_size/config_ibdmix.yaml'
MODEL=config['model']

chromosomes = config['chr']

paths = config['paths']
# substituted wildcards
yaml_wildcards = ['output_root', 'sample_name']
for wild in yaml_wildcards:
    for key, filename in paths.items():
        fmt = f'{{{wild}}}'
        if fmt in filename:
            paths[key] = filename.replace(fmt, paths[wild])
paths['exe'] = paths['exe_root'] + '{exe}'
exes = ['generate_gt', 'ibdmix']

if 'sample_file' in paths:
    populations = glob_wildcards(paths['sample_file']).population

if len(populations) == 0:
    paths['ibd_output'] = paths['ibd_output'].replace('{population}', 'ALL')
    paths['ibd_summary'] = paths['ibd_summary'].replace('{population}', 'ALL')
    populations = ['ALL']

include: '/projects/racimolab/people/gsd818/arcIntro/IBDmix/snakefiles/ibdmix.smk'


def all_input(wildcards):
    params = {'chrom': chromosomes, "model":MODEL}
    path = paths['ibd_output']

    if 'summary_lod' in config['IBDmix']:
        path = paths['combined_summary']
        params['LOD'] = config['IBDmix']['summary_lod']
        params['length'] = config['IBDmix']['summary_length']

    if '{population}' in path:
        params['population'] = populations

    return expand(path, **params)

def ibd_mask_bed_input(wildcards):
    return {
        'ibd': paths['ibd_output'],
        'mask': paths['mask_file']
    }

rule all:
    input:
        all_input


rule ibd_mask_bed:
    input:
        unpack(ibd_mask_bed_input)

    output:
        paths['ibd_bed']

    singularity:
        'docker://biocontainers/bedtools:2.25.0'

    shell:
        'awk \''
            'BEGIN{{OFS="\\t"}} '
            'NR==FNR{{'  # build chr,start,end to max and sum masked
                'sum[$1 "." $2 "." $3] += $7; '
                'if(max[$1 "." $2 "." $3] <= $7){{ '
                    'max[$1 "." $2 "." $3] = $7 }} '
                'next; }} '
            'FNR==1{{ print $0, "total_masked", "largest_mask" ; next}} '
            '{{ print $0, sum[$2 "." $3 "." $4], max[$2 "." $3 "." $4]}} \' '
            '<( zcat {input.ibd} | '
                'tail -n +2 | '
                'cut -f2-4 | '
                'sort -k2,3n | '
                'uniq | '
            'bedtools intersect '
                '-a stdin '
                '-b {input.mask} '
                '-wao '  # write all overlap
                '-sorted ) '
            '<( zcat {input.ibd} ) '
            '| gzip > {output} '

def summary_input(wildcards):
    ibd = paths['ibd_output']
    if 'mask_file' in paths and config['IBDmix']['mask_stats']:
        ibd = paths['ibd_bed']
    return {
        'exe': paths['source_root']+'src/summary.sh',
        'ibd': ibd
    }

rule summary:
    input:
        unpack(summary_input)

    output:
        paths['ibd_summary']

    shell:
        'zcat {input.ibd} | '
        '{input.exe} '
            '{wildcards.length} '
            '{wildcards.LOD} '
            '{wildcards.population} '
        '| gzip > {output} '

def combine_input(wildcards):
    return expand(paths['ibd_summary'],
                  chrom=wildcards.chrom,
                  LOD=wildcards.LOD,
                  model=wildcards.model,
                  length=wildcards.length,
                  population=populations)

checkpoint combine:
    input:
        combine_input

    output:
        paths['combined_summary']

    shell:
        "zcat {input} | awk "
            "'NR == 1 || !/^ID\tchr/ {{ print }}' "
        "| gzip -c > {output} "