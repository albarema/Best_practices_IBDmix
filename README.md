# Best_practices_IBDmix
Manual to perform good practices when looking for archaic segments using the state-of-the-art tools


# Demogrpahic model

We can use demes [Demes Github](https://popsim-consortium.github.io/demes-spec-docs/main/introduction.html) to visualise the demographic model of interest as follow:
````
demesdraw tubes --log-time modelA.yaml modelA_nlog.svg
````
[![INSERT YOUR GRAPHIC HERE](modelA.png)]()

The name of the model must be specified in the config.yaml and will be used for downstream analyses. 

# Workflow Simulations
### Step 0: Download software and packages 
We run the simulations using the Snakemake workflow management system. Please, download snakemake following the instructions here:
``` https://snakemake.readthedocs.io/en/stable/getting_started/installation.html ```

### Step 1: get your config.yaml ready
Check the toy config.yaml file. You will need to specify:
- name of the model ('modelA')
- path to environmental.yaml file ('path/to/environmetal.yaml')
- your output directory ('path/2/outputs')

### Step 2: run simulations and get tree sequences 
Step 2, uses the Snakefile named "rules/00.gen_simuls_step1.smk" which generates the simulation outputs. The pipeline includes a few jobs to generate the simulations for which we chose the package "slendr" [Slendr Manual](https://www.slendr.net/articles/vignette-05-tree-sequences.html). In our example, we use the simulator "slim" but similarly you can use "msprime". 

You can run the snakefile as follows:
````
snakemake --snakefile rules/00.gen_simuls_step1.smk --cores xx
````
This will generate several files in some directories under your output dir:
- {model}_{seed}/results
- {model}_{seed}/logs
- {model}_{seed}/model


STEP 2

# Evaluation
