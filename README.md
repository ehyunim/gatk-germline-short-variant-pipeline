## Snakemake pipeline for germline small-variant calls

- This pipeline implements the [GATK best-practices workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) uisng [Snakemake](https://snakemake.readthedocs.io/en/stable/) for germline small-variant calls.

- This pipeline follows the overall concept of the germline pipeline[(here)](https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling) which made by a develpoer of Snakemake but the difference is the way of reading input files, the annotation tool(Funcotator by GATK) and minor things of each rules......  

- GATK Germline short-variant discovery 

## Pipeline Description

This germline pipeline has 5  major steps:

- Mapping -> Calling -> Filteration -> Annotation -> Report



## Usage

#### 1) Set up

Create a copy of isg repository

```
git clone https://gitlab.insilicogen.com/dsc/isgrna.git

```


#### 2) Create conda environment 

Before gettting started, you have to install the Miniconda Python3 distribution [(here)](https://conda.io/en/latest/miniconda.html) 

 
* You can check the description in detail by clicking [Snakemake Installation Guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)  

```
$ conda install -c conda-forge mamba
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
$ conda install gatk4 bwa samtools trimmomatic gatktool fastqc multiqc

```
* Otherwise, Simply run:

```
conda env create -f isgrna/gatk_germline_pipeline/workflow/envs/environment.yml 

```
* Make sure you need to type your Miniconda path in the prefix on the bottom of the .yml file 

#### 3) Configure workflow

Make sure you should configure the workflow according to your needs via editing the config.yaml (you may edit Snakefile to select the modules separately)

#### 4) Sample name constraints

In germline case, you must follow the sample name format below


* In single-end reading:

```
{*SAMPLE NAME}_R1.fastq.gz 

```

*SAMPLE NAME is any charaters or numbers you would like to name (wildcard)

(e.g) brca1_R1.fastq.gz, brca1_contig1_R1.fastq.gz, ......

* if you run multiple samples, it can be 
 brca1_sample1_R1.fastq.gz,  brca1_sample2_R1.fastq.gz, brca1_sample3_R1.fastq.gz, ......

* In pair-end reading:

```
{SAMPLE NAME}_R1.fastq.gz, {SAMPLE NAME}_R2.fastq.gz

```

* if you run multiple samples, it can be 
brca1_sample1_R1.fastq.gz,  brca1_sample1_R2.fastq.gz, brca1_sample2_R1.fastq.gz, brca1_sample2_R2.fastq.gz, ......


#### 5) Execute workflow

* Active conda environment via:

```
conda activate snakemake 

```

In directory where connfig.yaml & Sankefile exist, please type the command below:

```
(snakemake) snakemake --cores $N

```
* you may dry-run for test via: 

```
snakemake -np 

```

## Results

while running the workflow, you can check the progress of creating **result/** directory and files which includes log files and result files for each process.

After successful execution, you can see a HTML report with all results via:

```
snakemake --report report.html

```









