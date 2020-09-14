import pandas as pd
import sys
import re
from collections import defaultdict
from snakemake.utils import validate
import glob
import os


###### Config file and sample sheets #####
configfile: "config.yaml"


##### Wildcard constraints #####


##### Variables #####
sample_dir=config["samples"]["path"]

adapter_seq=config["params"]["trimmomatic"]["adapter"]

ref_dir=config["ref"]["path"]
ref_list=os.listdir(ref_dir)
ref_filename=''.join([file for file in ref_list if file.endswith(".fasta")])


##### Helper functions #####
def get_trimmed_reads():
    """Get trimmed reads based on given readingtype""" 
    if str(config["readingtype"]) == "paired":
        return expand('results/trimmed/{sample}.2.unpaired.fastq.gz', sample=sample)
    else:
        return expand('results/trimmed/{sample}.unpaired.fastq.gz', sample=sample)


def get_sample_readingtype():
    """select samples based on given readingtype""" 
    if str(config["readingtype"]) == "paired":
        return {'r1_trim':'results/trimmed/{sample}.1.fastq.gz', 'r2_trim':'results/trimmed/{sample}.2.fastq.gz'}
    else:
        return {'r1_trim':'results/trimmed/{sample}.fastq.gz'}


def get_sample_bams(wildcards):
    return expand("results/recal/{sample}.bam", sample=wildcards.sample)


def get_gvcf_list(wildcards):
    path="/results/called/*"
    pwd=os.getcwd()
    gvcf_dir=pwd+path
    gvcf_list=glob.glob(gvcf_dir)
    """get gvcf list from gvcf files"""
    gvcfs="-V " +" -V ".join([file for file in gvcf_list if file.endswith(".g.vcf.gz")])
    return(gvcfs)


def main():
    get_trimmed_reads()
   # print(get_gvcf_list())

if __name__ == "__main__" :
    main()
