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
sample_list=os.listdir(sample_dir)
SAMPLES=[sample.split("_R")[0] for sample in sample_list if sample.endswith(".fastq.gz")]
#SAMPLES=["BRCA_NA12892"]
adapter_seq=config["params"]["trimmomatic"]["adapter"]
ref_dir=config["ref"]["path"]
ref_list=os.listdir(ref_dir)
ref_filename=''.join([file for file in ref_list if file.endswith(".fasta")])
bed_file=config["params"]["gatk"]["CollectHsMetrics"]["bedfile"]


##### Helper functions #####
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


def get_vartype(wildcards):
    return "--select-type-to-include {}".format(
            "SNP" if wildcards.vartype == "snvs" else "INDEL")


def get_vqsr_mode(wildcards):
    return  "--mode {}".format(
            "SNP" if wildcards.vartype == "snvs" else "INDEL")


def get_hardfilter(wildcards):
    return {"snv-hard-filter":
            config["filtering"]["hard"][wildcards.vartype]}#"snvs" or "indels"


def select_annotation(wildcards):
    annotation_list=config["params"]["gatk"]["VariantRecalibrator"]["annotation"]
    annotations="-an " +" -an ".join(annotation_list)
    return annotations

