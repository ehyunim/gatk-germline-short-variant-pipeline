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
"""Get trimmed reads based on given readingtype""" 

def get_trimmed_reads():
    if str(config["readingtype"]) == "paired":
        return expand('results/trimmed/{sample}.2.unpaired.fastq.gz', sample=sample)
    else:
        return expand('results/trimmed/{sample}.unpaired.fastq.gz', sample=sample)


def get_sample_readingtype():
    if str(config["readingtype"]) == "paired":
        return {'r1_trim':'results/trimmed/{sample}.1.fastq.gz', 'r2_trim':'results/trimmed/{sample}.2.fastq.gz'}
    else:
        return {'r1_trim':'results/trimmed/{sample}.fastq.gz'}



def main():
    get_trimmed_reads()


if __name__ == "__main__" :
    main()
