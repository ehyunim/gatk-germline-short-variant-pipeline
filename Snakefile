include: "workflow/rules/common.snake"

##### Target rules #####
# a target rule to define the desired final output
rule all:
    input:
        #'results/trimmed/BRCA_NA12892.fastq.gz',
        'results/trimmed/BRCA_NA12892.2.unpaired.fastq.gz'


##### Modules #####

include: "workflow/rules/mapping.snake"
