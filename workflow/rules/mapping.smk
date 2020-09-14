rule trim_reads_se:
    input:
        r1=sample_dir+"{sample}_R1.fastq.gz"
    output:
        r1_trim="results/trimmed/{sample}.fastq.gz",
        r1_unpair="results/trimmed/{sample}.unpaired.fastq.gz"
    params:
        pm=config["params"]["trimmomatic"]["trimmer"]
    threads:8
    log:
        "logs/trimmomatic/{sample}.log"
    shell:
        "trimmomatic SE -threads {threads} -phred33 \
        {input.r1} \
        {output.r1_trim} \
        ILLUMINACLIP:{adapter_seq}:2:30:10 {params.pm} \
        2> {log}"


rule trim_reads_pe:
    input:
        r1=sample_dir+"{sample}_R1.fastq.gz",
        r2=sample_dir+"{sample}_R2.fastq.gz"
    output:
        r1_trim="results/trimmed/{sample}.1.fastq.gz",
        r2_trim="results/trimmed/{sample}.2.fastq.gz",
        r1_unpair="results/trimmed/{sample}.1.unpaired.fastq.gz",
        r2_unpair="results/trimmed/{sample}.2.unpaired.fastq.gz"
    params:
        pm=config["params"]["trimmomatic"]["trimmer"]
    threads:8
    log:
        "logs/trimmomatic/{sample}.log"
    shell:
        "trimmomatic PE -threads {threads} -phred33 \
        {input.r1} {input.r2} \
        {output.r1_trim} {output.r1_unpair} \
        {output.r2_trim} {output.r2_unpair} \
        ILLUMINACLIP:{adapter_seq}:2:30:10 {params.pm} \
        2> {log}"

rule map_reads:
    input:
        ref=ref_dir+ref_filename,
        **get_sample_readingtype()

    output:
        "results/mapped/{sample}.bam"
    threads:8
    log:
        "logs/bwa_mem/{sample}.log"
    shell:
        "(bwa mem -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"


rule add_readgroups:
    input:
        "results/mapped/{sample}.bam"
    output:
        "results/mapped/{sample}_rg.bam"
    params:
        rg=config["params"]["bwa"]["readgroup"]
    log:
        "logs/bwa_mem/rg/{sample}.log"
    shell:
        "gatk AddOrReplaceReadGroups -I {input} -O {output} {params.rg} 2> {log}"


rule sort_bam:
    input:
        "results/mapped/{sample}_rg.bam"
    output:
        "results/sorted/{sample}.sorted.bam"
    log:
        "logs/samtools/sort/{sample}.log"
    shell:
        "samtools sort results/mapped/{wildcards.sample}_rg.bam"
        " -o {output} 2> {log}"


rule mark_duplicates:
    input:
        "results/sorted/{sample}.sorted.bam"
    output:
        bam="results/dedup/{sample}.dedup.bam",
        metrics="results/dedup/{sample}.metrics.txt"
    params:
        config["params"]["gatk"]["MarkDuplicates"]
    log:
        "logs/picard/dedup/{sample}.log"
    shell:
        "gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics} 2> {log}"

rule sortsam:
    input:
        "results/dedup/{sample}.dedup.bam"
    output:
        "results/dedup/{sample}.sortsam.bam"
    log:
        "logs/gatk/sortsam/{sample}.log"
    shell:
       "gatk SortSam -I {input} -O {output} -SO coordinate 2> {log}"


rule recalibrate_base_qualities:
    input:
        bam="results/dedup/{sample}.sortsam.bam",
        ref=ref_dir+ref_filename,
        known=config["knowndb"]["dbSNP"]
    output:
        recal_table="results/recal/{sample}.recal.table"
    params:
        config["params"]["gatk"]["BaseRecalibrator"]
    log:
        "logs/gatk/bqsr/{sample}.log"
    shell:
        "gatk BaseRecalibrator -I {input.bam} -R {input.ref} --known-sites {input.known} "
        "{params} -O {output.recal_table}"


rule apply_bqsr:
    apply_bqsr
        bam="results/dedup/{sample}.sortsam.bam",
        ref=ref_dir+ref_filename,
        recal_table="results/recal/{sample}.recal.table"
    output:
        final_bam=protected("results/recal/{sample}.bam")
    log:
        "logs/gatk/bqsr/applybqsr/{sample}.log"
    shell:
        "gatk ApplyBQSR -I {input.bam} -R {input.ref} --bqsr-recal-file {input.recal_table} "
        "-O {output.final_bam} 2> {log}"


rule index_bam:
    input:
        "results/recal/{sample}.bam"
    output:
        "results/recal/{sample}.bam.bai"
    shell:
        "samtools index {input}"


















