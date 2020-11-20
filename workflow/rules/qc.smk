rule fastqc:
    input:
        r1=sample_dir+"{sample}_R1.fastq.gz",
        r2=sample_dir+"{sample}_R2.fastq.gz"
    output:
        html_r1="results/qc/fastqc/{sample}_R1_fastqc.html",
        zip_r1="results/qc/fastqc/{sample}_R1_fastqc.zip",
        html_r2="results/qc/fastqc/{sample}_R2_fastqc.zip",
        zip_r2="results/qc/fastqc/{sample}_R2_fastqc.html"
    params: ""
    log:
        "results/logs/fastqc/{sample}.log"
    shell:
        "fastqc -o results/qc/fastqc/ {input.r1} "
        "fastqc -o results/qc/fastqc/ {input.r2} "
        "2> {log}" 

rule samtools_stats:
    input:
        "results/recal/{sample}.bam"
    output:
        "results/qc/samtools-stats/{sample}.txt"
    log:
        "results/logs/samtools-stats/{sample}.log"
    wrapper:
        "0.67.0-18-g93e9558/bio/samtools/stats"


rule collect_hs_metrics:
    input:
        bam="results/recal/{sample}.bam",
        reference=ref_dir+ref_filename,
        bait_intervals=bed_file,
        target_intervals=bed_file
    output:
        "results/qc/hs_metrics/{sample}.txt"
    log:
        "results/logs/hs_metrics/{sample}.log"
    wrapper:
        "0.67.0-22-gab7b31f/bio/picard/collecthsmetrics"


rule multiqc:
    input:
        expand(["results/qc/samtools-stats/{sample}.txt",
                "results/qc/hs_metrics/{sample}.txt",
                "results/qc/fastqc/{sample}_R1_fastqc.zip",
                "results/qc/fastqc/{sample}_R2_fastqc.zip",
                "results/qc/dedup/{sample}.metrics.txt"], sample=SAMPLES)
    output:
        report("results/qc/multiqc.html", caption="../report/multiqc.rst", category="Quality control")
    log:
        "results/logs/multiqc.log"
    wrapper:
        "0.67.0-18-g93e9558/bio/multiqc"
