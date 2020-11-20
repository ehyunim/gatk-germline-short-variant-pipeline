#Evaluate given VCF file with chm-eval (https://github.com/lh3/CHM-eval) for benchmarking variant calling.
rule chm_eval:
    input:
        kit="/tier4/DSC/GATK/DB/chm-eval-kit",
        vcf="results/funcotator/{sample}.vcf"
    output:
        summary="results/chm-eval/{sample}.summary", # summary statistics
        bed="results/chm-eval/{sample}.err.bed.gz" # bed file with errors
    params:
        extra="",
        build="37"
    log:
        "results/logs/chm-eval/{sample}.log"
    wrapper:
        "0.66.0-5-gfb6b623/bio/benchmark/chm-eval"
