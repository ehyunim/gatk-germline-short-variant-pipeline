rule call_variants:
    input:
        ref=ref_dir+ref_filename,
        bam=get_sample_bams
    output:
        "results/called/{sample}.g.vcf.gz"
    params:
        config["params"]["gatk"]["HaplotypeCaller"]
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    shell:
        "gatk HaplotypeCaller -R {input.ref} -I {input.bam} "
        "-O {output} {params} 2> {log}"


rule combine_calls:
    input:
        ref=ref_dir+ref_filename
    output:
        "results/called/all.g.vcf.gz"
    params:
        gvcfs=get_gvcf_list
    log:
        "logs/gatk/combinegvcfs/all.g.vcf.log"
    shell:
        "gatk CombineGVCFs -R {input.ref} " 
        "{params.gvcfs} -O {output} 2> {log}"


rule genotype_gvcfs:
    input:
        ref=ref_dir+ref_filename,
        gvcf="results/called/all.g.vcf.gz"
    output:
        "results/genotyped/all.vcf.gz"
    params:
        config["params"]["gatk"]["GenotypeGVCFs"]
    log:
        "logs/gatk/genotypegvcfs/all.vcf.log"
    shell:
        "gatk GenotypeGVCFs -R {input.ref} {params} " 
        "-V {input.gvcf} -O {output} 2> {log}"
