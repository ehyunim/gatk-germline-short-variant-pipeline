rule select_hard:
    input:
        ref=ref_dir+ref_filename,
        vcf="results/genotyped/all.{sample}.vcf.gz"
    output:
        vcf="results/filtered/all.{sample}.hardfiltered.vcf.gz"
    params:
        filters=get_hardfilter
    log:
        "logs/gatk/varaintfiltration/{sample}.log"
    wrapper:
         "0.65.0/bio/gatk/variantfiltration"


rule select_vqsr:
    input:
        vcf="results/genotyped/all.{sample}.vcf.gz",
        ref=ref_dir+ref_filename,
    output:
        vcf="results/filtered/all.{sample}.recal.vcf.gz",
        tranches="results/filtered/all.{sample}.tranches"
    params:
        hapmap=config["filtering"]["vqsr"]["resources"]["hapmap"],
        omni=config["filtering"]["vqsr"]["resources"]["omni"],
        g1k=config["filtering"]["vqsr"]["resources"]["g1k"],
        dbsnp=config["filtering"]["vqsr"]["resources"]["dbsnp"],
        mode=config["params"]["gatk"]["VariantRecalibrator"]["mode"],
        annotation=select_annotation,
        extra=config["params"]["gatk"]["VariantRecalibrator"]["extra"]
    log:
        "logs/gatk/variantrecalibrator/{sample}.log"
    shell:
        "gatk VariantRecalibrator -R {input.ref} -V {input.vcf} "
       # "--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} "
        "--resource:omni,known=false,training=true,truth=true,prior=12.0 {params.omni} "
        "--resource:g1k,known=false,training=true,truth=false,prior=10.0 {params.g1k} "
        "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} "
        "{params.annotation} "
        "{params.extra} "
        "-O {output.vcf} "
        "--tranches-file {output.tranches}"

