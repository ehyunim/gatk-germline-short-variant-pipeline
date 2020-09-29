rule select_variant:
    input:
        ref=ref_dir+ref_filename,
        all_vcf="results/genotyped/all.vcf.gz",
    output:
        vartype_vcf= "results/genotyped/vartype/all.{vartype}.vcf.gz"
    params:
        vartype_arg=get_vartype
    log:
        "logs/gatk/SelectVariants/{vartype}.log"
    shell:
        "gatk SelectVariants -R {input.ref} -V {input.all_vcf} "
        "{params.vartype_arg} "
        "-O {output.vartype_vcf} 2> {log}"


rule select_hard:
    input:
        ref=ref_dir+ref_filename,
        vcf="results/genotyped/vartype/all.{vartype}.vcf.gz"
    output:
        vcf="results/filtered/all.{vartype}.hardfiltered.vcf.gz"
    params:
        filters=get_hardfilter
    log:
        "logs/gatk/variantfiltration/{vartype}.log"
    wrapper:
         "0.65.0/bio/gatk/variantfiltration"


rule select_vqsr:
    input:
        vcf="results/genotyped/vartype/all.{vartype}.vcf.gz",
        ref=ref_dir+ref_filename
    output:
        vcf="results/filtered/all.{vartype}.recal.vcf.gz",
        tranches="results/filtered/all.{vartype}.tranches"
    params:
        hapmap=config["filtering"]["vqsr"]["resources"]["hapmap"],
        omni=config["filtering"]["vqsr"]["resources"]["omni"],
        g1k=config["filtering"]["vqsr"]["resources"]["g1k"],
        dbsnp=config["filtering"]["vqsr"]["resources"]["dbsnp"],
        mode=config["params"]["gatk"]["VariantRecalibrator"]["mode"],
        annotation=select_annotation,
        extra=config["params"]["gatk"]["VariantRecalibrator"]["extra"]
    log:
        "logs/gatk/variantrecalibrator/{vartype}.log"
    shell:
        "gatk VariantRecalibrator -R {input.ref} -V {input.vcf} "
        "--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} "
        "--resource:omni,known=false,training=true,truth=false,prior=12.0 {params.omni} "
        "--resource:g1k,known=false,training=true,truth=false,prior=10.0 {params.g1k} "
        "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} "
        "{params.annotation} "
        "{params.extra} "
        "-O {output.vcf} "
        "--tranches-file {output.tranches} 2> {log}"


rule merge_variants:
    input:
        vcfs=expand("results/filtered/all.{vartype}.{filtertype}.vcf.gz",
                    vartype=["snvs", "indels"],
                    filtertype="recal"
                                if config["filteroption"]=="vqsr"
                                else "hardfiltered")

    output:
        vcf="results/filtered/all.vcf.gz"
    params:
        extra=""
    log:
        "logs/gatk/merge-filtered-vcfs.log"
    wrapper:
        "0.66.0/bio/picard/mergevcfs"





#rule select_cnn_m1:
#    input:
#        vcf="results/genotyped/all.{sample}.vcf.gz",
#        ref=ref_dir+ref_filename
#    output:
#        "results/filtered/all.{sample}.annotated.vcf"
#    log:
#        "logs/gatk/CNNscoreVariants/{sample}.log"
#    shell:
#        "gatk CNNScoreVariants -V {input.vcf} -R {input.ref} -O {output} 2> {log}"
