rule funcotator_annotate_variants:
    input:
        vcf="results/filtered/all.vcf.gz",
        ref=ref_dir+ref_filename
    output:
        report("results/annotated/funcotator/all.vcf.gz", caption="../report/vcf.rst",  category="Calls")
        #"results/annotated/funcotator/all.vcf.gz"
    params:
        db=config["params"]["gatk"]["Funcotator"]["path"],
        refversion=config["ref"]["version"],
        outputformat= "VCF"
    threads: 4
    log:
        "results/logs/funcotator/annotate.log"
    shell:
        "gatk Funcotator -R {input.ref} -V {input.vcf} -O {output} "
        "--data-sources-path {params.db} "
        "--output-file-format {params.outputformat} "
        "--ref-version {params.refversion} 2> {log}"
