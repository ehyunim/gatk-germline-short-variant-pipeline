#rule vep_annotate_variants:
#    input:
#       calls="results/filtered/all.vcf.gz",
#	cache="../resources/vep/cache",
#	plugins="../resources/vep/plugins"
#    output:
#        calls=report("results/annotated/vep/all.vcf.gz", caption="./report/vcf.rst", category="Calls"),
#        stats=report("results/stats/all.stats.html", caption="./report/stats.rst", category="Calls")
#    params:
#        plugins=["LoFtool","CADD"],
#	#plugins=config["params"]["vep"]["plugins"],
#	extra=config["params"]["vep"]["extra"]
#    log:
#        "logs/vep/annotate.log"
#    threads: 4
#    wrapper:
#        "0.66.0/bio/vep/annotate"


rule funcotator_annotate_variants:
    input:
        vcf="results/filtered/all.vcf.gz",
        ref=ref_dir+ref_filename
    output:
        "results/annotated/funcotator/all.vcf.gz"
    params:
        db=config["params"]["gatk"]["Funcotator"]["path"],
        refversion=config["ref"]["version"],
        outputformat= "VCF"
    threads: 4
    log:
        "logs/funcotator/annotate.log"
    shell:
        "gatk Funcotator -R {input.ref} -V {input.vcf} -O {output} "
        "--data-sources-path {params.db} "
        "--output-file-format {params.outputformat} "
        "--ref-version {params.refversion} 2> {log}"
