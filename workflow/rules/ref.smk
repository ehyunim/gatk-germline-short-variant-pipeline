rule download_vep_plugins:
    output:
        directory("../resources/vep/plugins")
    params:
        release=config["ref"]["release"]
    wrapper:
        "0.66.0/bio/vep/plugins"


rule get_vep_cache:
    output:
        directory("../resources/vep/tmp/cache")
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"]
    wrapper:
        "0.65.0/bio/vep/cache"
