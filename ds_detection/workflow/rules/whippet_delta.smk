rule whippet_delta:
    input:
        unpack(delta_input)
    output:
        config["BASE_PATH"]+"/results/delta_psi/{experiment}.diff.gz"
    params:
        grp1 = lambda wildcards: grp1(wildcards),
        grp2 = lambda wildcards: grp2(wildcards),
        out = config["BASE_PATH"]+"/results/delta_psi/{experiment}",
        julia = config["julia"],
        whippet_bin = config["whippet_bin"]
    shell:
        "{params.julia} {params.whippet_bin}whippet-delta.jl -a {params.grp1} -b {params.grp2} -o {params.out}"
