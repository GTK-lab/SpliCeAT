rule whippet_delta:
    input:
        unpack(delta_input)
    output:
        config["base_path"]+"/results/delta_psi/{experiment}.diff.gz"
    params:
        grp1 = lambda wildcards: grp1(wildcards),
        grp2 = lambda wildcards: grp2(wildcards),
        out = config["base_path"]+"/results/delta_psi/{experiment}"
    shell:
        "~/julia-1.7.2/bin/julia ~/Whippet.jl/bin/whippet-delta.jl -a {params.grp1} -b {params.grp2} -o {params.out}"