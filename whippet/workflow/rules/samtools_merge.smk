rule samtools_merge:
    input:
        bams=lambda wildcards: samtools_input(wildcards)
    output:
        bam=config["bam_dir"]+"{experiment}_treatment_merged.bam",
        bai=config["bam_dir"]+"{experiment}_treatment_merged.bam.bai"
    threads:
        16
    shell:
        "samtools merge {input.bams} -@ {threads} -o {output.bam} && samtools index {output.bam}"