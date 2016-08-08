SAMPLES = config['samples']
REFERENCE = "data/reference/Egrandis_v2.fasta"

# bash safe mode
shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")



rule all:
    input:
        expand("data/reads/qc/{sample}_il.fastq.gz", sample=SAMPLES)


rule qcreads:
    input:
        r1="data/reads/raw/{sample}_R1.fastq.gz",
        r2="data/reads/raw/{sample}_R2.fastq.gz",
    output:
        r="data/reads/qc/{sample}_il.fastq.gz",
        y="data/qcrep/{sample}.yml",
    log:
        "data/logs/qc/{sample}.log"
    threads: 1
    shell:
        "( pairs join {input} "
        "| trimit"
        "   -b"
        "   -Q"
        "   -q 20"
        "   -y {output.y}"
        "   -" # stdin
        "| gzip > {output.r}"
        " )2>{log}"

rule align:
    input:
        "data/reads/qc/{sample}_il.fastq.gz"
    output:
        "data/alignments/{sample}.bam"
    log:
        "data/logs/align/{sample}.log"
    threads: 8
    params:
        ref=REFERENCE,
    shell:
        "( bwa mem"
        "   -p" # detect pairs in IL file
        "   -R '@RG\\tID:{wildcards.sample}\\tPL:ILLUMINA\\tSM:{wildcards.sample}'"
        "   -t {threads}"
        "   {params.ref}"
        "   {input}"
        " | samtools view -Suh -"
        " | samtools sort"
        "    -o {output}"
        "    -" # stdin
        " && samtools index {output}"
        " ) 2>{log}"

rule mergebam:
    input:
        expand("data/alignments/{sample}.bam", sample=SAMPLES)
    output:
        "data/merged.bam"
    log:
        "data/logs/merge.log"
    threads: 8
    shell:
        "( samtools merge"
        "   -@ {threads}"
        "   {output}"
        "   {input}"
        " && samtools index {output}"
        " ) 2>{log}"

rule mpileup:
    input:
        "data/merged.bam"
    output:
        "data/variants/samplewise-mpileup-raw.vcf"
    log:
        "data/logs/mpileup.log"
    threads: 1
    shell:
        "( samtools mpileup"
        "   -C 50" # Cap MQ
        "   -d 5000" # Max depth
        "   -u" # uncompressed bcf
        "   -E" # redo baq calc
        "   -t AD" # allele depth tag
        "   -t INFO/AD" # allele depth info
        "   -t DP" # read depth tag
        "   {input}"
        "| bcftools call -vm -"
        "   >{output}"
        " ) 2>{log}"

rule mpileup_all:
    input:
        "data/merged.bam"
    output:
        "data/variants/concat-mpileup-raw.bcf"
    log:
        "data/logs/mpileup-all.log"
    threads: 1
    params:
        ref=REFERENCE,
    shell:
        "( samtools mpileup"
        "   -R" # ignore rg, treat all as one sample
        "   -C 50" # Cap MQ
        "   -f {params.ref}"
        "   -d 5000" # Max depth
        "   -u" # uncompressed bcf
        "   -E" # redo baq calc
        "   -t AD" # allele depth tag
        "   -t INFO/AD" # allele depth info
        "   -t DP" # read depth tag
        "   {input}"
        "| bcftools call"
        "    -mv"
        "    -Ob"
        "    -o {output}"
        "    -" # stdin
        " ) 2>{log}"

rule freebayes:
    input:
        "data/merged.bam"
    output:
        "data/variants/freebayes-raw.vcf"
    log:
        "data/logs/freebayes.log"
    threads: 1
    params:
        ref=REFERENCE,
    shell:
        "( freebayes"
        "   -f {params.ref}"
        "   {input}"
        "   >{output}"
        " ) 2>{log}"

rule freebayes_par:
    input:
        "data/merged.bam"
    output:
        "data/variants/par-freebayes-raw.vcf"
    log:
        "data/logs/freebayes.log"
    threads: 12
    params:
        ref=REFERENCE,
    shell:
        "( freebayes-parallel"
        "   {params.ref}.freebayes_regions"
        "   {threads}"
        "   -f {params.ref}"
        "   {input}"
        "   >{output}"
        " ) 2>{log}"

rule vcffilt:
    input:
        "data/variants/{vcf}.bcf"
    output:
        "data/variants_filtered/{vcf}.vcf.gz"
    log:
        "data/logs/vcffilt-{vcf}.log"
    threads: 1
    params:
        ref=REFERENCE,
    shell:
        "( bcftools filter"
        " -o {output}"
        " -Oz"
        " -e 'QUAL < 40'"
        " {input}"
        " ) 2>{log}"

