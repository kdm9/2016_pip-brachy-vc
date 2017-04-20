SAMPLES = config['samples']
REFERENCE = "data/reference/Bdistachyon_314_v3.0.fa"
VARCALLERS = ['mpileup', 'freebayes']

# bash safe mode
shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

onsuccess:
    shell("mail -s 'Workflow complete' kevin.murray@anu.edu.au pip.wilson@anu.edu.au < {log}")
onerror:
    shell("mail -s 'Workflow error' kevin.murray@anu.edu.au < {log}")


rule all:
    input:
        expand("data/reads/qc/{sample}_il.fastq.gz", sample=SAMPLES),
        expand("data/variants_filtered/{vcf}_strict.vcf.gz", vcf=VARCALLERS),


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
    threads: 4
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
        "data/variants/mpileup.bcf"
    log:
        "data/logs/mpileup.log"
    params:
        ref=REFERENCE,
    threads: 1
    shell:
        "( samtools mpileup"
        "   -d 5000" # Max depth
        "   -t AD" # allele depth tag
        "   -t INFO/AD" # allele depth info
        "   -t DP" # read depth tag
        "   -gu"
        "   -f {params.ref}"
        "   {input}"
        "| bcftools call"
        "   -vm"
        "   -Ob"
        "   -o {output}"
        "   -" # stdin
        " ) 2>{log}"

rule freebayes:
    input:
        "data/merged.bam"
    output:
        "data/variants/freebayes.bcf"
    log:
        "data/logs/freebayes.log"
    threads: 1
    params:
        ref=REFERENCE,
    shell:
        "( freebayes"
        "   -f {params.ref}"
        "   {input}"
        " | bcftools view"
        "   -o {output}"
        "   -O b"
        "   -"  # STDIN
        " ) 2>{log}"

rule vcffilt:
    input:
        "data/variants/{vcf}.bcf"
    output:
        "data/variants_filtered/{vcf}_strict.vcf.gz"
    log:
        "data/logs/vcffilt-{vcf}.log"
    threads: 1
    params:
        ref=REFERENCE,
    shell:
        "( bcftools view"
        " -o {output}"
        " -Oz"
        " -e 'QUAL < 20 || MAF <= 0.02'" # exclude sites
        " -g '^het'" # filter out hets
        " {input}"
        " ) 2>{log}"
