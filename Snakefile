from os.path import join
configfile: "config.yaml"
ref_genome=config["BOWTIEINDEX"]
thread=config["THREAD"]
trim_length=config["TRIM"]

import pandas as pd
from pathlib import Path
df = pd.read_csv('sample.tsv', sep='\t', index_col='Sample')
ID = df.index


for sample in ID:
  print("Sample " + sample + " will be processed")

rule all:
    input:
        expand("log/atac_{sample}.txt", sample=ID)

rule run_trim:
    input:
        "raw_fastq/{sample}_R1_001.fastq.gz",
        "raw_fastq/{sample}_R2_001.fastq.gz"
    output:
        "raw_fastq/{sample}_R1_001.trim.fastq.gz",
        "raw_fastq/{sample}_R2_001.trim.fastq.gz"
    shell:
        """
         /home/fgao/software/seqtk/seqtk trimfq -L {trim_length} {input[0]} | gzip > {output[0]}
         /home/fgao/software/seqtk/seqtk trimfq -L {trim_length} {input[1]} | gzip > {output[1]}
        """

rule run_bowtie:
    input:
        "raw_fastq/{sample}_R1_001.trim.fastq.gz",
        "raw_fastq/{sample}_R2_001.trim.fastq.gz",
        ref_genome
    output:
        "alignment/{sample}.sam"
    shell:
        """
         bowtie2 -x {input[2]}/genome -1 {input[0]} -2 {input[1]} -S {output} \
        --end-to-end --very-sensitive --no-mixed --no-discordant -q --phred33 -I 10 -X 700 -p {thread}
        """

rule run_bam:
    input:
        "alignment/{sample}.sam"
    output:
        "log/atac_{sample}.txt"
    params:
        sample="{sample}"
    shell:
        """
         samtools view -bS {input} > alignment/{params.sample}.bam
         samtools sort alignment/{params.sample}.bam > alignment/{params.sample}.sort.bam
         samtools index alignment/{params.sample}.sort.bam
         echo "ATACseq read alignment is complete" > {output}
        """
