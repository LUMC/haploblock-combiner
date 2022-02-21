include: "common.smk"


pepfile: config["pepfile"]


# Apply the settings from the pepfile, overwriting the default ones
default.update(pep.config.get("haploblock-combiner", dict()))

# Apply the options specified to snakemake, overwriting the default settings
# and the settings from the PEP file
default.update(config)

# Set the updated dict as the configuration for the pipeline
config = default


rule all:
    input:
        fasta=expand("{sample}/fasta_files.txt", sample=pep.sample_table["sample_name"]),


rule exclude_homref:
    input:
        vcf=get_vcf,
    output:
        vcf="{sample}/{sample}_no_homref.vcf",
    log:
        "log/{sample}_exclude_homref.txt",
    container:
        containers["bcftools"]
    shell:
        """
        bcftools view \
            --include 'GT[*]="alt"' \
            {input.vcf} > {output.vcf} 2> {log}
        """


checkpoint vcf_combo:
    input:
        vcf=rules.exclude_homref.output.vcf,
        src=srcdir("scripts/vcf_combos.py"),
    output:
        folder=directory("{sample}/combinations"),
    log:
        "log/{sample}_vcf_combo.txt",
    container:
        containers["pysam"]
    shell:
        """
        mkdir -p {output.folder}
        python {input.src} \
            {input.vcf} {output.folder} 2>&1 > {log}
        """


rule apply_variants:
    input:
        vcf="{sample}/combinations/{i}_{ab}.vcf.gz",
        ref=config["reference"],
    output:
        fasta="{sample}/fasta/{i}_{ab}_all.fasta",
        tbi="{sample}/combinations/{i}_{ab}.vcf.gz.tbi",
    log:
        "log/{sample}_fasta_{i}_{ab}.txt",
    container:
        containers["bcftools"]
    shell:
        """
        # Create the output folder
        mkdir -p $(dirname {output.fasta}) 2> {log}

        # Index the vcf input file
        tabix --force --preset vcf {input.vcf}

        cat {input.ref} | \
        bcftools consensus \
            --haplotype A \
            {input.vcf} > {output.fasta} 2>> {log}
        """


rule gather_fasta:
    input:
        gather_final_output,
    output:
        "{sample}/fasta_files.txt",
    log:
        "log/{sample}_gather.txt",
    container:
        containers["bcftools"]
    shell:
        """
        rm -f {output}
        for file in {input}; do
            echo $file >> {output}
        done
        """
