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
        final_files=expand(
            "{sample}/final_files.txt", sample=pep.sample_table["sample_name"]
        ),


rule exclude_homref:
    """ Exclude homref calls, and restrict the output VCF file to the specified
    region. This saves runtime on larg VCF files
    """
    input:
        vcf=get_vcf,
    output:
        vcf="{sample}/{sample}_no_homref.vcf",
    params:
        region=f"--regions {config['region']}" if "region" in config else "",
    log:
        "log/{sample}_exclude_homref.txt",
    container:
        containers["bcftools"]
    shell:
        """
        bcftools view \
            --include 'GT[*]="alt"' \
            {params.region} \
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
    params:
        region=config.get("region", ""),
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

        if [ -z {params.region} ]; then
            cat {input.ref} |
            bcftools consensus \
                --haplotype A \
                {input.vcf} > {output.fasta} 2>> {log}
        else
            samtools faidx {input.ref} {params.region} |
            bcftools consensus \
                --haplotype A \
                {input.vcf} > {output.fasta} 2>> {log}
        fi
        """


rule fasta_to_seq:
    input:
        fasta="{sample}/fasta/{i}_{ab}_all.fasta",
        src=srcdir("scripts/fasta_to_seq.py"),
    output:
        seq="{sample}/seq/{i}_{ab}_all.seq",
    params:
        "--reverse-complement",
    log:
        "log/{sample}/seq_{i}_{ab}_all.txt",
    container:
        containers["pyfasta"]
    shell:
        """
        mkdir -p $(dirname {output.seq})

        python {input.src} \
            --fasta {input.fasta} \
            {params} > {output.seq} 2> {log}
        """


rule gather_final_outputs:
    input:
        gather_final_output,
    output:
        "{sample}/final_files.txt",
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
