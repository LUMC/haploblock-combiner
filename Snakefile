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
        dummy_files=expand(
            "{sample}/dummy_file.txt", sample=pep.sample_table["sample_name"]
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
        region=f"--regions {config['region']}",
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


checkpoint haploblock_shuffler:
    input:
        vcf=rules.exclude_homref.output.vcf,
    output:
        folder=directory("{sample}/combinations"),
    params:
        config["max_blocks"],
    log:
        "log/{sample}_haploblock_shuffler.txt",
    container:
        containers["haploblock-shuffler"]
    shell:
        """
        haploblock-shuffler \
            --max-blocks {params} \
            {input.vcf} \
            {output.folder} 2>&1 > {log}

        shopt -s nullglob
        for vcf in {output.folder}/out_*.vcf; do
            bgzip $vcf && tabix -p vcf $vcf.gz
        done 2>&1 >> {log}
        """


rule apply_variants:
    input:
        vcf="{sample}/combinations/out_{i}.vcf.gz",
        ref=config["reference"],
    output:
        hap1="{sample}/fasta/{i}_hap1.fasta",
        hap2="{sample}/fasta/{i}_hap2.fasta",
    params:
        region=config["region"],
    log:
        "log/{sample}_apply_variants_{i}.txt",
    container:
        containers["bcftools"]
    shell:
        """
        # Create the output folder
        mkdir -p $(dirname {output.hap1}) 2> {log}

        samtools faidx {input.ref} {params.region} |
        bcftools consensus --haplotype 1 {input.vcf} > {output.hap1} 2>> {log}

        samtools faidx {input.ref} {params.region} |
        bcftools consensus --haplotype 2 {input.vcf} > {output.hap2} 2>> {log}
        """


rule fasta_to_seq:
    input:
        hap1=rules.apply_variants.output.hap1,
        hap2=rules.apply_variants.output.hap2,
        src=srcdir("scripts/fasta_to_seq.py"),
    output:
        seq1="{sample}/seq/{i}_hap1.seq",
        seq2="{sample}/seq/{i}_hap2.seq",
    params:
        "--reverse-complement",
    log:
        "log/{sample}/fasta_to_seq_{i}.txt",
    container:
        containers["pyfasta"]
    shell:
        """
        mkdir -p $(dirname {output.seq1})

        python {input.src} --fasta {input.hap1} {params} > {output.seq1} 2> {log}
        python {input.src} --fasta {input.hap2} {params} > {output.seq2} 2>> {log}
        """


rule gather_final_outputs:
    """ Dummy rule that gathers all checkpoint output """
    input:
        gather_final_output,
    output:
        "{sample}/dummy_file.txt",
    log:
        "log/{sample}_gather.txt",
    container:
        containers["bcftools"]
    shell:
        """
        touch {output}
        """
