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
        samples=expand("{sample}/filelist.txt", sample=pep.sample_table["sample_name"]),


# TODO: Fix bcftools filter to use --include 'GT[*]="alt"'
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
            --genotype ^hom \
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
        gather_vcf_combo,
    output:
        "{sample}/filelist.txt",
    log:
        "log/{sample}_apply_variants.txt",
    container:
        containers["bcftools"]
    shell:
        """
        for file in {input}; do
            echo $file >> {output}
        done
        """
