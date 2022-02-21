containers = {
    "pysam": "docker://quay.io/biocontainers/pysam:0.18.0--py39h20405f9_0",
    "bcftools": "docker://quay.io/biocontainers/bcftools:1.14--hde04aa1_1",
}
default = {}


def get_vcf(wildcards):
    return pep.sample_table.loc[wildcards.sample, "vcf"]


def gather_final_output(wildcards):
    checkpoint_output = checkpoints.vcf_combo.get(**wildcards).output[0]

    return expand(
        "{sample}/fasta/{i}_{ab}_all.fasta",
        sample=wildcards.sample,
        ab=["A", "B"],
        i=glob_wildcards(os.path.join(checkpoint_output, "{i}_A.vcf.gz")).i,
    )
