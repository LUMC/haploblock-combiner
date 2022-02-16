containers = {
    "pysam": "docker://quay.io/biocontainers/pysam:0.18.0--py39h20405f9_0",
    "bcftools": "docker://quay.io/biocontainers/bcftools:1.14--hde04aa1_1",
}
default = {}


def get_vcf(wildcards):
    return pep.sample_table.loc[wildcards.sample, "vcf"]
