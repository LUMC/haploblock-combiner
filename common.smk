containers = {"pysam": "docker://quay.io/biocontainers/pysam:0.18.0--py39h20405f9_0"}
default = {}


def get_vcf(wildcards):
    return pep.sample_table.loc[wildcards.sample, "vcf"]
