containers = {"pysam": "docker://quay.io/biocontainers/pysam:0.18.0--py39h20405f9_0"}
default = {}


def get_vcf(wildcards):
    print(f"Function get_vcf: {wildcards}")
    print(pep.sample_table.loc[wildcards.sample, "vcf"])
    return pep.sample_table.loc[wildcards.sample, "vcf"]


def get_forward(wildcards):
    forward = pep.sample_table.loc[wildcards.sample, "forward"]

    # If a single fastq file is specified, forward will be a string
    if isinstance(forward, str):
        return [forward]
    # If multiple fastq files were specified, forward will be a list
    else:
        return forward


def get_reverse(wildcards):
    reverse = pep.sample_table.loc[wildcards.sample, "reverse"]

    # If a single fastq file is specified, reverse will be a string
    if isinstance(reverse, str):
        return [reverse]
    # If multiple fastq files were specified, reverse will be a list
    else:
        return reverse
