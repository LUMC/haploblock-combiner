containers = {
    "pysam": "docker://quay.io/biocontainers/pysam:0.18.0--py39h20405f9_0",
    # bcftools 1.10.2, samtools 1.10, both using htslib 1.10.2
    "bcftools": "docker://quay.io/biocontainers/mulled-v2-03d30cf7bcc23ba5d755e498a98359af8a2cd947:40ff43e422729149fe4c282ed29f2513644857f0-0",
    "pyfasta": "docker://quay.io/biocontainers/pyfasta:0.5.2--py_1",
}
default = {}


def get_vcf(wildcards):
    return pep.sample_table.loc[wildcards.sample, "vcf"]


def gather_final_output(wildcards):
    checkpoint_output = checkpoints.vcf_combo.get(**wildcards).output[0]

    return expand(
        "{sample}/seq/{i}_{ab}_all.seq",
        sample=wildcards.sample,
        ab=["A", "B"],
        i=glob_wildcards(os.path.join(checkpoint_output, "{i}_A.vcf.gz")).i,
    )
