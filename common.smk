containers = {
    "haploblock-shuffler": "docker://quay.io/redmar_van_den_berg/haploblock-shuffler:0.0.6",
    # bcftools 1.10.2, samtools 1.10, both using htslib 1.10.2
    "bcftools": "docker://quay.io/biocontainers/mulled-v2-03d30cf7bcc23ba5d755e498a98359af8a2cd947:40ff43e422729149fe4c282ed29f2513644857f0-0",
    "pyfasta": "docker://quay.io/biocontainers/pyfasta:0.5.2--py_1",
}
default = {"max_blocks": 11}


def get_vcf(wildcards):
    return pep.sample_table.loc[wildcards.sample, "vcf"]


def gather_final_output(wildcards):
    checkpoint_output = checkpoints.haploblock_shuffler.get(**wildcards).output[0]

    return expand(
        "{sample}/seq/{i}_hap{hap}.seq",
        sample=wildcards.sample,
        hap=[1, 2],
        i=glob_wildcards(os.path.join(checkpoint_output, "out_{i}.vcf.gz")).i,
    )
