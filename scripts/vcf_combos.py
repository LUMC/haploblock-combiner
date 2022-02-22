from pysam import VariantFile
import argparse
import logging
import pprint


# TODO: how to deal with gene deletion and duplication?!


def walk(alleles, stack, results):
    logging.debug(f"Walk input:{alleles}")
    logging.debug(f"Stack:{stack}")
    if not alleles:
        logging.debug("No remaining input")
        logging.debug(f"Storing stack to results:{stack}")
        results.append(stack.copy())
        if stack:
            p = stack.pop()
            logging.debug(f"Popped: '{p}' at no input")
        return

    pos, suffix = alleles[0].split("_")

    logging.debug(f"Put element '{pos}_{suffix}' on stack")
    stack.append(alleles[0])

    if len(alleles) > 1:
        pos_next, suffix_next = alleles[1].split("_")
    else:
        logging.debug("We are at last input")
        logging.debug(f"Storing stack to results: {stack}")
        results.append(stack.copy())
        p = stack.pop()
        logging.debug(f"Popped: '{p}' at last input")
        return

    if pos == pos_next:
        # If next is paired allele, skip over that
        logging.debug(f"Going to skip (pos={pos}) ...")
        walk(alleles[2:], stack, results)
        logging.debug(f"After skip (pos={pos}) ...")

        if stack and stack[-1] == f"{pos}_{suffix}":
            p = stack.pop()
            logging.debug(f"Popped: '{p}' at {pos}")

    logging.debug(f"Going to walk next (pos={pos}) ...")
    walk(alleles[1:], stack, results)

    # Only remove elements from stack that lay beyond
    if stack and stack[-1].split("_")[0] >= pos:
        p = stack.pop()
        logging.debug(f"Popped: '{p}' at end of element '{pos}_{suffix}'")


def wrap(alleles):

    results = []
    walk(alleles, [], results)

    lengths = set()
    for path in results:
        logging.debug(f"path {results.index(path)}")
        logging.debug(f"{len(path)} {path}")
        lengths.add(len(path))

    logging.debug(len(results))

    # Make sure all paths have the same number of phase sets
    assert len(lengths) == 1

    return results


def work(vcf_in, prefix):
    output = {}

    # Detect all PS's
    for rec in vcf_in.fetch():

        assert len(rec.samples) == 1

        sample = rec.samples[0]
        if "PS" in rec.format:
            ps = sample["PS"]
            if ps not in output:
                for allele in ["A", "B"]:
                    ps_allele = f"{ps}_{allele}"
                    output[ps_allele] = []
        # Homozygous
        elif rec.info["AC"][0] == 2:
            output[f"{rec.pos}_HOM"] = []
        else:
            output[f"{rec.pos}_HET"] = []
            output[f"{rec.pos}_X"] = []

    for rec in vcf_in.fetch():
        logging.debug(f"{rec.chrom}:{rec.pos}")
        if rec.info["AC"] not in [(1,), (2,)]:
            msg = f"{rec.chrom}:{rec.pos} has AC={rec.info['AC']}, skipping variant completely"
            logging.warning(msg)
            continue

        assert rec.info["AC"] in [(1,), (2,)], msg

        sample = rec.samples[0]

        # Assume these formats
        assert sample["GT"] in [(0, 1), (1, 0), (1, 1)]

        # Variants in a Phase Set
        if "PS" in rec.format:

            # Assumed to be phased
            assert sample.phased

            # Assume these formats
            assert sample["GT"] in [(0, 1), (1, 0)]

            ps = sample["PS"]

            if sample["GT"] == (1, 0):
                allele = "A"
            elif sample["GT"] == (0, 1):
                allele = "B"

            ps_allele = f"{ps}_{allele}"
            output[ps_allele].append(rec)

        # Heterozygous unphased
        elif rec.info["AC"][0] != 2:
            assert not sample.phased

            # "create" new phase set consisting of only this heterozygous
            # variant and all homozygous variants
            output[f"{rec.pos}_HET"].append(rec)

        # Homozygous
        elif rec.info["AC"][0] == 2:
            assert not sample.phased
            assert "PS" not in rec.format
            assert sample["GT"] in [(1, 1)]
            output[f"{rec.pos}_HOM"].append(rec)

    output_list = []
    for key in output:
        output_list.append(key)

    # Make sure the output list is still sorted
    assert output_list == sorted(output_list, key=lambda x: int(x.split("_")[0]))

    results = wrap(output_list)

    unique_pairs = []
    for path in results:

        logging.debug(path)
        alt_path = []
        for ps in path:
            m = {"A": "B", "B": "A", "HET": "X", "X": "HET", "HOM": "HOM"}
            pos, suffix = ps.split("_")
            alt = f"{pos}_{m[suffix]}"
            alt_path.append(alt)
            logging.debug(f"{ps} {alt}")

        if (path, alt_path) not in unique_pairs and (
            alt_path,
            path,
        ) not in unique_pairs:
            unique_pairs.append((path, alt_path))

    for index, (path, alt_path) in enumerate(unique_pairs):
        vcf1_out = VariantFile(f"{prefix}/{index}_A.vcf.gz", "w", header=vcf_in.header)
        vcf2_out = VariantFile(f"{prefix}/{index}_B.vcf.gz", "w", header=vcf_in.header)

        for ps in path:
            for v in output[ps]:
                vcf1_out.write(v)

        for alt in alt_path:
            for v in output[alt]:
                vcf2_out.write(v)

        vcf1_out.close()
        vcf2_out.close()
    logging.info(f"Created {index} file pairs in {prefix}")


def test():
    assert wrap([]) == [[]]
    assert wrap(["1_A"]) == [["1_A"]]
    assert wrap(["1_A", "2_HOM"]) == [["1_A", "2_HOM"]]
    assert wrap(["1_A", "1_B", "2_HOM"]) == [["1_A", "2_HOM"], ["1_B", "2_HOM"]]


def main(args):
    test()
    vcf_in = VariantFile(args.filename)
    work(vcf_in, args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate pairs.")
    parser.add_argument("filename", help="vcf file")
    parser.add_argument("output", help="output dir file")
    parser.add_argument(
        "--log-level", choices=["DEBUG", "INFO", "WARN", "ERROR"], default="INFO"
    )
    args = parser.parse_args()

    log_levels = {
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARN": logging.WARN,
        "ERROR": logging.ERROR,
    }
    logging.basicConfig(level=log_levels[args.log_level])

    main(args)
