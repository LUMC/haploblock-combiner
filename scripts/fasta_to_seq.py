import argparse
import pyfasta


def fasta2seq(fn, rc):

    fasta = pyfasta.Fasta(fn)

    # Make sure the fasta file contains a single chromosome
    entries = fasta.keys()
    assert len(entries) == 1

    for e in entries:
        # Set the strand we are interested in
        strand = "-" if rc else "+"

        seq = fasta.sequence(
            {"chr": e, "start": 0, "stop": -1, "strand": strand}, one_based=False
        )
        print(seq)


def main():
    parser = argparse.ArgumentParser(
        description="Mangling star alleles for fun and profit."
    )
    parser.add_argument("--fasta", help="Input fasta file to convert to seq",
            required=True)
    parser.add_argument(
        "--reverse-complement",
        default=False,
        action="store_true",
        help="Output the reverse complement of the fasta file",
    )
    args = parser.parse_args()

    fasta2seq(args.fasta, args.reverse_complement)


if __name__ == "__main__":
    main()
