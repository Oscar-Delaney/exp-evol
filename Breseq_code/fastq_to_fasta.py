import os
import random
import gzip

def main():
    input_file = os.path.expanduser("~/liv/breseq_output/SD8826_S39/data/SD8826_S39_R1_001.unmatched.fastq")
    output_file = os.path.expanduser("~/liv/breseq_output/SD8826_S39/data/SD8826_S39_R1_001_subset.fasta")

    with open(input_file, 'rt') as fastq_file:
        reads = parse_fastq(fastq_file)

    subset = random.sample(reads, 100)
    write_fasta(subset, output_file)


def parse_fastq(fastq_file):
    reads = []
    for line in fastq_file:
        if line.startswith("@"):
            header = line.strip()
            sequence = fastq_file.readline().strip()
            plus_line = fastq_file.readline()
            quality = fastq_file.readline().strip()
            reads.append((header, sequence))

    return reads


def write_fasta(reads, output_file):
    with open(output_file, 'wt') as fasta_file:
        for header, sequence in reads:
            fasta_file.write(f">{header[1:]}\n{sequence}\n")


if __name__ == "__main__":
    main()

