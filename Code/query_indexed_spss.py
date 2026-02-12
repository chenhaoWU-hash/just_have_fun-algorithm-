import argparse
import tracemalloc
import fmi
from spss import rev_complement
import time


def shared_kmers(sequence: str, my_fmi: fmi.FmIndex, k: int) -> float:
    """Returns the number of shared kmers between a sequence and a fmi."""
    nb_kmers = len(sequence) - k + 1
    nb_shared_kmers = 0
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i : i + k]
        if my_fmi.contains(kmer):
            nb_shared_kmers += 1
        else:
            kmer = rev_complement(kmer)
            if my_fmi.contains(kmer):
                nb_shared_kmers += 1
    return round(nb_shared_kmers / nb_kmers, 3)


def parse_input_fasta_file(fasta_file_name: str, k: int, my_fmi: fmi.FmIndex, output_file_name: str):
    """For each sequence in the fasta file, calculate the number of shared kmers with the fmi."""
    tracemalloc.start()  # Start memory tracking
    sequence_count = 0

    with open(output_file_name, encoding="utf-8", mode="w") as output_stream, open(
        fasta_file_name, encoding="utf-8"
    ) as fasta_file:
        snapshot1 = tracemalloc.take_snapshot()  # Initial snapshot

        while True:
            line = fasta_file.readline()
            if not line:
                break
            if not line[0] == ">":
                print(f"Fasta format error, line {line}")
            comment = line.strip(">").strip()
            line = fasta_file.readline()
            line = line.strip()

            # Measure memory and processing for each sequence
            start_time = time.time()
            ratio = shared_kmers(line, my_fmi, k)
            end_time = time.time()

            sequence_count += 1
            output_stream.write(f"{comment}\t{ratio}\n")

            # Capture memory snapshot every 10 sequences
            if sequence_count % 10 == 0:
                snapshot2 = tracemalloc.take_snapshot()
                stats = snapshot2.compare_to(snapshot1, "lineno")
                total_memory_change = sum(stat.size_diff for stat in stats)
                print(f"[Sequence {sequence_count}] Memory change: {total_memory_change / 1024:.2f} KiB")
                print(f"[Sequence {sequence_count}] Processing time: {end_time - start_time:.2f} seconds")
                snapshot1 = snapshot2  # Update snapshot

    tracemalloc.stop()  # Stop memory tracking


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="QueryIndexedSPSS with optimized memory tracking")
    parser.add_argument("-q", "--query", help="query file name", required=True)
    parser.add_argument("-i", "--index", help="index file name", required=True)
    parser.add_argument("-k", "--kmer", help="kmer size", type=int, required=True)
    parser.add_argument("-o", "--output", help="output file name", required=True)
    args = parser.parse_args()

    tracemalloc.start()  # Start memory tracking for the whole script

    # Load the FM-index
    print("Loading FM-index...")
    snapshot1 = tracemalloc.take_snapshot()
    start_time = time.time()
    fmi_index = fmi.load_fm_index(args.index)
    end_time = time.time()
    snapshot2 = tracemalloc.take_snapshot()
    stats = snapshot2.compare_to(snapshot1, "lineno")
    total_memory_change = sum(stat.size_diff for stat in stats)
    print(f"[FM-index loading] Memory change: {total_memory_change / 1024:.2f} KiB")
    print(f"[FM-index loading] Processing time: {end_time - start_time:.2f} seconds")

    # Parse the input FASTA file
    print("Processing query sequences...")
    parse_input_fasta_file(args.query, args.kmer, fmi_index, args.output)

    tracemalloc.stop()  # Stop memory tracking


if __name__ == "__main__":
    main()
