#!/usr/bin/env python3

import argparse
import magnumopus
import re

def format_sequence(sequence):
    # Remove any non-alphabet characters and replace spaces with dashes
    return re.sub(r'[^a-zA-Z]+', '-', sequence)

def remove_sequence_headers(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    sequences.append(sequence)
                sequence = ''
            else:
                sequence += line
        if sequence:
            sequences.append(sequence)
    return sequences

def main():
    # Create a command-line argument parser
    parser = argparse.ArgumentParser(description="Perform in-silico PCR on two assemblies and align the amplicons")

    parser.add_argument('-1', dest='assembly1', required=True, help="Path to the first assembly file")
    parser.add_argument('-2', dest='assembly2', required=True, help="Path to the second assembly file")
    parser.add_argument('-p', dest='primer_file', required=True, help="Path to the primer file")
    parser.add_argument('-m', dest='max_amplicon_size', type=int, required=True, help="Maximum amplicon size for isPCR")
    parser.add_argument('--match', type=int, required=True, help="Match score to use in alignment")
    parser.add_argument('--mismatch', type=int, required=True, help="Mismatch penalty to use in alignment")
    parser.add_argument('--gap', type=int, required=True, help="Gap penalty to use in alignment")

    args = parser.parse_args()

    # Remove headers and format sequences from the assembly and primer files
    sequences1 = remove_sequence_headers(args.assembly1)
    sequences2 = remove_sequence_headers(args.assembly2)
    primer_sequences = remove_sequence_headers(args.primer_file)

    formatted_sequences1 = [format_sequence(seq) for seq in sequences1]
    formatted_sequences2 = [format_sequence(seq) for seq in sequences2]

    # Perform in-silico PCR on the provided sequences
    amplicons1 = magnumopus.ispcr(primer_sequences, formatted_sequences1, args.max_amplicon_size)
    amplicons2 = magnumopus.ispcr(primer_sequences, formatted_sequences2, args.max_amplicon_size)

    # Align the amplicons and get the best alignment
    alignment, score = magnumopus.needleman_wunsch(amplicons1, amplicons2, args.match, args.mismatch, args.gap)

    # Print the formatted alignment and alignment score
    print("Alignment:")
    print(alignment[0])
    print(alignment[1])
    print("Alignment Score:", score)

if __name__ == "__main__":
    main()

