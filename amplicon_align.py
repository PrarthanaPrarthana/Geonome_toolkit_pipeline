#!/usr/bin/env python
import argparse
import magnumopus

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

    # Getting amplicons
    amplicons1 = magnumopus.ispcr(args.primer_file, args.assembly1, args.max_amplicon_size)
    amplicons2 = magnumopus.ispcr(args.primer_file, args.assembly2, args.max_amplicon_size)
  
    # Align the amplicons
    alignment1, score1 = magnumopus.needleman_wunsch(amplicons1.split("\n")[1], amplicons2.split("\n")[1], args.match, args.mismatch, args.gap)
    
    # reverse complement of amplicon 2
    rev_comp_amplicon2 = reverse_complement(amplicons2.split("\n")[1])
    
    # Align amplicon 1 with the reverse complement of amplicon 2
    alignment2, score2 = magnumopus.needleman_wunsch(amplicons1.split("\n")[1], rev_comp_amplicon2, args.match, args.mismatch, args.gap)

    if score1 >= score2:
        print(alignment1[0])
        print(alignment1[1])
        print(score1)
    else:
        print(alignment2[0])
        print(alignment2[1])
        print(score2)

def reverse_complement(sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[base] for base in sequence[::-1])

if __name__ == "__main__":
    main()

