#!/usr/bin/env python
import subprocess
import pandas as pd
import io
import itertools
import tempfile

def ispcr(primer_file: str, assembly_file: str, max_amplicon_size: int) -> str:
    # Function to extract sequence data from a FASTA file without headers
    def extract_sequence_from_fasta(file):
        sequences = []
        current_sequence = ''
        with open(file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    if current_sequence:
                        sequences.append(current_sequence)
                        current_sequence = ''
                else:
                    current_sequence += line.strip()
        if current_sequence:
            sequences.append(current_sequence)
        return sequences

    # Extract sequences from primer and assembly files
    primer_sequences = extract_sequence_from_fasta(primer_file)
    assembly_sequences = extract_sequence_from_fasta(assembly_file)

    # Run blastn
    result = subprocess.run(['blastn', '-task', 'blastn-short', '-query', primer_file, '-subject', assembly_file,
                            '-outfmt', '6 std qlen'],
                           capture_output=True, text=True)
    df = pd.read_csv(io.StringIO(result.stdout), sep='\t', header=None)

    # Filter blastn results
    filtered = df[(df[3] == df[12]) & (df[2] >= 80)]

    filtered = filtered.astype(str)

    # Find matching primer pairs
    matching_primer_pairs = []

    for hit1, hit2 in itertools.combinations(filtered.values.tolist(), 2):
        qseqid1, sstart1, send1 = hit1[1], int(hit1[8]), int(hit1[9])
        qseqid2, sstart2, send2 = hit2[1], int(hit2[8]), int(hit2[9])

        # Check if primers anneal to the same sequence
        if qseqid1 != qseqid2:
            continue

        # Check if primers anneal close enough
        if abs(sstart1 - sstart2) <= max_amplicon_size:
            if sstart1 < sstart2:
                matching_primer_pairs.append([hit1, hit2])
            else:
                matching_primer_pairs.append([hit2, hit1])

    # Step 4: Create BED file and extract amplicon sequences
    bed_contents = []
    prev_end = None

    for pair in matching_primer_pairs:
        start, end = int(pair[0][9]), int(pair[1][9]) - 1

        if prev_end is not None and start < prev_end:
            # Overlapping pairs, adjust the start position
            start = prev_end + 1

        if start <= end:
            # Append to the BED content
            bed_contents.append(f"{pair[0][1]}\t{start}\t{end}")

        prev_end = end

    # Create a temporary BED file
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as bed_file:
        bed_file.write('\n'.join(bed_contents))

    # Run seqtk to extract amplicon sequences
    seqtk_cmd = f'seqtk subseq {assembly_file} {bed_file.name}'
    amplicon_sequences = subprocess.check_output(seqtk_cmd, shell=True, text=True)

    # Clean up the temporary BED file
    subprocess.run(['rm', bed_file.name])

    return amplicon_sequences

