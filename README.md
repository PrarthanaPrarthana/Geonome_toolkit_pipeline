# Geonome_toolkit_pipeline

## Overview
The Genome Toolkit is a comprehensive collection of Python scripts designed to facilitate various tasks in bioinformatics, particularly in genome analysis. The toolkit provides a pipeline to perform in-silico PCR procedure, global sequence alignment using the Needleman-Wunsch algorithm, and alignment of amplicons generated from assembly files.

## Components
### 1. ispcr.py
This module implements the in-silico PCR functionality. It includes a unified function named `ispcr` that takes a primer file, an assembly file, and a maximum amplicon size as inputs, and returns predicted amplicons.

### 2. nw.py
The Needleman-Wunsch algorithm is implemented in this module for global sequence alignment. The `needleman_wunsch` function takes two sequences along with match, mismatch, and gap scores as inputs and returns the aligned sequences along with the alignment score.

### 3. amplicon_align.py
This script combines the functionalities of `ispcr.py` and `nw.py` to perform in-silico PCR on assembly files and align the resulting amplicons. It accepts command-line inputs for two assembly files, a primer file, maximum amplicon size, match, mismatch, and gap scores, and prints the best alignment along with the alignment score to the terminal.

## Usage
To use the Genome Toolkit, follow these steps:
1. Clone the repository to your local machine.
2. Navigate to the project directory.
3. Run the desired script with appropriate command-line arguments.

## Example
```bash
python amplicon_align.py -1 assembly1.fna -2 assembly2.fna -p primers.fna -m 2000 --match 1 --mismatch -1 --gap -1

./1.py
./2.py 

