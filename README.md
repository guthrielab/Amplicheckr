# Amplichekr: a qPCR Primer Checking Tool

## Introduction
Amplichekr is a python based primer checking tool for influenza and hantavirus genomes with the capability for sequences and primers of any origin. 

It achieves this through local alignment of primers and genomes with the Smith-Waterson algorithm followed by scoring using a set criteria.

## Installation and Dependencies
Download the repository:
```
git clone https://github.com/guthrielab/Amplicheckr
```
Amplichekr uses numerous python packages detailed in [amplichekr.toml](https://github.com/guthrielab/Amplicheckr/blob/main/amplichekr.toml) that are not installed with python. 
This includes:

 - numpy
 - pandas
 - biopython
 - jinja2

Ensure the packages are installed by running 
```
pip list
```
in your command terminal.

Amplichekr should be functional on all systems with python installed.

## Workflow and Scoring
The general workflow is as follows:
1. Parse primer csv and genome fasta
2. Create index of kmers for all genomes

For each primer set:

3. Each primer sequence is split into identical length kmers and used as a query in the index (reverse primers are transformed to reverse complement)
4. Once a match is found in the index, a window is set at the matched genome and a matrix is built between the primer and genome window
5. If maximum score meets threshold (no more than 4 mismatches) perform traceback
6. Pass alignment info to either align output or primer scoring and filtering
7. Output findings to stdout or write in html report

### Initial alignment scoring
Alignment matches are scored with the following matrix, using IUPAC nucleotide names for degenerate bases and gap score of -6.
```
    A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N
A   2  -3  -3  -3  -3   2   2  -3  -3   2  -3   2   2   2   2
T  -3   2  -3  -3  -3   2  -3   2   2  -3   2  -3   2   2   2
G  -3  -3   2  -3   2  -3   2  -3   2  -3   2   2  -3   2   2
C  -3  -3  -3   2   2  -3  -3   2  -3   2   2   2   2  -3   2
S  -3  -3   2   2   2  -3   2   2   2   2   2   2   2   2   2
W   2   2  -3  -3  -3   2   2   2   2   2   2   2   2   2   2
R   2  -3   2  -3   2   2   2  -3   2   2   2   2   2   2   2
Y  -3   2  -3   2   2   2  -3   2   2   2   2   2   2   2   2
K  -3   2   2  -3   2   2   2   2   2  -3   2   2   2   2   2
M   2  -3  -3   2   2   2   2   2  -3   2   2   2   2   2   2
B  -3   2   2   2   2   2   2   2   2   2   2   2   2   2   2
V   2  -3   2   2   2   2   2   2   2   2   2   2   2   2   2
H   2   2  -3   2   2   2   2   2   2   2   2   2   2   2   2
D   2   2   2  -3   2   2   2   2   2   2   2   2   2   2   2
N   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2
```
### Primer mismatch scoring
Scoring is conducted similarly to GISAID's primerchecker tool and follows guidelines from [Deustch et. al 2024](https://academic.oup.com/bioinformatics/article/40/11/btae657/7876261) ([Viral Primer](https://viralprimer.elte.hu/)).

These guidelines categorize primers into three grades:
#### 1. High Risk:
- At least 1 mismatch in critical regions (Forward and reverse: 5 bases from 3' end; Probe: central area 5 bases in from both sides)
- and/or 3 mismatches outside critical regions or 2 worst-case mismatches (purine-purine, ex. A-A, A-G, R-R pairings)
#### 2. Moderate Risk:
- No mismatch in critical regions
- At least 1 mismatch of moderate and above severity (purine-purine, pyrmidine-pyrmidine) outside critical regions
#### 3. No Risk:
- No mismatch in critical regions
- Only 1 mismatch of low severity (purine-pyrmidine, ex. A-C) or perfect match

### Filtering and warnings
A set is considered "failing" for a genome when at least 1 primer fails to align to the genome, or the primers align in the incorrect order (correct order: forward < probe < reverse).

A warning message is included in the output detailing which primer in which set caused the failure.

However, for cases with multiple alignments to the genome, the duplicate alignment is dropped and the set continues as normal, but a warning message is included detailing which primer has multiple alignments on which genome.

## Input 
### Required parameters
| Parameter | Description |
|------------|-------------|
| `Primers`  | Name of primer file, including absolute/relative path if not in parent folder |
| `Genomes` | File name of single or directory of multiple genome fasta sequences |
### Additional parameters
| Parameter | Description |
|------------|-------------|
| `-f`, `--files` | Use `Genomes` as a directory input and common suffix to search (ex. fna) |
| `-k`, `--kmer` | Specify kmer length for searching and indexing (default k=7) |
| `-v`, `--virus` | Name of target virus species (influenza/i or hantavirus/h) |
| `-a`, `--align` | Perform alignment of input only without primer scoring. Outputs all alignments above 4 mismatches |
| `-H`, `--html` | Activate html report generator, reroutes from stdout |
| `-o`, `--output` | To be used with --html, directory for output file |
### Primers
Primers should be in file with commas as delimiter (csv). Data must include 5 columns, named "name", "target", "sequence", "set", and "type" (as F/R/P). These columns can appear in any order.

Each primer set should have exactly one forward primer, one reverse primer, and one probe. If a set has multiple probes, duplicate the forward and reverse primer and create a new set with each probe.
### Genomes
Any file in FASTA format is able to be parsed. However, Genbank/NCBI formats work better than GISAID for recording metadata, as the name of each sequence will be just the accession for genbank/NCBI but will be the entire ">" line for GISAID sequences. 

The program will additionally record the segment type for influenza ("HA", "MP", "NA) and hantavirus ("segment S", "segment M", "segment L") and record the year of the entry if present. If not found, the entry is replaced by "unknown" for the segment and/or year.

## Output
In contrast with other primer checking tools, this tool groups matches by the matching sequence rather by each matching genome. 

For example, if a primer matches to the sequence `GCAGCTGTGTCTACATTGGAGAC`, instead of outputting the matching genomes one by one, the tool groups by which genomes contain this sequence and instead provides the name and prevalence in the total uploaded sets. If, for example, 5 genomes had a single mismatch at the 3' terminal region while 10 genomes were perfect matches, the tool will output two entries, with "high risk" at 5 (33%) and "no risk" at 10 (67%).

The prevalence count accounts for prevalence of the segment matched, but if the segment name could not be parsed or genomes from an origin other than influenza or hantavirus are used it will instead count over total entries.

## Test Data
The sequences under [testdata_hantavirus](https://github.com/guthrielab/Amplicheckr/blob/main/testdata_hantavirus/info.md) are the 27 most recent andes hantavirus genomes on NCBI. 


### Standard grading
```
python3 main.py -v hantavirus .\testdata_hantavirus\tprimerhanta.csv .\testdata_hantavirus\tgenome.fna
```
Alternate, using -f
```
python3 main.py -f fna -v hantavirus .\testdata_hantavirus\tprimerhanta.csv .\testdata_hantavirus
```
### HTML report mode
```
python3 main.py -v hantavirus .\testdata_hantavirus\tprimerhanta.csv .\testdata_hantavirus\tgenome.fna -H
```
Find this report by opening "dashboard.html" in the pwd or specified output location
### Align only mode 
```
python3 main.py -v hantavirus .\testdata_hantavirus\tprimerhanta.csv .\testdata_hantavirus\tgenome.fna -a
```
(Warning: long output)


#### Note: if on linux/macOS systems, path directories must be written with the forward slash (" / ") instead of backslash (" \\ "). 



